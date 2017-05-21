"""Script to download or update SIMBAD-related databases"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "17 May 2017"
__version__ = "0.2"

import argparse
import base64
import datetime
import glob
import numpy as np
import os
import shutil
import sys
import tarfile
import tempfile
import time
import urllib2
import zlib

import cctbx.crystal
import simbad.constants
import simbad.command_line
import simbad.exit
import simbad.rotsearch.amore_search
import simbad.util.simbad_util
import simbad.util.workers_util

logger = None

# The space groups in the list below cannot be recognized by CCTBX, so we convert them
# to similar ones understandle by the library
CCTBX_ERROR_SG = {
    'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422', 'I21': 'I2', 'I1211': 'I2', 
    'P21212A': 'P212121', 'R3': 'R3:R', 'C4212': 'P422',
}

SYS_PLATFORM = sys.platform
CUSTOM_PLATFORM = "linux" if SYS_PLATFORM in ["linux", "linux2"] \
                   else "mac" if SYS_PLATFORM in ["mac"] \
                   else "windows"


def download_morda():
    """Download the MoRDa database

    Returns
    -------
    str
       The path to the downloaded MoRDa database

    """
    logger.info("Downloading MoRDa database from CCP4")
    url = "http://www.ccp4.ac.uk/morda/MoRDa_DB.tar.gz"
    tmp_db = os.path.basename(url)
    query = urllib2.urlopen(url)

    # Chunk size writes data as it's read
    # http://stackoverflow.com/a/34831866/3046533
    chunk_size = 1 << 20
    with open(tmp_db, "wb") as f_out:
        while True:
            chunk = query.read(chunk_size)
            if not chunk:
                break
            f_out.write(chunk)
    # Extract relevant files from the tarball
    with tarfile.open(tmp_db) as tar:
        members = [
            tarinfo for tarinfo in tar.getmembers()
            if tarinfo.path.startswith("MoRDa_DB/home/ca_DOM")
            or tarinfo.path.startswith("MoRDa_DB/home/ca_DB")
            or tarinfo.path.startswith("MoRDa_DB/pdb_DB_gz")
            or tarinfo.path.startswith("MoRDa_DB/" + "bin_" + CUSTOM_PLATFORM)
            or tarinfo.path.startswith("MoRDa_DB/list/domain_list.dat")
        ]
        tar.extractall(members=members)
    # Remove the database file
    os.remove(tmp_db)
    return os.path.abspath("MoRDa_DB")


def create_lattice_db(database):
    """Create a lattice search database

    Parameters
    ----------
    database : str
       The path to the database file

    """
    logger.info('Querying the RCSB Protein DataBank')

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=lengthOfUnitCellLatticeA,'\
          + 'lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,unitCellAngleAlpha,unitCellAngleBeta,' \
            'unitCellAngleGamma,spaceGroup,experimentalTechnique&service=wsfile&format=csv'

    crystal_data, error_count = [], 0
    for line in urllib2.urlopen(url):
        if line.startswith('structureId'):
            continue
        pdb_code, rest = line.replace('"', "").split(',', 1)
        unit_cell, space_group, exp_tech = rest.rsplit(',', 2)
        # Ignore non-xtal structures
        if exp_tech.strip().upper() != "X-RAY DIFFRACTION":
            continue
        # Some entries do not have stored unit cell parameters
        try:
            unit_cell = map(float, unit_cell.split(','))
        except ValueError as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            error_count += 1
            continue
        space_group = space_group.replace(" ", "").strip()
        space_group = CCTBX_ERROR_SG.get(space_group, space_group)   
        try:
            symmetry = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group)
        except Exception as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            error_count += 1
            continue
        crystal_data.append((pdb_code, symmetry))
    logger.debug('Error processing %d pdb entries', error_count)

    logger.info('Calculating the Niggli cells')
    niggli_data = np.zeros((len(crystal_data), 10))
    for i, xtal_data in enumerate(crystal_data):
        niggli_data[i][:4] = np.fromstring(xtal_data[0], dtype='uint8').astype(np.float64)
        niggli_data[i][4:] = np.asarray(xtal_data[1].niggli_cell().unit_cell().parameters())
    logger.info("Total Niggli cells loaded: %d", len(crystal_data))

    if not database.endswith('.npz'):
        database += ".npz"
    logger.info('Storing database in file: %s', database)
    np.savez_compressed(database, niggli_data)


def create_morda_db(database, nproc=2, submit_cluster=False, submit_qtype=None, 
                    submit_queue=False, submit_array=None, submit_max_array=None):
    """Create the MoRDa search database

    Parameters
    ----------
    database : str
       The path to the database folder
    nproc : int, optional
       The number of processors [default: 2]
    submit_cluster : bool
       Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system [default: False]
    submit_qtype : str
       The cluster submission queue type - currently support SGE and LSF
    submit_queue : str
       The queue to submit to on the cluster
    submit_array : st
       Submit SGE jobs as array jobs
    submit_max_array : str
       The maximum number of jobs to run concurrently with SGE array job submission
    
    Raises
    ------
    RuntimeError
       Windows is currently not supported

    """
    if CUSTOM_PLATFORM == "windows":
        msg = "Windows is currently not supported"
        raise RuntimeError(msg)

    # Download the MoRDa database
    os.environ['MRD_DB'] = download_morda()

    # Find all relevant dat files in the MoRDa database and check which are new
    morda_dat_path = os.path.join('MoRDa_DB', 'home', 'ca_DOM', '*.dat')
    simbad_dat_path = os.path.join(database, '**', '*.dat')
    morda_dat_files = set([os.path.basename(f) for f in glob.glob(morda_dat_path)])
    simbad_dat_files = set([os.path.basename(f) for f in glob.glob(simbad_dat_path)])
    dat_files = list(morda_dat_files - simbad_dat_files)

    # Check if we even have a job
    if len(dat_files) < 1:
        logger.info('SIMBAD database up-to-date')
        shutil.rmtree(os.environ['MRD_DB'])
        return
    else:
        logger.info("%d new entries were found in the MoRDa database, update SIMBAD database", len(dat_files))

    # Get the "get_model" script to extract the xyz coordinates
    exe = os.path.join(os.environ['MRD_DB'], "bin_" + CUSTOM_PLATFORM, "get_model")
    
    # Creating temporary output directory
    ccp4_scr = os.environ["CCP4_SCR"]
    os.environ["CCP4_SCR"] = tempfile.mkdtemp(dir=os.getcwd())
    logger.debug("$CCP4_SCR environment variable changed to %s", os.environ["CCP4_SCR"])

    # Create the database files
    what_to_do = []
    for f in dat_files:
        code = os.path.basename(f).rsplit('.', 1)[0]
        final_file = os.path.join(database, code[1:3], code + ".dat")
        # We need a temporary directory within because "get_model" uses non-unique file names
        tmp_dir = tempfile.mkdtemp(dir=os.environ["CCP4_SCR"])
        get_model_output = os.path.join(tmp_dir, code + ".pdb")

        # Prepare script for multiple submissions
        script = simbad.util.simbad_util.tmp_file_name(delete=False, directory=tmp_dir,
                                                       suffix=simbad.util.simbad_util.SCRIPT_EXT)
        log = script.rsplit('.', 1)[0] + '.log'
        with open(script, 'w') as f_out:
            f_out.write(simbad.util.simbad_util.SCRIPT_HEADER + os.linesep)
            f_out.write("export MRD_DB=" + os.environ['MRD_DB'] + os.linesep)
            f_out.write(" ".join([exe, "-c", code, "-m", "d"]) + os.linesep)
        os.chmod(script, 0o777)
        what_to_do += [(script, log, tmp_dir, (get_model_output, final_file))]

    # Run the scripts
    scripts, logs, tmps, files = zip(*what_to_do)
    simbad.util.workers_util.run_scripts(
        job_scripts=scripts,
        job_name='morda_db', chdir=True, nproc=nproc,
        submit_cluster=submit_cluster, submit_qtype=submit_qtype,
        submit_queue=submit_queue, submit_array=submit_array,
        submit_max_array=submit_max_array,
    )

    # Create PDB-like database subdirectories
    sub_dir_names = set([os.path.basename(f).rsplit('.', 1)[0][1:3] for f in dat_files])
    for sub_dir_name in sub_dir_names:
        sub_dir = os.path.join(database, sub_dir_name)
        if os.path.isdir(sub_dir):
            continue
        os.makedirs(sub_dir)

    # Move created files to database
    for output, final in files:
        with open(final, 'wb') as f_out:
            compr = zlib.compress(open(output, 'r').read())
            f_out.write(base64.b64encode(compr))

    # Remove temporary files
    shutil.rmtree(os.environ['MRD_DB'])
    shutil.rmtree(os.environ["CCP4_SCR"])
    os.environ["CCP4_SCR"] = ccp4_scr

    # Leave a timestamp
    leave_timestamp('simbad_morda.txt')


def create_sphere_db(database, shres=3, nproc=2, submit_cluster=False, submit_qtype=None, 
                     submit_queue=False, submit_array=None, submit_max_array=None, chunk_size=5000):
    """Create the spherical harmonics search database

    Parameters
    ----------
    database : str
       The path to the database folder
    shres : int, optional
       Spherical harmonic resolution [default 3.0]
    nproc : int, optional
       The number of processors [default: 2]
    submit_cluster : bool, optional
       Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system [default: False]
    submit_qtype : str, optional
       The cluster submission queue type - currently support SGE and LSF
    submit_queue : str, optional
       The queue to submit to on the cluster
    submit_array : str, optional
       Submit SGE jobs as array jobs
    submit_max_array : str, optional
       The maximum number of jobs to run concurrently with SGE array job submission
    chunk_size : int, optional
       The number of jobs to submit at the same time
    
    Raises
    ------
    RuntimeError
       Windows is currently not supported
    ValueError
       The provided SIMBAD-MoRDa database does not seem to be in the correct format

    """
    # Create and/or update the SIMBAD ".dat" database - this is essential to make sure we are up-to-date
    create_morda_db(database, nproc=nproc, submit_cluster=submit_cluster, submit_qtype=submit_qtype,
                    submit_queue=submit_queue, submit_array=submit_array, submit_max_array=submit_max_array)

    # Find all relevant files in the SIMBADa database and check which are new
    simbad_dat_files = glob.glob(os.path.join(database, '**', '*.dat'))
    simbad_hkl_files = glob.glob(os.path.join(database, '**', '*.hkl.tar.gz'))
    simbad_tab_files = glob.glob(os.path.join(database, '**', '*.tab.tar.gz'))
    simbad_clmn_files = glob.glob(os.path.join(database, '**', '*.clmn.tar.gz'))
    simbad_dat_codes = set([os.path.basename(f).replace('.dat', '') for f in simbad_dat_files])
    simbad_hkl_codes = set([os.path.basename(f).replace('_search.hkl.tar.gz', '') for f in simbad_hkl_files]) 
    simbad_tab_codes = set([os.path.basename(f).replace('_search.clmn.tar.gz', '') for f in simbad_tab_files]) 
    simbad_clmn_codes = set([os.path.basename(f).replace('_search-sfs.tab.tar.gz', '') for f in simbad_clmn_files]) 
    # Find the difference between all and the intersection of the other three file types
    dat_codes = list(simbad_dat_codes - (simbad_hkl_codes & simbad_tab_codes & simbad_clmn_codes))
    dat_files = [os.path.abspath(f) for f in simbad_dat_files if os.path.basename(f).replace('.dat', '') in dat_codes]

    # Check if we even have a job
    if len(dat_files) < 1:
        logger.info("SIMBAD database up-to-date")
        return
    else:
        logger.info("%d new entries were found in the MoRDa database, update SIMBAD database", len(dat_files))

    # Creating temporary output directory
    ccp4_scr = os.environ["CCP4_SCR"]
    os.environ["CCP4_SCR"] = tempfile.mkdtemp(dir=os.getcwd())
    logger.debug("$CCP4_SCR environment variable changed to %s", os.environ["CCP4_SCR"])

    # Which AMORE executable we use
    amore_exe = os.path.join(os.environ["CCP4"], "bin", "amoreCCB2.exe")

    # Submit in chunks, so we don't take too much disk space
    # and can terminate without loosing the processed data
    for i in range(0, len(dat_files), chunk_size):
        # Take a chunk
        chunk_dat_files = dat_files[i:i+chunk_size]

        # ============================
        # First round of tabfun
        everything_1 = []
        for xyzin1_dat in chunk_dat_files:
            xyzin1 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.pdb')
            xyzout1 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.pdb')
            table1 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.car')
            script = simbad.util.simbad_util.tmp_file_name(delete=False, directory=os.environ["CCP4_SCR"],
                                                           suffix=simbad.util.simbad_util.SCRIPT_EXT)
            # Convert the dat file to pdb
            with open(xyzin1_dat, 'rb') as f_in, open(xyzin1, 'w') as f_out:
                content = zlib.decompress(base64.b64decode(f_in.read()))
                f_out.write(content)
            # Construct run scripts and log files
            cmd, stdin = simbad.rotsearch.amore_search.AmoreRotationSearch.tabfun(amore_exe, xyzin1, xyzout1, table1)
            log = script.rsplit('.', 1)[0] + '.log'
            with open(script, 'w') as f_out:
                f_out.write(simbad.util.simbad_util.SCRIPT_HEADER + os.linesep)
                f_out.write(" ".join(map(str, cmd)) + " << eof" + os.linesep)
                f_out.write(stdin + os.linesep + "eof" + os.linesep)
            os.chmod(script, 0o777)
            everything_1 += [(script, log, table1, xyzout1)]

        scripts, logs, table1s, _ = zip(*everything_1)
        simbad.util.workers_util.run_scripts(
            job_scripts=scripts,
            job_name='sphere_db_1_{0}'.format(i), chdir=True, nproc=nproc,
            submit_cluster=submit_cluster, submit_qtype=submit_qtype,
            submit_queue=submit_queue, submit_array=submit_array,
            submit_max_array=submit_max_array,
        )
        # Remove some files to clear disk space
        amore_tmps = glob.glob(os.path.join(os.environ["CCP4_SCR"], 'amoreCCB2_*'))
        for f in list(scripts) + list(logs) + list(table1s) + list(amore_tmps):
            if os.path.isfile(f):
                os.remove(f)

        # Get output for next step
        _, _, _, xyzout1s = zip(*everything_1)

        # ============================
        # Second round of tabfun
        everything_2 = []
        for xyzout1 in xyzout1s:
            # Some files fail in step 1, skip those from now on
            if not os.path.isfile(xyzout1):
                continue
            # Get some data and do step 2 for rest
            x, y, z, intrad = simbad.rotsearch.amore_search.AmoreRotationSearch.calculate_integration_box(xyzout1)
            xyzout2 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.pdb')
            table2 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.car')
            script = simbad.util.simbad_util.tmp_file_name(delete=False, directory=os.environ["CCP4_SCR"],
                                                           suffix=simbad.util.simbad_util.SCRIPT_EXT)
            cmd, stdin = simbad.rotsearch.amore_search.AmoreRotationSearch.tabfun(
                amore_exe, xyzout1, xyzout2, table2, x=x, y=y, z=z
            )
            log = script.rsplit('.', 1)[0] + '.log'
            with open(script, 'w') as f_out:
                f_out.write(simbad.util.simbad_util.SCRIPT_HEADER + os.linesep)
                f_out.write(" ".join(map(str, cmd)) + " << eof" + os.linesep)
                f_out.write(stdin + os.linesep + "eof" + os.linesep)
            os.chmod(script, 0o777)
            everything_2 += [(script, log, xyzout2, table2, intrad)]

        scripts, logs, xyzout2s, _, _ = zip(*everything_2)
        simbad.util.workers_util.run_scripts(
            job_scripts=scripts,
            job_name='sphere_db_2_{0}'.format(i), chdir=True, nproc=nproc,
            submit_cluster=submit_cluster, submit_qtype=submit_qtype,
            submit_queue=submit_queue, submit_array=submit_array,
            submit_max_array=submit_max_array,
        )
        # Remove some files to clear disk space
        amore_tmps = glob.glob(os.path.join(os.environ["CCP4_SCR"], 'amoreCCB2_*'))
        for f in list(scripts) + list(logs) + list(xyzout1s) + list(xyzout2s) + list(amore_tmps):
            if os.path.isfile(f):
                os.remove(f)

        # Get output for next step
        _, _, _, table2s, intrads = zip(*everything_2)

        # ============================
        # First round of rotfun
        everything_3 = []
        for table2, intrad in zip(table2s, intrads):
            # Skip those that failed in step 2
            if not os.path.isfile(table2):
                continue
            # Do processing on the rest
            hklpck1 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.car')
            clmn1 = simbad.util.simbad_util.tmp_file_name(directory=os.environ["CCP4_SCR"], suffix='.cof')
            script = simbad.util.simbad_util.tmp_file_name(delete=False, directory=os.environ["CCP4_SCR"],
                                                           suffix=simbad.util.simbad_util.SCRIPT_EXT)
            cmd, stdin = simbad.rotsearch.amore_search.AmoreRotationSearch.rotfun(
                amore_exe, table2, hklpck1, clmn1, shres, intrad
            )
            log = script.rsplit('.', 1)[0] + '.log'
            with open(script, 'w') as f_out:
                f_out.write(simbad.util.simbad_util.SCRIPT_HEADER + os.linesep)
                f_out.write(" ".join(map(str, cmd)) + " << eof" + os.linesep)
                f_out.write(stdin + os.linesep + "eof" + os.linesep)
            os.chmod(script, 0o777)
            everything_3 += [(script, log, hklpck1, clmn1)]

        scripts, logs, _, _ = zip(*everything_3)
        simbad.util.workers_util.run_scripts(
            job_scripts=scripts,
            job_name='sphere_db_3_{0}'.format(i), chdir=True, nproc=nproc,
            submit_cluster=submit_cluster, submit_qtype=submit_qtype,
            submit_queue=submit_queue, submit_array=submit_array,
            submit_max_array=submit_max_array,
        )
        # Remove some files to clear disk space
        amore_tmps = glob.glob(os.path.join(os.environ["CCP4_SCR"], 'amoreCCB2_*'))
        for f in list(scripts) + list(logs) + list(amore_tmps):
            if os.path.isfile(f):
                os.remove(f)

        # Get output for next step
        _, _, hklpck1s, clmn1s = zip(*everything_3)

        # ============================
        # Update the database

        # Create PDB-like database subdirectories
        sub_dir_names = set([os.path.basename(f).rsplit('.', 1)[0][1:3] for f in chunk_dat_files])
        for sub_dir_name in sub_dir_names:
            sub_dir = os.path.join(database, sub_dir_name)
            if os.path.isdir(sub_dir):
                continue
            os.makedirs(sub_dir)

        # Move files
        for dat, table2, hklpck1, clmn1 in zip(chunk_dat_files, table2s, hklpck1s, clmn1s):
            name = os.path.basename(dat).rsplit('.', 1)[0]
            file_combos = [
                (hklpck1, name + "_search.hkl"),
                (clmn1, name + "_search.clmn"),
                (table2, name + "_search-sfs.tab"),
            ]
            for f, new_f in file_combos:
                tarball = new_f + ".tar.gz"
                with tarfile.open(tarball, "w:gz") as tar:
                    shutil.move(f, new_f)
                    tar.add(new_f)
                shutil.move(tarball, os.path.join(database, name[1:3]))
                os.unlink(new_f)

    # Remove the large temporary tmp directory
    shutil.rmtree(os.environ["CCP4_SCR"])
    os.environ["CCP4_SCR"] = ccp4_scr

    # Leave a timestamp
    leave_timestamp('simbad_sphere.txt')


def create_db_argparse():
    """Argparse function database creationg"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = p.add_subparsers(help='Database-specific commands')

    pa = sp.add_parser('lattice', help='lattice database')
    pa.set_defaults(which="lattice")
    simbad.command_line._argparse_lattice_options(pa)

    pb = sp.add_parser('morda', help='morda database')
    pb.set_defaults(which="morda")
    simbad.command_line._argparse_job_submission_options(pb)
    pb.add_argument('simbad_db', type=str, help='Path to local copy of the SIMBAD database')
    
    pc = sp.add_parser('sphere', help='sphere database')
    pc.set_defaults(which="sphere")
    simbad.command_line._argparse_job_submission_options(pc)
    pc.add_argument('-chunk_size', default=5000, type=int,
                    help='Max jobs to submit at any given time [disk space dependent')
    pc.add_argument('simbad_db', type=str, help='Path to local copy of the SIMBAD database')

    return p


def leave_timestamp(f):
    """Write the current date & time to a file"""
    with open(f, 'w') as f_out:
        f_out.write(str(datetime.datetime.now()))


def main():
    """SIMBAD database creation function"""
    args = create_db_argparse().parse_args()

    # Logger setup
    global logger
    args.debug_lvl = 'info'
    logger = simbad.command_line.setup_logging(level=args.debug_lvl)

    # Print a fancy header
    simbad.command_line.print_header()

    # Take a time snapshot
    time_start = time.time()

    # Create the requested database
    if args.which == "lattice":
        create_lattice_db(args.latt_db)
    elif args.which == "morda":
        create_morda_db(args.simbad_db, nproc=args.nproc, submit_cluster=args.submit_cluster,
                        submit_qtype=args.submit_qtype, submit_queue=args.submit_queue, 
                        submit_array=args.submit_array, submit_max_array=args.submit_max_array)
    elif args.which == "sphere":
        create_sphere_db(args.simbad_db, shres=3, nproc=args.nproc,
                         submit_cluster=args.submit_cluster, submit_qtype=args.submit_qtype,
                         submit_queue=args.submit_queue, submit_array=args.submit_array,
                         submit_max_array=args.submit_max_array, chunk_size=args.chunk_size)

    # Calculate and display the runtime 
    days, hours, mins, secs = simbad.command_line.calculate_runtime(time_start, time.time())
    logger.info("Database creation completed in %d days, %d hours, %d minutes, and %d seconds", days, hours, mins, secs)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
