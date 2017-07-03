"""Script to download or update SIMBAD-related databases"""

__author__ = "Felix Simkovic & Adam Simpkin"
__date__ = "17 May 2017"
__version__ = "1.0"

import argparse
import base64
import datetime
import glob
import numpy as np
import os
import shutil
import sys
import tarfile
import urllib2
import zlib

from pyjob import Job
from pyjob.misc import StopWatch, make_script, tmp_dir 

import cctbx.crystal

import simbad.command_line
import simbad.exit
import simbad.rotsearch.amore_search

logger = None

# The space groups in the list below cannot be recognized by CCTBX, so we convert them
# to similar ones understandable by the library
CCTBX_ERROR_SG = {
    'A1': 'P1', 'B2': 'B112', 'C1211': 'C2', 'F422': 'I422', 'I21': 'I2', 'I1211': 'I2', 
    'P21212A': 'P212121', 'R3': 'R3:R', 'C4212': 'P422',
}

SYS_PLATFORM = sys.platform
CUSTOM_PLATFORM = "linux" if SYS_PLATFORM in ["linux", "linux2"] \
                   else "mac" if SYS_PLATFORM in ["darwin"] \
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
    local_db = os.path.basename(url)
    # http://stackoverflow.com/a/34831866/3046533
    chunk_size = 1 << 20
    with open(local_db, "wb") as f_out:
        query = urllib2.urlopen(url)
        while True:
            chunk = query.read(chunk_size)
            if not chunk:
                break
            f_out.write(chunk)
    # Extract relevant files from the tarball
    with tarfile.open(local_db) as tar:
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
    os.remove(local_db)
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
        pdb_code, rest = line[1:-1].split('","', 1)
        unit_cell, space_group, exp_tech = rest.rsplit('","', 2)
        unit_cell = unit_cell.replace('","', ',')
        space_group = space_group.replace(" ", "").strip()
        
        # Ignore non-xtal structures
        if not "X-RAY DIFFRACTION" in exp_tech.strip().upper():
            continue
        
        # Some entries do not have stored unit cell parameters
        try:
            unit_cell = map(float, unit_cell.split(','))
        except ValueError as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            error_count += 1
            continue
        space_group = CCTBX_ERROR_SG.get(space_group, space_group)   
        try:
            symmetry = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group, 
                                              correct_rhombohedral_setting_if_necessary=True)
        except Exception as e:
            logger.debug('Skipping pdb entry %s\t%s', pdb_code, e)
            error_count += 1
            continue
        crystal_data.append((pdb_code, symmetry))
    logger.debug('Error processing %d pdb entries', error_count)

    logger.info('Calculating the Niggli cells')
    niggli_data = np.zeros((len(crystal_data), 11))
    # Leave this as list, .append is faster than .vstack
    alt_niggli_data = []
    for i, xtal_data in enumerate(crystal_data):
        niggli_data[i][:4] = np.fromstring(xtal_data[0], dtype='uint8').astype(np.float64)
        niggli_data[i][4] = ord('\x00')
        niggli_data[i][5:] = np.asarray(xtal_data[1].niggli_cell().unit_cell().parameters())
        a, b, c, alpha, beta, gamma = niggli_data[i][5:]

        # Add alternate niggli cell where a and b may be flipped
        if np.allclose(a, b, atol=0.1) and not np.allclose(alpha, beta, atol=5):
            alt_niggli_data += [np.concatenate((niggli_data[i][:4], np.array([ord('*'), b, a, c, beta, alpha, gamma])))]
            
        # Add alternate niggli cell where b and c may be flipped
        if np.allclose(b, c, atol=0.1) and not np.allclose(beta, gamma, atol=5):
            alt_niggli_data += [np.concatenate((niggli_data[i][:4], np.array([ord('*'), a, c, b, alpha, gamma, beta])))]
        
    niggli_data = np.vstack([niggli_data, np.asarray(alt_niggli_data)])
    logger.info("Total Niggli cells loaded: %d", niggli_data.shape[0])

    if not database.endswith('.npz'):
        database += ".npz"

    logger.info('Storing database in file: %s', database)
    np.savez_compressed(database, niggli_data)


def create_morda_db(database, nproc=2, submit_qtype=None, submit_queue=False, chunk_size=5000):
    """Create the MoRDa search database

    Parameters
    ----------
    database : str
       The path to the database folder
    nproc : int, optional
       The number of processors [default: 2]
    submit_qtype : str
       The cluster submission queue type - currently support SGE and LSF
    submit_queue : str
       The queue to submit to on the cluster
    chunk_size : int, optional
       The number of jobs to submit at the same time [default: 5000]
    
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
        logger.info('SIMBAD dat database up-to-date')
        shutil.rmtree(os.environ['MRD_DB'])
        leave_timestamp(os.path.join(database, 'simbad_morda.txt'))
        return
    else:
        logger.info("%d new entries were found in the MoRDa database, update SIMBAD database", len(dat_files))

    # Get the "get_model" script to extract the xyz coordinates
    exe = os.path.join(os.environ['MRD_DB'], "bin_" + CUSTOM_PLATFORM, "get_model")
    
    # Creating temporary output directory
    run_dir = tmp_dir(directory=os.getcwd())

    # Submit in chunks, so we don't take too much disk space
    # and can terminate without loosing the processed data
    total_chunk_cycles = len(dat_files) // chunk_size + (len(dat_files) % 5 > 0)
    for cycle, i in enumerate(range(0, len(dat_files), chunk_size)):
        logger.info("Working on chunk %d out of %d", cycle + 1, total_chunk_cycles)
        chunk_dat_files = dat_files[i:i + chunk_size]

        # Create the database files
        what_to_do = []
        for f in chunk_dat_files:
            code = os.path.basename(f).rsplit('.', 1)[0]
            final_file = os.path.join(database, code[1:3], code + ".dat")
            # We need a temporary directory within because "get_model" uses non-unique file names
            tmp_d = tmp_dir(directory=run_dir)
            get_model_output = os.path.join(tmp_d, code + ".pdb")
            # Prepare script for multiple submissions
            script = make_script(
                [["export CCP4_SCR=", tmp_d],
                 ["export MRD_DB=" + os.environ['MRD_DB']],
                 ["cd", tmp_d],
                 [exe, "-c", code, "-m", "d"]],
                directory=tmp_d
            )
            log = script.rsplit('.', 1)[0] + '.log'
            what_to_do += [(script, log, tmp_d, (get_model_output, final_file))]

        # Run the scripts
        scripts, _, tmps, files = zip(*what_to_do)
        j = Job(submit_qtype)
        j.submit(scripts, name='morda_db', nproc=nproc, submit_queue=submit_queue)
        j.wait()

        # Create PDB-like database subdirectories
        sub_dir_names = set([os.path.basename(f).rsplit('.', 1)[0][1:3] for f in chunk_dat_files])
        for sub_dir_name in sub_dir_names:
            sub_dir = os.path.join(database, sub_dir_name)
            if os.path.isdir(sub_dir):
                continue
            os.makedirs(sub_dir)

        # Move created files to database
        for output, final in files:
            with open(final, 'wb') as f_out:
                content = base64.b64encode(zlib.compress(open(output, 'r').read()))
                f_out.write(content)

        # Remove temporary directories to avoid "too many links" OS Error
        for d in tmps:
            shutil.rmtree(d)

    # Remove temporary files
    shutil.rmtree(run_dir)
    shutil.rmtree(os.environ['MRD_DB'])

    # Leave a timestamp
    leave_timestamp(os.path.join(database, 'simbad_morda.txt'))


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
    pb.add_argument('-chunk_size', default=5000, type=int,
                    help='Max jobs to submit at any given time [disk space dependent')
    pb.add_argument('simbad_db', type=str, help='Path to local copy of the SIMBAD database')
    
    return p


def leave_timestamp(f):
    """Write the current date & time to a file"""
    with open(f, 'w') as f_out:
        f_out.write(str(datetime.datetime.now()))


def main():
    """SIMBAD database creation function"""
    p = create_db_argparse()
    p.add_argument('-debug_lvl', type=str, default='info',
                   help='The console verbosity level < notset | info | debug | warning | error | critical > ')
    args = p.parse_args()

    # Logger setup
    global logger
    logger = simbad.command_line.setup_logging(level=args.debug_lvl)

    # Print a fancy header
    simbad.command_line.print_header()

    # Start taking time
    stopwatch = StopWatch()
    stopwatch.start()

    # Create the requested database
    if args.which == "lattice":
        create_lattice_db(args.latt_db)
    elif args.which == "morda":
        create_morda_db(args.simbad_db, nproc=args.nproc, submit_qtype=args.submit_qtype, submit_queue=args.submit_queue, 
                        chunk_size=args.chunk_size)

    # Calculate and display the runtime in hours
    stopwatch.stop()
    logger.info("Database creation completed in %d days, %d hours, %d minutes, and %d seconds", *stopwatch.time_pretty)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        simbad.exit.exit_error(*sys.exc_info())
