"""Commad line facility for SIMBAD scripts"""

__author__ = "Felix Simkovic"
__date__ = "14 Apr 2017"
__version__ = "0.1"

from distutils.version import StrictVersion

import datetime
import logging
import os
import platform
import sys
import time

from mbkit.apps import EXE_EXT 
from mbkit.dispatch.cexectools import cexec

import simbad.constants
import simbad.lattice.search
import simbad.mr
import simbad.rotsearch.amore_search
import simbad.util.mtz_util
import simbad.version


def _argparse_core_options(p):
    """Add core options to an already existing parser"""
    sg = p.add_argument_group('Basic options')
    sg.add_argument('-ccp4_jobid', type=int,
                    help='Set the CCP4 job id - only needed when running from the CCP4 GUI')
    sg.add_argument('-debug_lvl', type=str, default='info',
                    help='The console verbosity level < notset | info | debug | warning | error | critical > ')
    sg.add_argument('-early_term', default=True,
                    help="Terminate the program early if a solution is found")
    sg.add_argument('-max_to_keep', type=int, default=20,
                    help="The maximum number of results to return")
    sg.add_argument('-name', type=str, default="simbad",
                    help='4-letter identifier for job [simb]')
    sg.add_argument('-run_dir', type=str, default=os.getcwd(),
                    help='Directory where the SIMBAD work directory will be created [current dir]')
    sg.add_argument('-work_dir', type=str,
                    help='Path to the directory where SIMBAD will run (will be created if it doesn\'t exist)')
    sg.add_argument('-webserver_uri',
                    help='URI of the webserver directory - also indicates we are running as a webserver')
    sg.add_argument('--version', action='version', version='SIMBAD v{0}'.format(simbad.version.__version__),
                    help='Print the SIMBAD version')
    sg.add_argument('-no_gui', default=False,
                    help="No simbad GUI")


def _argparse_job_submission_options(p):
    """Add the options for submission to a cluster queuing system"""
    sg = p.add_argument_group('Cluster queue submission options')
    sg.add_argument('-nproc', type=int, default=1,
                    help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors. "\
                         "For cluster submission, this should be the number of processors on a node.")
    sg.add_argument('-submit_qtype', type=str, default='local',
                    help='The job submission queue type [ local | sge ]')
    sg.add_argument('-submit_queue', type=str, default=None,
                    help='The queue to submit to on the cluster.')


def _argparse_contaminant_options(p):
    """Contaminant search specific options"""
    sg = p.add_argument_group('Contaminant search specific options')
    sg.add_argument('-cont_db', type=str, default=simbad.constants.CONTAMINANT_MODELS,
                    help='Path to local copy of the contaminant database')


def _argparse_morda_options(p):
    """Morda search specific options"""
    sg = p.add_argument_group('Morda search specific options')
    sg.add_argument('-morda_db', type=str,
                    help='Path to local copy of the MoRDa database')


def _argparse_lattice_options(p):
    """Lattice search specific options"""
    sg = p.add_argument_group('Lattice search specific options')
    sg.add_argument('-latt_db', type=str, default=simbad.constants.SIMBAD_LATTICE_DB,
                    help='Path to local copy of the lattice database')


def _argparse_rot_options(p):
    """Rotation search specific options"""
    sg = p.add_argument_group('AMORE Rotation search specific options')
    sg.add_argument("-npic", type=int, default=50,
                    help="Number of peaks to output from the translation function map for each orientation")
    sg.add_argument('-min_solvent_content', type=int, default=30,
                    help="The minimum solvent content present in the unit cell with the input model")
    sg.add_argument('-pklim', type=float, default=0.5,
                    help="Peak limit, output all peaks above <pklim>")
    sg.add_argument('-rotastep', type=float, default=1.0,
                    help="Size of rotation step")
    sg.add_argument('-shres', type=float, default=3.0,
                    help="Spherical harmonic resolution")


def _argparse_mr_options(p):
    sg = p.add_argument_group('Molecular Replacement specific options')
    sg.add_argument('-enan', default=False,
                    help='Check enantiomorphic space groups')
    sg.add_argument('-mr_keywords', type=str,
                    help='Path to file containing keywords for MR program')
    sg.add_argument('-refine_keywords', type=str,
                    help='Path to file containing keywords for the refinement program')
    sg.add_argument('-amore_exe', type=str, default=os.path.join(os.environ["CCP4"], 'bin', 'amoreCCB2.exe'),
                    help='Path to amore executable')
    sg.add_argument('-mr_program', type=str, default="molrep",
                    help='Path to the MR program to use. Options: < molrep | phaser >')
    sg.add_argument('-refine_program', type=str, default="refmac5",
                    help='Path to the refinement program to use. Options: < refmac5 >')
    sg.add_argument('-pdb_db', type=str,
                    help='Path to local copy of the PDB, this is needed if there is no internet access')


def _argparse_mtz_options(p):
    """Add MTZ specific options"""
    sg = p.add_argument_group('MTZ specific options')
    sg.add_argument('-F', type=str,
                    help='Flag for F column in the MTZ')
    sg.add_argument('-SIGF', type=str,
                    help='Flag for SIGF column in the MTZ')
    sg.add_argument('-FREE', type=str,
                    help='Flag for FREE column in the MTZ')
    sg.add_argument('-DANO', type=str,
                    help='Flag for the DANO column in the MTZ')
    sg.add_argument('-SIGDANO', type=str,
                    help='Flag for the SIGDANO column in the MTZ')


# This function is looking for a new home, any suggestions? - I suggest rotsearch.__init__
# Hold out for now until we find a better solution for args
def _simbad_contaminant_search(args):
    """A wrapper function to run the SIMBAD contaminant search
    
    Parameters
    ----------
    args : obj
       A :obj:`ArgumentParser <argparse.ArgumentParser>` object 

    Returns
    -------
    bool
       Successful or not

    """
    logger = logging.getLogger(__name__)
    stem = os.path.join(args.work_dir, 'cont')
    contaminant_model_dir = os.path.join(stem, 'contaminant_input_models')
    os.makedirs(contaminant_model_dir)
    rotation_search = simbad.rotsearch.amore_search.AmoreRotationSearch(
        args.amore_exe, args.mtz, stem, args.max_to_keep
    )
    rotation_search.sortfun()
    rotation_search.run_pdb(
        args.cont_db, output_model_dir=contaminant_model_dir, nproc=args.nproc, shres=args.shres,
        pklim=args.pklim, npic=args.npic, rotastep=args.rotastep, min_solvent_content=args.min_solvent_content,
        submit_qtype=args.submit_qtype, submit_queue=args.submit_queue,
    )
    if rotation_search.search_results:
        rot_summary_f = os.path.join(stem, 'rot_search.csv')
        logger.debug("Contaminant search summary file: %s", rot_summary_f)
        rotation_search.summarize(rot_summary_f)
        # Create directories for the contaminant search MR
        contaminant_output_dir = os.path.join(stem, 'mr_contaminant')
        # Run MR on results
        molecular_replacement = simbad.mr.MrSubmit(
            args.mtz, args.mr_program, args.refine_program, contaminant_model_dir, contaminant_output_dir,
            enant=args.enan
        )
        molecular_replacement.submit_jobs(rotation_search.search_results, nproc=args.nproc, early_term=args.early_term,
                                          submit_qtype=args.submit_qtype, submit_queue=args.submit_queue)
        mr_summary_f = os.path.join(stem, 'cont_mr.csv')
        logger.debug("Contaminant MR summary file: %s", mr_summary_f)
        molecular_replacement.summarize(mr_summary_f)
        if simbad.mr.mr_succeeded_csvfile(mr_summary_f):
            return True
    return False


def _simbad_morda_search(args):
    """A wrapper function to run the SIMBAD morda search
    
    Parameters
    ----------
    args : obj
       A :obj:`ArgumentParser <argparse.ArgumentParser>` object 

    Returns
    -------
    bool
       Successful or not

    """
    logger = logging.getLogger(__name__)
    stem = os.path.join(args.work_dir, 'morda')
    morda_model_dir = os.path.join(stem, 'morda_input_models')
    os.makedirs(morda_model_dir)

    if not args.morda_db:
        msg = "Failed to specify the location of the MoRDa database for the SIMBAD MoRDa run"
        logger.critical(msg)
        raise RuntimeError(msg)

    rotation_search = simbad.rotsearch.amore_search.AmoreRotationSearch(
        args.amore_exe, args.mtz, stem, 200
    )
    rotation_search.sortfun()
    rotation_search.run_pdb(
        args.morda_db, output_model_dir=morda_model_dir, nproc=args.nproc, shres=args.shres,
        pklim=args.pklim, npic=args.npic, rotastep=args.rotastep, min_solvent_content=args.min_solvent_content,
        submit_qtype=args.submit_qtype, submit_queue=args.submit_queue,
    )
    
    if rotation_search.search_results:
        rot_summary_f = os.path.join(stem, 'rot_search.csv')
        logger.debug("MoRDa search summary file: %s", rot_summary_f)
        rotation_search.summarize(rot_summary_f)
        # Create directories for the morda search MR
        morda_output_dir = os.path.join(stem, 'mr_morda')
        # Run MR on results
        molecular_replacement = simbad.mr.MrSubmit(
            args.mtz, args.mr_program, args.refine_program, morda_model_dir, morda_output_dir, enant=args.enan
        )
        molecular_replacement.submit_jobs(rotation_search.search_results, nproc=args.nproc, early_term=args.early_term,
                                          submit_qtype=args.submit_qtype, submit_queue=args.submit_queue)
        mr_summary_f = os.path.join(stem, 'morda_mr.csv')
        logger.debug("MoRDa search MR summary file: %s", mr_summary_f)
        molecular_replacement.summarize(mr_summary_f)
        if simbad.mr.mr_succeeded_csvfile(mr_summary_f):
            return True

    return False


# This function should really be moved to simbad/lattice/__init__.py
# but we hold out for now until we find a better solution for args
def _simbad_lattice_search(args):
    """A wrapper function to run the SIMBAD lattice search
    
    Parameters
    ----------
    args : obj
       A :obj:`ArgumentParser <argparse.ArgumentParser>` object 

    Returns
    -------
    bool
       Successful or not

    """
    logger = logging.getLogger(__name__)
    space_group, _, cell_parameters = simbad.util.mtz_util.crystal_data(args.mtz)
    stem = os.path.join(args.work_dir, 'latt')
    lattice_in_mod = os.path.join(stem, 'lattice_input_models')
    lattice_mr_dir = os.path.join(stem, 'mr_lattice')
    os.makedirs(lattice_in_mod)
    os.makedirs(lattice_mr_dir)
    lattice_search = simbad.lattice.search.LatticeSearch(
        cell_parameters, space_group, args.latt_db, max_to_keep=args.max_to_keep
    )
    lattice_search.search()
    if lattice_search.search_results:
        latt_summary_f = os.path.join(stem, 'lattice_search.csv')
        logger.debug("Lattice search summary file: %s", latt_summary_f)
        lattice_search.summarize(latt_summary_f)
        if args.pdb_db:
            lattice_search.copy_results(args.pdb_db, lattice_in_mod)
        else:
            lattice_search.download_results(lattice_in_mod)
        # Run MR on results
        molecular_replacement = simbad.mr.MrSubmit(
            args.mtz, args.mr_program, args.refine_program, lattice_in_mod, lattice_mr_dir, enant=args.enan
        )
        molecular_replacement.submit_jobs(lattice_search.search_results, nproc=args.nproc, early_term=args.early_term,
                                          submit_qtype=args.submit_qtype, submit_queue=args.submit_queue)
        mr_summary_f = os.path.join(stem, 'lattice_mr.csv')
        logger.debug("Lattice search MR summary file: %s", mr_summary_f)
        molecular_replacement.summarize(mr_summary_f)
        if simbad.mr.mr_succeeded_csvfile(mr_summary_f):
            return True
    return False


def calculate_runtime(start, stop):
    """Calculate the runtime in hours

    Parameters
    ----------
    start : float
       The start time
    stop : float
       The stop time

    Returns
    -------
    tuple
       The runtime in (weeks, days, hours, minutes, seconds)

    """
    diff = stop - start
    d = datetime.datetime(1970, 1, 1) + datetime.timedelta(seconds=diff)
    # Leave -1 in day as we start on the first day of the year
    return (d.day-1, d.hour, d.minute, d.second)


def ccp4_root():
    """Run CCP4 specific checks to check if it was setup correctly

    Returns
    -------
    str
       The path to the CCP4 root directory

    Raises
    ------
    KeyError
       Cannot find CCP4 installation
    KeyError
       $CCP4_SCR environment variable not set
    ValueError
       Cannot find the $CCP4_SCR directory

    """
    logger = logging.getLogger(__name__)
    if "CCP4" not in os.environ:
        msg = "Cannot find CCP4 installation - please make sure CCP4 " \
              "is installed and the setup scripts have been run!"
        logger.critical(msg)
        raise KeyError(msg)
    elif "CCP4_SCR" not in os.environ:
        msg = "$CCP4_SCR environment variable not set - please make sure " \
              "CCP4 is installed and the setup scripts have been run!"
        logger.critical(msg)
        raise KeyError(msg)
    elif not os.path.isdir(os.environ['CCP4_SCR']):
        msg = "Cannot find the $CCP4_SCR directory: {0}".format(os.environ["CCP4_SCR"])
        logger.critical(msg)
        raise ValueError(msg)
    return os.environ['CCP4']


def ccp4_version():
    """Identify the CCP4 version

    Returns
    -------
    obj
       A :obj:`StrictVersion` object containing the CCP4 version
    
    Raises
    ------
    RuntimeError
       Cannot determine CCP4 version

    """
    # Currently there seems no sensible way of doing this other then running a program and grepping the output
    stdout = cexec(['pdbcur' + EXE_EXT], permit_nonzero=True)
    tversion = None
    for line in stdout.split(os.linesep):
        if line.startswith(' ### CCP4'):
            tversion = line.split()[2].rstrip(':')
    if tversion is None:
        raise RuntimeError("Cannot determine CCP4 version")
    # Create the version as StrictVersion to make sure it's valid and allow for easier comparisons
    return StrictVersion(tversion)


def make_workdir(run_dir, ccp4_jobid=None, rootname='SIMBAD_'):
    """Make a work directory rooted at work_dir and return its path

    Parameters
    ----------
    run_dir : str
       The path to a run directory
    ccp4_jobid : int, optional
       CCP4-assigned job identifier
    rootname : str, optional
       Base name of the SIMBAD directory [default: \'SIMBAD_\']

    Returns
    -------
    work_dir : str
       The path to the working directory

    Raises
    ------
    ValueError
       There is an existing SIMBAD CCP4 work directory

    """
    if ccp4_jobid:
        dname = os.path.join(run_dir, rootname + str(ccp4_jobid))
        if os.path.exists(dname):
            raise ValueError("There is an existing SIMBAD CCP4 work directory: {0}\n"
                             "Please delete/move it aside.")
        os.mkdir(dname)
        return dname
    
    run_inc = 0
    while True:
        work_dir = os.path.join(run_dir, rootname + str(run_inc))
        if os.path.isdir(work_dir):
            run_inc += 1
        else:
            os.mkdir(work_dir)
            return work_dir


def print_header():
    """Print the header information at the top of each script"""
    logger = logging.getLogger(__name__)
    # When changing the `line` text make sure it does not exceed 118 characters, otherwise adjust nhashes
    nhashes = 120
    logger.info("%(sep)s%(hashish)s%(sep)s%(hashish)s%(sep)s%(hashish)s%(sep)s#%(line)s#%(sep)s%(hashish)s%(sep)s",
                {'hashish': '#' * nhashes, 'sep': os.linesep,
                'line': 'SIMBAD - Sequence Independent Molecular '
                        'replacement Based on Available Database'.center(nhashes-2, ' ')}
    )
    logger.info("SIMBAD version: %s", simbad.version.__version__)
    logger.info("Running with CCP4 version: %s from directory: %s", ccp4_version(), ccp4_root())
    logger.info("Running on host: %s", platform.node())
    logger.info("Running on platform: %s", platform.platform())
    logger.info("Job started at: %s", time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
    logger.info("Invoked with command-line:\n%s\n", " ".join(map(str, sys.argv)))


def setup_logging(level='info', logfile=None):
    """Set up logging to the console for the root logger.

    Parameters
    ----------
    level : str, optional
       The console logging level to be used [default: info]
       To change, use one of 
           [ notset | info | debug | warning | error | critical ]
    logfile : str, optional
       The path to a full file log

    Returns
    -------
    logger
       Instance of a :obj:`logger <logging.Logger>`

    """
    # Reset any Handlers or Filters already in the logger to start from scratch
    # https://stackoverflow.com/a/16966965
    map(logging.getLogger().removeHandler, logging.getLogger().handlers[:])
    map(logging.getLogger().removeFilter, logging.getLogger().filters[:])

    logging_levels = {
        'notset': logging.NOTSET, 'info': logging.INFO, 'debug': logging.DEBUG,
        'warning': logging.WARNING, 'error': logging.ERROR, 'critical': logging.CRITICAL
    }
    
    # Create logger and default settings
    logging.getLogger().setLevel(logging.DEBUG)
    
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging_levels.get(level, logging.INFO))
    ch.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger().addHandler(ch)

    # create file handler which logs even debug messages
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.NOTSET)
        fh.setFormatter(
            logging.Formatter('%(asctime)s\t%(name)s [%(lineno)d]\t%(levelname)s\t%(message)s')
        )
        logging.getLogger().addHandler(fh)

    logging.getLogger().debug('Console logger level: %s', logging_levels.get(level, logging.INFO))
    logging.getLogger().debug('File logger level: %s', logging.NOTSET)

    return logging.getLogger()
