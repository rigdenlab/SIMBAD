# -*- coding: utf-8 -*-
"""Commad line facility for SIMBAD scripts"""

__author__ = "Felix Simkovic"
__date__ = "14 Apr 2017"
__version__ = "0.1"

from distutils.version import StrictVersion

import argparse
import logging
import os
import platform
import sys
import time

from enum import Enum

from pyjob import cexec
from pyjob.factory import TASK_PLATFORMS
from pyjob.script import EXE_EXT
from simbad.mr.options import MrPrograms, RefPrograms, SGAlternatives
from simbad.util import SIMBAD_DIRNAME

import simbad.db
import simbad.util.mtz_util
import simbad.version

if os.name != "nt":
    if "SSL_CERT_FILE" not in os.environ:
        os.environ["SSL_CERT_FILE"] = os.path.join(os.environ["CCP4"], "etc", "ssl", "cacert.pem")


class CCP4(object):
    """Wrapper class for CCP4 installation"""

    def __init__(self):
        self.root = CCP4RootDirectory()
        self.version = CCP4Version()


class CCP4RootDirectory(object):
    """The CCP4 root directory"""

    def __init__(self):
        if "CCP4" not in os.environ:
            raise KeyError("Cannot find CCP4 installation - please make sure CCP4 " + "is installed and the setup scripts have been run!")
        elif "CCP4_SCR" not in os.environ:
            raise KeyError("$CCP4_SCR environment variable not set - please make sure " + "CCP4 is installed and the setup scripts have been run!")
        elif not os.path.isdir(os.environ["CCP4_SCR"]):
            raise ValueError("Cannot find the $CCP4_SCR directory: {0}".format(os.environ["CCP4_SCR"]))
        else:
            self._root = os.environ["CCP4"]

    def __str__(self):
        return self._root

    def __repr__(self):
        return "{}: {}".format(self.__class__.__name__, self._root)


class CCP4Version(StrictVersion):
    """The CCP4 version class"""

    def __init__(self):
        StrictVersion.__init__(self)
        ccp4_major_minor = os.path.join(os.environ["CCP4"], "lib", "ccp4", "MAJOR_MINOR")
        if os.path.isfile(ccp4_major_minor):
            with open(ccp4_major_minor, "r") as f_in:
                tversion = f_in.read().strip()
        else:
            logger = logging.getLogger(__name__)
            logger.debug("Detecting CCP4 version via executing pdbcur")
            stdout = cexec(["pdbcur" + EXE_EXT], permit_nonzero=True)
            tversion = None
            for line in stdout.split(os.linesep):
                if line.startswith(" ### CCP4"):
                    tversion = line.split()[2].rstrip(":")
            if tversion is None:
                raise RuntimeError("Cannot determine CCP4 version")
        self.parse(tversion)


def is_valid_file(parser, arg):
    if os.path.exists(arg):
        return arg
    else:
        parser.error("The file %s does not exist!" % arg)


def is_valid_dir(parser, arg):
    if os.path.isdir(arg):
        return arg
    else:
        parser.error("The directory %s does not exist!" % arg)


class LogColors(Enum):
    """Color container for log messages"""

    CRITICAL = 31
    DEBUG = 34
    DEFAULT = 0
    ERROR = 31
    WARNING = 33


class LogLevels(Enum):
    """Log level container"""

    DEBUG = logging.DEBUG
    ERROR = logging.ERROR
    INFO = logging.INFO
    NOTSET = logging.NOTSET
    WARNING = logging.WARNING


class LogColorFormatter(logging.Formatter):
    """Formatter for log messages"""

    def format(self, record):
        if record.levelname in LogColors.__members__:
            prefix = "\033[1;{}m".format(LogColors[record.levelname].value)
            postfix = "\033[{}m".format(LogColors["DEFAULT"].value)
            record.msg = os.linesep.join([prefix + msg + postfix for msg in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)


class LogController(object):
    """Controller class for log messaging"""

    def __init__(self, reset=True):
        logging.getLogger().setLevel(logging.NOTSET)
        self._custom_added = False

    def add_console(self, level="info", format="%(message)s", stream=sys.stdout):
        levelname = self.get_levelname(level)
        ch = logging.StreamHandler(stream=stream)
        ch.setLevel(levelname)
        ch.setFormatter(LogColorFormatter(format))
        if not self._custom_added:
            self.reset()
        logging.getLogger().addHandler(ch)
        self._custom_added = True

    def add_logfile(self, file, level="info", format="%(message)s"):
        levelname = self.get_levelname(level)
        fh = logging.FileHandler(file)
        fh.setLevel(levelname)
        fh.setFormatter(logging.Formatter(format))
        if not self._custom_added:
            self.reset()
        logging.getLogger().addHandler(fh)
        self._custom_added = True

    def get_levelname(self, level):
        level_uc = level.upper()
        if LogController.level_valid(level_uc):
            return LogLevels[level_uc].value
        else:
            raise ValueError("Please provide a valid log level - %s is not!" % level)

    def get_logger(self):
        return logging.getLogger()

    def close(self):
        for h in logging.getLogger().handlers[:]:
            h.close()
            logging.getLogger().removeHandler(h)

    def reset(self):
        map(logging.getLogger().removeHandler, logging.getLogger().handlers[:])
        map(logging.getLogger().removeFilter, logging.getLogger().filters[:])
        self._custom_added = False

    @staticmethod
    def level_valid(level):
        return level in LogLevels.__members__


def _argparse_core_options(p):
    """Add core options to an already existing parser"""
    sg = p.add_argument_group("Basic options")
    sg.add_argument("-ccp4_jobid", type=int, help="Set the CCP4 job id - only needed when running from the CCP4 GUI")
    sg.add_argument("-ccp4i2_xml", help=argparse.SUPPRESS)
    sg.add_argument("-chunk_size", default=0, type=int, help="Max jobs to submit at any given time")
    sg.add_argument(
        "-debug_lvl", type=str, default="info", choices=["info", "debug", "warning", "error", "critical"], help="The console verbosity level"
    )
    sg.add_argument("-name", type=str, default="simbad", help="The identifier for each job [simbad]")
    sg.add_argument("-output_pdb", type=str, help="Path to the output PDB for the best result")
    sg.add_argument("-output_mtz", type=str, help="Path to the output MTZ for the best result")
    sg.add_argument("-run_dir", type=str, default=".", help="Directory where the SIMBAD work directory will be created")
    sg.add_argument("-results_to_display", type=int, default=10, help="The number of results to display in the GUI")
    sg.add_argument("-tmp_dir", type=str, help="Directory in which to put temporary files from SIMBAD")
    sg.add_argument("-work_dir", type=str, help="Path to the directory where SIMBAD will run (will be created if it doesn't exist)")
    sg.add_argument("-webserver_uri", help="URI of the webserver directory - also indicates we are running as a webserver")
    sg.add_argument("-rvapi_document", help=argparse.SUPPRESS)
    sg.add_argument("-tab_prefix", type=str, default="", help=argparse.SUPPRESS)
    sg.add_argument("--benchmark", default=False, action="store_true", help="Store all data for benchmarking")
    sg.add_argument("--cleanup", default=False, action="store_true", help="Delete all data not reported by the GUI")
    sg.add_argument("--display_gui", default=False, action="store_true", help="Show the SIMBAD GUI")
    sg.add_argument("--process_all", default=False, action="store_true", help="Trial all search models")
    sg.add_argument("--skip_mr", default=False, action="store_true", help="Skip Molecular replacement step")
    sg.add_argument("--version", action="version", version="SIMBAD v{0}".format(simbad.version.__version__), help="Print the SIMBAD version")


def _argparse_job_submission_options(p):
    """Add the options for submission to a cluster queuing system"""
    sg = p.add_argument_group("Cluster queue submission options")
    sg.add_argument(
        "-nproc",
        type=int,
        default=1,
        help="Number of processors. For local, serial runs the jobs will be split across nproc "
        "processors. For cluster submission, this should be the number of processors on a node.",
    )
    sg.add_argument(
        "-submit_nproc",
        type=int,
        default=1,
        help="For cluster submission, the number of processors to use on head node when creating " "submission scripts",
    )
    sg.add_argument("-submit_qtype", type=str, default="local", choices=TASK_PLATFORMS.keys(), help="The job submission queue type")
    sg.add_argument("-submit_queue", type=str, default=None, help="The queue to submit to on the cluster.")


def _argparse_contaminant_options(p):
    """Contaminant search specific options"""
    sg = p.add_argument_group("Contaminant search specific options")
    sg.add_argument(
        "-cont_db", type=lambda x: is_valid_dir(sg, x), default=simbad.CONTAMINANT_MODELS, help="Path to local copy of the contaminant database"
    )
    sg.add_argument("-max_contaminant_results", type=int, default=20, help="The maximum number of contaminant results to return")
    sg.add_argument("-organism", type=str, help="Select a specific host organism using the UniProt mnemonic organism" "identification code")


def _argparse_morda_options(p):
    """Morda search specific options"""
    sg = p.add_argument_group("Morda search specific options")
    sg.add_argument("-morda_db", type=lambda x: is_valid_dir(sg, x), default=simbad.MORDA_MODELS, help="Path to local copy of the MoRDa database")
    sg.add_argument("-max_morda_results", type=int, default=200, help="The maximum number of contaminant results to return")


def _argparse_lattice_options(p):
    """Lattice search specific options"""
    sg = p.add_argument_group("Lattice search specific options")
    sg.add_argument("-latt_db", type=lambda x: is_valid_file(sg, x), default=simbad.LATTICE_DB, help="Path to local copy of the lattice database")
    sg.add_argument("-max_lattice_results", type=int, default=20, help="The maximum number of lattice results to return")
    sg.add_argument("-max_penalty_score", type=int, default=7, help="The maximum lattice penalty score allowed")


def _argparse_rot_options(p):
    """Rotation search specific options"""
    sg = p.add_argument_group("AMORE Rotation search specific options")
    sg.add_argument("-amore_exe", type=str, default=os.path.join(os.environ["CCP4"], "libexec", "amore-rs"), help="Path to amore executable")
    sg.add_argument("-npic", type=int, default=50, help="Number of peaks to output from the translation function map for each orientation")
    sg.add_argument("-min_solvent_content", type=int, default=30, help="The minimum solvent content present in the unit cell with the input model")
    sg.add_argument("-pklim", type=float, default=0.5, help="Peak limit, output all peaks above <pklim>")
    sg.add_argument("-rotastep", type=float, default=1.0, help="Size of rotation step")
    sg.add_argument("-rot_program", type=str, default="amore", choices=["amore", "phaser"], help="Program with which to perform rotation search")
    sg.add_argument("-shres", type=float, default=3.0, help="Spherical harmonic resolution")


def _argparse_mr_options(p):
    sg = p.add_argument_group("Molecular Replacement specific options")
    sg.add_argument("-sga", "--sgalternative", default="none", choices=SGAlternatives.__members__.keys(), help="Check alternative space groups")
    sg.add_argument("-mr_program", type=str, default="molrep", choices=MrPrograms.__members__.keys(), help="Path to the MR program to use.")
    sg.add_argument(
        "-refine_program", type=str, default="refmac5", choices=RefPrograms.__members__.keys(), help="Path to the refinement program to use."
    )
    sg.add_argument("-refine_cycles", type=int, help="The number of refinement cycles to run [default: 30]")
    sg.add_argument("-pdb_db", type=str, help="Path to local copy of the PDB, this is needed if there is no internet access")
    sg.add_argument("-phaser_kill", type=int, default=30, help="Time in minutes after which phaser will be killed (0 to leave running)")


def _argparse_mtz_options(p):
    """Add MTZ specific options"""
    sg = p.add_argument_group("MTZ specific options")
    sg.add_argument("-F", type=str, help="Flag for F column in the MTZ")
    sg.add_argument("-SIGF", type=str, help="Flag for SIGF column in the MTZ")
    sg.add_argument("-FREE", type=str, help="Flag for FREE column in the MTZ")
    sg.add_argument("-DANO", type=str, help="Flag for the DANO column in the MTZ")
    sg.add_argument("-SIGDANO", type=str, help="Flag for the SIGDANO column in the MTZ")


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
    stem = os.path.join(args.work_dir, "cont")
    os.makedirs(stem)

    if args.organism:
        organism_cont_db = os.path.join(args.cont_db, args.organism.upper())
        if not os.path.isdir(organism_cont_db):
            msg = "Cannot find organism {0} in contaminant database. Running full contaminant database".format(args.organism)
            logging.debug(msg)
        else:
            args.cont_db = organism_cont_db

    temp_mtz = os.path.join(args.work_dir, "input.mtz")
    if os.path.isfile(temp_mtz):
        pass
    else:
        simbad.util.mtz_util.ctruncate(args.mtz, temp_mtz)

    from simbad.rotsearch import rotation_search_factory

    rotation_obj = rotation_search_factory(args.rot_program)
    rotation_search = rotation_obj(
        temp_mtz,
        args.mr_program,
        args.tmp_dir,
        stem,
        amore_exe=args.amore_exe,
        max_to_keep=args.max_contaminant_results,
        skip_mr=args.skip_mr,
        process_all=args.process_all,
    )
    rotation_search.run(
        os.path.abspath(args.cont_db),
        nproc=args.nproc,
        shres=args.shres,
        pklim=args.pklim,
        npic=args.npic,
        rotastep=args.rotastep,
        min_solvent_content=args.min_solvent_content,
        submit_nproc=args.submit_nproc,
        submit_qtype=args.submit_qtype,
        submit_queue=args.submit_queue,
        chunk_size=args.chunk_size,
    )

    if rotation_search.search_results:
        rot_summary_f = os.path.join(stem, "rot_search.csv")
        logger.debug("Contaminant search summary file: %s", rot_summary_f)
        rotation_search.summarize(rot_summary_f)

        if args.benchmark:
            rot_all_f = os.path.join(stem, "rot_search_all.csv")
            logger.debug("Contaminant search file containing all results: %s", rot_summary_f)
            rotation_search.summarize_all(rot_all_f)

    if rotation_search.search_results and not args.skip_mr:
        from simbad.mr import mr_succeeded_csvfile

        contaminant_mr_dir = os.path.join(stem, "mr_search")

        if args.refine_cycles:
            refine_cycles = args.refine_cycles
        else:
            refine_cycles = 30

        molecular_replacement = submit_mr_jobs(temp_mtz, contaminant_mr_dir, rotation_search.search_results, None, refine_cycles, args)
        mr_summary_f = os.path.join(stem, "cont_mr.csv")
        logger.debug("Contaminant MR summary file: %s", mr_summary_f)
        molecular_replacement.summarize(mr_summary_f)

        output_dir = os.path.join(args.work_dir, "output_files")
        make_output_dir(stem, output_dir, mr_summary_f, args.mr_program)

        if args.cleanup:
            cleanup(stem, os.path.join(args.work_dir, "output_files"), mr_summary_f, args.results_to_display)

        if mr_succeeded_csvfile(mr_summary_f):
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
    stem = os.path.join(args.work_dir, "morda")
    os.makedirs(stem)

    if not args.morda_db:
        msg = "Failed to specify the location of the MoRDa database for the SIMBAD MoRDa run"
        logger.critical(msg)
        raise RuntimeError(msg)

    temp_mtz = os.path.join(args.work_dir, "input.mtz")
    if os.path.isfile(temp_mtz):
        pass
    else:
        simbad.util.mtz_util.ctruncate(args.mtz, temp_mtz)

    from simbad.rotsearch import rotation_search_factory

    rotation_obj = rotation_search_factory(args.rot_program)
    rotation_search = rotation_obj(
        temp_mtz,
        args.mr_program,
        args.tmp_dir,
        stem,
        amore_exe=args.amore_exe,
        max_to_keep=args.max_morda_results,
        skip_mr=args.skip_mr,
        process_all=args.process_all,
    )
    rotation_search.run(
        args.morda_db,
        nproc=args.nproc,
        shres=args.shres,
        pklim=args.pklim,
        npic=args.npic,
        rotastep=args.rotastep,
        min_solvent_content=args.min_solvent_content,
        submit_nproc=args.submit_nproc,
        submit_qtype=args.submit_qtype,
        submit_queue=args.submit_queue,
        chunk_size=args.chunk_size,
    )

    if rotation_search.search_results:
        rot_summary_f = os.path.join(stem, "rot_search.csv")
        logger.debug("MoRDa search summary file: %s", rot_summary_f)
        rotation_search.summarize(rot_summary_f)

        if args.benchmark:
            rot_all_f = os.path.join(stem, "rot_search_all.csv")
            logger.debug("MoRDa search file containing all results: %s", rot_summary_f)
            rotation_search.summarize_all(rot_all_f)

    if rotation_search.search_results and not args.skip_mr:
        from simbad.mr import mr_succeeded_csvfile

        morda_mr_dir = os.path.join(stem, "mr_search")

        if args.refine_cycles:
            refine_cycles = args.refine_cycles
        else:
            refine_cycles = 30

        molecular_replacement = submit_mr_jobs(temp_mtz, morda_mr_dir, rotation_search.search_results, "jelly_body", refine_cycles, args)
        mr_summary_f = os.path.join(stem, "morda_mr.csv")
        logger.debug("MoRDa search MR summary file: %s", mr_summary_f)
        molecular_replacement.summarize(mr_summary_f)

        output_dir = os.path.join(args.work_dir, "output_files")
        make_output_dir(stem, output_dir, mr_summary_f, args.mr_program)

        if args.cleanup:
            cleanup(stem, os.path.join(args.work_dir, "output_files"), mr_summary_f, args.results_to_display)

        if mr_succeeded_csvfile(mr_summary_f):
            return True

    return False


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
    from simbad.lattice.lattice_search import LatticeSearch

    MTZ_AVAIL = args.mtz is not None
    temp_mtz = None

    logger = logging.getLogger(__name__)
    if MTZ_AVAIL:
        temp_mtz = os.path.join(args.work_dir, "input.mtz")
        simbad.util.mtz_util.ctruncate(args.mtz, temp_mtz)
        space_group, _, cell_parameters = simbad.util.mtz_util.crystal_data(temp_mtz)
    else:
        space_group = args.space_group
        cell_parameters = args.unit_cell.replace(",", " ")
        cell_parameters = (float(i) for i in cell_parameters.split())

    stem = os.path.join(args.work_dir, "latt")
    lattice_mr_dir = os.path.join(stem, "mr_search")
    lattice_mod_dir = os.path.join(lattice_mr_dir, "mr_models")
    os.makedirs(lattice_mod_dir)

    ls = LatticeSearch(args.latt_db, lattice_mod_dir)
    ls.search(space_group, cell_parameters, max_to_keep=args.max_lattice_results, max_penalty=args.max_penalty_score)

    if len(ls.results) > 0:
        latt_summary_f = os.path.join(stem, "lattice_search.csv")
        ls.summarize(latt_summary_f)

        if args.pdb_db:
            ls.copy_results(args.pdb_db, lattice_mod_dir)
        else:
            ls.download_results(lattice_mod_dir)

        # Check so we don't attempt MR when download/copying failed
        if len(ls.results) < 1:
            return False

        if MTZ_AVAIL and not args.skip_mr:
            from simbad.mr import mr_succeeded_csvfile

            if args.refine_cycles:
                refine_cycles = args.refine_cycles
            else:
                refine_cycles = 0

            molecular_replacement = submit_mr_jobs(temp_mtz, lattice_mr_dir, ls.results, None, refine_cycles, args)
            mr_summary_f = os.path.join(stem, "lattice_mr.csv")
            logger.debug("Lattice search MR summary file: %s", mr_summary_f)
            molecular_replacement.summarize(mr_summary_f)

            output_dir = os.path.join(args.work_dir, "output_files")
            make_output_dir(stem, output_dir, mr_summary_f, args.mr_program)

            if args.cleanup:
                cleanup(stem, os.path.join(args.work_dir, "output_files"), mr_summary_f, args.results_to_display)

            if mr_succeeded_csvfile(mr_summary_f):
                return True

    return False


def cleanup(directory, output_dir, csv, results_to_keep):
    """Function to clean up working directory results not reported by GUI"""
    import pandas as pd
    import shutil

    df = pd.read_csv(csv)
    data = df.pdb_code.tolist()

    if len(data) > results_to_keep:
        for i in data[results_to_keep:]:
            shutil.rmtree(os.path.join(output_dir, i))

    mr_dir = os.path.join(directory, "mr_search")
    if os.path.isdir(mr_dir):
        shutil.rmtree(mr_dir)


def get_work_dir(run_dir, work_dir=None, ccp4_jobid=None, ccp4i2_xml=None):
    """Figure out the relative working directory by provided options"""
    if work_dir:
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)
    elif run_dir and os.path.isdir(run_dir):
        work_dir = make_workdir(run_dir, ccp4_jobid=ccp4_jobid, ccp4i2_xml=ccp4i2_xml)
    elif run_dir:
        os.mkdir(run_dir)
        work_dir = make_workdir(run_dir, ccp4_jobid=ccp4_jobid)
    else:
        raise RuntimeError("Not entirely sure what has happened here " + "but I should never get to here")
    return os.path.abspath(work_dir)


def make_output_dir(run_dir, output_dir, csv, mr_program):
    """Make a directory containing output files

    Parameters
    ----------
    run_dir : str
        The path to the run directory
    output_dir : str
        The path to the output directory to create
    csv : file
        CSV file containing results
    mr_program : str
        The MR program used
    """
    import pandas as pd
    import shutil

    df = pd.read_csv(csv)
    data = df.pdb_code.tolist()

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for pdb_code in data:
        pdb_output_path = os.path.join(output_dir, pdb_code)
        if not os.path.isdir(pdb_output_path):
            os.mkdir(pdb_output_path)
            mr_workdir = os.path.join(run_dir, "mr_search", pdb_code, "mr", mr_program)

            files_to_copy = [
                os.path.join(mr_workdir, "{0}_mr.log".format(pdb_code)),
                os.path.join(mr_workdir, "refine", "{0}_refinement_output.pdb".format(pdb_code)),
                os.path.join(mr_workdir, "refine", "{0}_refinement_output.mtz".format(pdb_code)),
                os.path.join(mr_workdir, "refine", "{0}_ref.log".format(pdb_code)),
            ]
            for f in files_to_copy:
                shutil.copy(f, pdb_output_path)


def make_workdir(run_dir, ccp4_jobid=None, ccp4i2_xml=None, rootname=SIMBAD_DIRNAME + "_"):
    """Make a work directory rooted at work_dir and return its path

    Parameters
    ----------
    run_dir : str
       The path to a run directory
    ccp4_jobid : int, optional
       CCP4-assigned job identifier
    ccp4i2_xml : str, optional
       Path to CCP4 I2 XML output (just used to indicate running under CCP4 I2)
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
    if ccp4i2_xml:
        dname = os.path.join(run_dir, SIMBAD_DIRNAME)
        os.mkdir(dname)
        return dname
    elif ccp4_jobid:
        dname = os.path.join(run_dir, rootname + str(ccp4_jobid))
        if os.path.exists(dname):
            raise ValueError("There is an existing SIMBAD CCP4 work directory: {0}\n" "Please delete/move it aside.")
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
    nhashes = 120
    ccp4 = CCP4()
    logger.info(
        "%(sep)s%(hashish)s%(sep)s%(hashish)s%(sep)s%(hashish)s%(sep)s#%(line)s#%(sep)s%(hashish)s%(sep)s",
        {
            "hashish": "#" * nhashes,
            "sep": os.linesep,
            "line": "SIMBAD - Sequence Independent Molecular " "replacement Based on Available Database".center(nhashes - 2, " "),
        },
    )
    logger.info("SIMBAD version: %s", simbad.version.__version__)
    logger.info("Running with CCP4 version: %s from directory: %s", ccp4.version, ccp4.root)
    logger.info("Running on host: %s", platform.node())
    logger.info("Running on platform: %s", platform.platform())
    logger.info("Job started at: %s", time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
    script_name = os.path.basename(sys.argv[0]).replace(".py", "").replace("_", "-").replace("-main", "")
    logger.info("Invoked with command-line:\n%s\n", " ".join(map(str, [script_name] + sys.argv[1:])))


def submit_mr_jobs(mtz, mr_dir, search_results, refine_type, refine_cycles, args):
    """Function to submit molecular replacement jobs

    Parameters
    ----------
    mtz : str
        Path to input mtz file
    mr_dir : str
        Path to input models directory
    mr_dir : str
        Path to directory where MR will be run
    search_results : list
        list of results from SIMBAD search
    refine_type : str
        The type of refinement to run
    refine_cycles : int
        The number of refinement cycles

    Returns
    -------
    object
        MrSubmit object containing results from MR
    """
    from simbad.mr import MrSubmit

    molecular_replacement = MrSubmit(
        mtz=mtz,
        mr_program=args.mr_program,
        refine_program=args.refine_program,
        refine_type=refine_type,
        refine_cycles=refine_cycles,
        output_dir=mr_dir,
        sgalternative=args.sgalternative,
        tmp_dir=args.tmp_dir,
        timeout=args.phaser_kill,
    )
    molecular_replacement.submit_jobs(
        search_results,
        nproc=args.nproc,
        process_all=args.process_all,
        submit_nproc=args.submit_nproc,
        submit_qtype=args.submit_qtype,
        submit_queue=args.submit_queue,
    )
    return molecular_replacement
