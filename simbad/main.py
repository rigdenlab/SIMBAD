#!/usr/bin/env ccp4-python
'''
This is SIMBAD

Sequence Independent Molecular replacement Based on Available Database

@author hlasimpk
'''

# imports
import argparse
import logging
import os
import platform
import sys
import time

from simbad.util import argparse_util
from simbad.util import config_util
from simbad.parsers import database_parser
from simbad.util import exit_util
from simbad.util import lattice_util
from simbad.util import logging_util
from simbad.util import mtz_util
from simbad.util import options_processor
from simbad.util import simbad_util
from simbad.util import version
from simbad.util import workers_util

__main_author__ = "Adam Simpkin"
__contributing_authors__ = "Jens Thomas, Felix Simkovic and Ronan Keegan"
__credits__ = "Daniel Rigden, William Shepard, Charles Ballard, Villi Uski, Andrey Lebedev"
__email__ = "hlasimpk@liv.ac.uk"
__version__ = version.__version__

LOGGER = logging_util.setup_console_logging()
monitor = None


class SIMBAD(object):
    """
    Class identify candidate structures for use in molecular replacement
    independent of sequence
    """

    def __init__(self):

        self.sopt = None
        return

    def setup(self, optd):

        # Check if work directory exists or make it
        if optd['work_dir'] and os.path.exists(optd['work_dir']):
            pass
        else:
            try:
                os.mkdir(optd['work_dir'])
            except:
                msg = "Cannot create work_dir {0}".format(optd['work_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])

        # Go to the work directory
        os.chdir(optd['work_dir'])

        # Check mandatory/exclusive options
        optd = options_processor.check_mandatory_options(optd)

        # Set up logging
        simbad_log = os.path.join(optd['work_dir'], 'SIMBAD.log')
        debug_log = os.path.join(optd['work_dir'], 'debug.log')
        optd['simbad_log'] = simbad_log

        # Set up ample output file and debug log file.
        logging_util.setup_file_logging(simbad_log, level=logging.INFO)
        logging_util.setup_file_logging(debug_log, level=logging.DEBUG)

        # Make sure the CCP4 environment is set up properly
        ccp4_home = self.setup_ccp4(optd)
        ccp4_version = ".".join([str(x) for x in optd['ccp4_version']])

        # Print out Version and invocation
        LOGGER.info(simbad_util.header)
        LOGGER.info("SIMBAD version: {0}".format(version.__version__))
        LOGGER.info("Running with CCP4 version: {0} from directory: {1}".format(ccp4_version, ccp4_home))
        LOGGER.info("Running on host: {0}".format(platform.node()))
        LOGGER.info("Running on platform: {0}".format(platform.platform()))
        LOGGER.info("Job started at: {0}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
        LOGGER.info("Invoked with command-line:\n{0}\n".format(" ".join(sys.argv)))
        LOGGER.info("Running in directory: {0}\n".format(optd['work_dir']))

        ################################################################################################################
        # SHOULD ADD IN CODE HERE TO REPORT THE OUTPUT OF SIMBAD ON THE FLY
        ################################################################################################################

        optd = options_processor.process_options(optd)

        LOGGER.info('All needed programs are found, continuing...')

        return optd

    def main(self, args=None):
        """
        Main SIMBAD routine
        :param args:
        :return:
        """

        # Set command line options
        argso = self.process_command_line(args=args)
        self.sopt = sopt = config_util.SIMBADConfigOptions()
        sopt.populate(argso)

        print self.amopt.d

        # Setup things like logging, file structure, etc...
        self.setup(sopt.d)

        # Display the parameters used
        LOGGER.debug(sopt.prettify_parameters())

        sopt.write_config_file()

        print sopt.d

        if sopt.d['lattice'] == "True":
            
            ### NEEDS TESTING ###
            os.chdir(sopt['work_dir'])
            os.mkdir('lattice_input_models')
            sopt = lattice_util.Lattice_search(sopt.d)


        exit()
        ########################################################################
        # SCRIPT PROPER STARTS HERE    

        time_start = time.time()

        os.chdir(self.working_dir)
        os.mkdir('logs')
        os.mkdir('scripts')
        os.chdir('logs')

        if self.mtz_in_file and self.pdb_database:
            # Run AMORE
            self.amore()

        # Timing data
        time_stop = time.time()
        elapsed_time = time_stop - time_start
        run_in_min = elapsed_time / 60
        run_in_hours = run_in_min / 60
        msg = os.linesep + 'All processing completed  (in {0:6.2F} hours)'.format(run_in_hours) + os.linesep
        msg += '----------------------------------------' + os.linesep
        logging.info(msg)
        msg += 'Results can be viewed in {0}'.format(self.simbad_log) + os.linesep
        sys.stdout.write(msg)

    def process_command_line(self, args=None, anomalous=True):
        """
        Process command line.

        :param args: optional argument that can hold the command-line arguments if called within python for testing
        :return:
        """
        parser = argparse.ArgumentParser(
            description="SIMBAD: Sequence Independent Molecular replacement Based on Available Database",
            prefix_chars="-")

        argparse_util.add_general_options(parser)
        argparse_util.add_cluster_submit_options(parser)

        if anomalous: argparse_util.add_anomalous_options(parser)

        return vars(parser.parse_args(args))

    ### ALL PROGRAMS FROM HERE NEEDS TO BE STORED IN SEPARATE UTIL FILES TO BE CALLED BY THE MAIN ROUTINE

    def amore(self):
        """
        Run Amore sortfun and rotfun run on a pdb database of pre-calculated
        spherical harmonics
        Note: database requires:
        pdb file ([name].pdb), 
        spherical harmonics file ([name].clmn)
        hkl file ([name].hkl)
        integration radius file ([name]_intrad.txt)
        table file ([name]_search-sfs.tab]
        molecular weight file ([name]_MW.txt)
        Formatted in the same way as the MORDA database

        :return:
        """

        # Give a header output
        sys.stdout.write("###################################################\n")
        sys.stdout.write("Running amore SORTFUN run...\n")
        sys.stdout.write("###################################################\n")
        sys.stdout.write("\n")

        self.sortfun()

        # Give a header output
        sys.stdout.write("###################################################\n")
        sys.stdout.write("Running amore ROTFUN run...\n")
        sys.stdout.write("###################################################\n")
        sys.stdout.write("\n")

        # Create a dictionary of all the files in the database
        pdb_dictionary = database_parser.directory_information(self.pdb_database)

        sys.stdout.write('Generating scripts')
        # Run amore on each PDB in the database in parallel
        script_list = []
        pdb_names = database_parser.get_pdb_names(self.pdb_database)
        for name in pdb_names:
            info = pdb_dictionary[name]
            tab = None
            hkl = None
            clmn = None
            intrad = None
            for f in info:
                database_file = os.path.join(self.pdb_database, name[1:3], f)
                if 'search-sfs.tab' in database_file:
                    tab = database_file
                elif 'search.hkl' in database_file:
                    hkl = database_file
                elif 'search.clmn' in database_file:
                    clmn = database_file
                elif 'intrad' in database_file:
                    intrad = database_file
                elif "MW.txt" in database_file:
                    molecular_weight = database_file
                else:
                    continue
                if self.matt_coef(name, molecular_weight):
                    script_list.append(self.rotfun(tab, hkl, clmn, intrad, name))
                else:
                    continue

        os.mkdir('rotfun_output')

        sys.stdout.write('Running scripts')
        # Run the jobs
        workers_util.run_scripts(job_scripts=script_list,
                                 monitor=None,
                                 chdir=False,
                                 nproc=self.nproc,
                                 job_time=7200,
                                 job_name='SIMBAD',
                                 submit_cluster=self.submit_cluster,
                                 submit_qtype=self.submit_qtype,
                                 submit_queue=self.submit_queue,
                                 submit_array=self.submit_array)

        return

    def sortfun(self):

        # Need a way to find FP and SIGFP names
        cmd = ["{0}  hklin {1}".format(self.amore_installation, self.mtz_in_file),
               "hklpck0 {0}".format(os.path.join(self.working_dir, "spmipch.hkl"))]
        command_line = os.linesep.join(map(str, cmd))
        key = """TITLE   ** spmi  packing h k l F for crystal**
        SORTFUN RESOL 100.  2.5
        LABI FP={0}  SIGFP={1}
        """.format(self.fp, self.sigf)

        logfile = os.path.join(self.working_dir, 'SORTFUN.log')
        simbad_util.run_job(command_line, logfile, key)

    def matt_coef(self, name, molecular_weight):

        cmd = ["matthews_coef"]
        key = """CELL {0}
symm {1}
molweight {2}
auto""".format(self.cell_parameters, self.space_group, molecular_weight)
        logfile = os.path.join(self.working_dir, 'matt_coef_{0}.log'.format(name))
        simbad_util.run_job(cmd, logfile, key)

        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith('  1'):
                    solvent_content = float(line.split()[2])
                    if solvent_content >= 30:
                        result = True
                    else:
                        result = False

        os.remove(logfile)

        return result

    def rotfun(self, tab, hkl, clmn, intrad, name):
        # Make a directory for each PDB
        run_dir = os.path.join(self.working_dir, 'rotfun_output', name)
        os.mkdir(run_dir)

        cmd = ["{0}  table1 {1}".format(self.amore_installation, tab),
               "       HKLPCK1 {0}".format(hkl),
               "       hklpck0 {0}".format(os.path.join(self.working_dir, "spmipch.hkl")),
               "       clmn1 {0}".format(clmn),
               "       clmn0 {0}".format(os.path.join(run_dir, "spmipch.clmn")),
               "       MAPOUT {0}".format(os.path.join(run_dir, "amore_cross.map"))]
        command_line = os.linesep.join(map(str, cmd))
        stdin_txt = """ROTFUN
VERB
TITLE : Generate HKLPCK1 from MODEL FRAGMENT   1
CLMN CRYSTAL ORTH  1 RESO  20.0  {0}  SPHERE   {1}
ROTA  CROSS  MODEL 1  PKLIM {3}  NPIC {4}""".format(self.SHRES, intrad, self.PKLIM, self.NPIC)

        stdin_txt = """{0} << EOF
{1}
EOF
""".format(command_line, stdin_txt)
        script_path = os.path.join(self.working_dir, 'scripts', "{0}.sh".format(name))
        with open(script_path, "w") as job_script:
            # Header
            script_header = '#!/bin/sh\n'
            job_script.write(script_header)
            job_script.write(stdin_txt)

        # Make executable
        os.chmod(script_path, 0o777)

        return script_path

    def setup_ccp4(self, amoptd):
        # type: (object) -> object
        """Check CCP4 is available and return the top CCP4 directory"""
        # Make sure CCP4 is around
        if not "CCP4" in os.environ:
            msg = "Cannot find CCP4 installation - please make sure CCP4 is installed and the setup scripts have been run!"
            exit_util.exit_error(msg)

        if not "CCP4_SCR" in os.environ:
            msg = "$CCP4_SCR environement variable not set - please make sure CCP4 is installed and the setup scripts have been run!"
            exit_util.exit_error(msg)

        if not os.path.isdir(os.environ['CCP4_SCR']):
            msg = "*** WARNING ***\n"
            msg += "Cannot find the $CCP4_SCR directory: {0}\n".format(os.environ['CCP4_SCR'])
            msg += "The directory will be created, but it should have already been created by the CCP4 startup scripts\n"
            msg += "Please make sure CCP4 is installed and the setup scripts have been run."
            LOGGER.critical(msg)
            os.mkdir(os.environ['CCP4_SCR'])
            # exit_util.exit_error(msg)

        # Record the CCP4 version we're running with  - also required in pyrvapi_results
        amoptd['ccp4_version'] = simbad_util.ccp4_version()

        return os.environ['CCP4']


if __name__ == "__main__":
    try:
        SIMBAD().main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
