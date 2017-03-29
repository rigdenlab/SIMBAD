#!/usr/bin/env ccp4-python
'''
This is SIMBAD

Sequence Independent Molecular replacement Based on Available Database

@author hlasimpk
'''


import argparse
import logging
import os
import platform
import sys
import time

from constants import CONTAMINANT_MODELS, SIMBAD_LATTICE_DB
from simbad.rotsearch import amore_search
from simbad.util import argparse_util
from simbad.util import config_util
from simbad.util import exit_util
from simbad.util import logging_util
from simbad.util import mr_util
from simbad.util import options_processor
from simbad.util import simbad_util
from simbad.util import version

import simbad.lattice.search

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

        self.finished = False
        self.sopt = None

        return

    def setup(self, optd):

        if optd['work_dir']:
            LOGGER.info('Making a named work directory: {0}'.format(optd['work_dir']))
            try:
                os.mkdir(optd['work_dir'])
            except:
                msg = "Cannot create work_dir {0}".format(optd['work_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])
        else:
            if not os.path.exists(optd['run_dir']):
                msg = "Cannot find run directory: {0}".format(optd['run_dir'])
                exit_util.exit_error(msg, sys.exc_info()[2])
            LOGGER.info('Making a run_directory: checking for previous runs...')
            optd['work_dir'] = simbad_util.make_workdir(optd['run_dir'],
                                                        ccp4_jobid=optd['ccp4_jobid'])

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

        ########################################################################
        # SCRIPT PROPER STARTS HERE
        ########################################################################

        time_start = time.time()

        # Set command line options
        argso = self.process_command_line(args=args)
        self.sopt = sopt = config_util.SIMBADConfigOptions()
        sopt.populate(argso)

        # Setup things like logging, file structure, etc...
        self.setup(sopt.d)

        # Display the parameters used
        LOGGER.debug(sopt.prettify_parameters())

        sopt.write_config_file()

        sopt.d['solution'] = False

        ########################################################################
        # Lattice search
        ########################################################################

        if sopt.d['lattice'] == "True":
            # Perform the lattice search
            lattice_search = simbad.lattice.search.LatticeSearch(
                sopt.d['cell_parameters'], sopt.d['space_group'],
                SIMBAD_LATTICE_DB, max_to_keep=sopt.d['njob_lattice']
            )
            lattice_search.search()

            # Only create lattice directories if lattice search produced results
            search_results = lattice_search.search_results
            if search_results:
                lattice_search.summarize()
                if sopt.d['pdb_db']:
                    lattice_search.copy_results(sopt.d['pdb_db'], sopt.d['work_dir'])
                else:
                    lattice_search.download_results(sopt.d['work_dir'])

                # Create directories for lattice search MR
                os.mkdir('MR_LATTICE')

                # Run MR on results
                molecular_replacement = mr_util.MrSubmit(sopt.d['mtz'], sopt.d['MR_program'], sopt.d['refine_program'],
                                                         os.path.join(sopt.d['work_dir'], 'lattice_input_models'),
                                                         os.path.join(sopt.d['work_dir'], 'MR_LATTICE'))
                molecular_replacement.multiprocessing(search_results, nproc=sopt.d['nproc'])
                molecular_replacement.summarize()

                # Check if a solution was found
                for model in search_results:
                    if not sopt.d['solution']:
                        try:
                            sopt.d['solution'] = molecular_replacement.solution_found(model)
                        except:
                            pass

            else:
                LOGGER.info("No results found - lattice search was unsuccessful")

        elif sopt.d["lattice"] != "True":
            LOGGER.info("Lattice run set to {0}: Skipping...".format(sopt.d["lattice"]))

        if sopt.d['solution'] and sopt.d['early_term'] or sopt.d['contaminant'] == 'False' and sopt.d['full'] == 'False':
            self.finished = True

        ########################################################################
        # Contaminant search
        ########################################################################

        if sopt.d['contaminant'] == "True" and not (sopt.d['early_term'] and sopt.d['solution']):
            # Create work directories
            os.mkdir(os.path.join(sopt.d['work_dir'], 'output'))
            os.mkdir(os.path.join(sopt.d['work_dir'], 'clogs'))

            rotation_search = amore_search.AmoreRotationSearch(os.path.join(os.environ["CCP4"], 'bin', 'amoreCCB2.exe'),
                                                               sopt.d['mtz'], sopt.d['work_dir'], sopt.d['njob_contam'])
            rotation_search.sortfun()
            rotation_search.amore_run(CONTAMINANT_MODELS, os.path.join(sopt.d['work_dir'], 'clogs'), sopt.d['nproc'],
                                      sopt.d['SHRES'], sopt.d['PKLIM'], sopt.d['NPIC'], sopt.d['ROTASTEP'],
                                      sopt.d['min_solvent_content'])
            rotation_search.summarize()
            search_results = rotation_search.search_results

            if search_results:
                # Create directories for the contaminant search MR
                os.mkdir('MR_CONTAMINANT')

                # Run MR on results
                molecular_replacement = mr_util.MrSubmit(sopt.d['mtz'], sopt.d['MR_program'], sopt.d['refine_program'],
                                                         CONTAMINANT_MODELS,
                                                         os.path.join(sopt.d['work_dir'], 'MR_CONTAMINANT'))
                molecular_replacement.multiprocessing(search_results, nproc=sopt.d['nproc'])
                molecular_replacement.summarize()

                # Check if a solution was found
                for model in search_results:
                    if not sopt.d['solution']:
                        try:
                            sopt.d['solution'] = molecular_replacement.solution_found(model)
                        except:
                            pass
            else:
                LOGGER.info("No results found - Contaminant search was unsuccessful")

        elif sopt.d["contaminant"] != "True":
            LOGGER.info("Contaminant run set to {0}: Skipping...".format(sopt.d["contaminant"]))

        if sopt.d['solution'] and sopt.d['early_term'] or sopt.d['full'] == 'False':
            self.finished = True

        ########################################################################
        # Full search
        ########################################################################

        if sopt.d['full'] == True and not (sopt.d['early_term'] and sopt.d['solution']):
            # Create work directories
            os.mkdir(os.path.join(sopt.d['work_dir'], 'output'))
            os.mkdir(os.path.join(sopt.d['work_dir'], 'logs'))
            sopt.d['mode'] = "FULL_ROT"


            first = True
            count = 0
            if sopt.d['morda_db']:
                for e in os.walk(sopt.d['morda_db']):
                    if first:
                        first = False
                        pass
                    else:
                        if count < 1:
                            pdb_dir = e[0]
                            count += 1


            elif sopt.d['sphere_database']:
                pass



            self.finished = True

        elif sopt.d["full"] != "True" and not self.finished:
            LOGGER.info("Full run set to {0}: Skipping...".format(sopt.d["full"]))


        ########################################################################
        # Finish up
        ########################################################################

        if self.finished:

            #  Add a clean up function here to remove files from output

            # Timing data
            time_stop = time.time()
            elapsed_time = time_stop - time_start
            run_in_min = elapsed_time / 60
            run_in_hours = run_in_min / 60
            msg = os.linesep + 'All processing completed  (in {0:6.2F} hours)'.format(run_in_hours) + os.linesep
            msg += '----------------------------------------' + os.linesep
            msg += 'Results can be viewed in {0}'.format(os.path.join(sopt.d['work_dir'], 'SIMBAD.log')) + os.linesep
            LOGGER.info(msg)


            exit()

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
