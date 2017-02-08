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
import sys
import time

from simbad.util import config_util
from simbad.parsers import database_parser
from simbad.util import exit_util
from simbad.util import mtz_util
from simbad.util import simbad_util
from simbad.util import workers_util

def setup_console_logging():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # First create console logger for outputting stuff
    # create file handler and set level to debug
    # Seems they changed the api in python 2.6->2.7
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(strm=sys.stdout)
    cl.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s\n') # Always add a blank line after every print
    cl.setFormatter(formatter)
    logger.addHandler(cl)
    return logger

def setup_file_logging(main_logfile, debug_logfile):
    """
    Set up the various log files/console logging and return the logger

    :param main_logfile:
    :param debug_logfile:
    :return:
    """

    logger = logging.getLogger()

    # create file handler for debug output
    fl = logging.FileHandler(debug_logfile)
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s [%(lineno)d] - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)

    # Finally create the main logger
    fl = logging.FileHandler(main_logfile)
    fl.setLevel(logging.INFO)
    fl.setFormatter(formatter) # Same formatter as screen
    logger.addHandler(fl)
    
    return logger

LOGGER = setup_console_logging()
monitor = None

class ArgParser(object):
    """
    Class to add command line arguments
    """

    def main(self, args=None):
        parser = argparse.ArgumentParser(description="SIMBAD: Sequence Independent Molecular replacement Based on Available Database", 
                                         prefix_chars="-")
        self._add_GENERAL(parser)
        
    def _add_GENERAL(self, parser):
        parser.add_argument('-mtz', metavar='MTZ in', type=str, nargs=1,
                           help='The MTZ file with the reflection data.')
        
        parser.add_argument('-nproc', type=int, nargs=1,
                           help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors. " + \
                            "For cluster submission, this should be the number of processors on a node.")
        
        parser.add_argument('-submit_array', metavar='True/False', type=str, nargs=1,
                           help='Submit SGE jobs as array jobs')
        
        parser.add_argument('-submit_cluster', metavar='True/False', type=str, nargs=1,
                           help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')
        
        parser.add_argument('-submit_qtype', type=str, nargs=1,
                           help='cluster submission queue type - currently support SGE and LSF')
        
        parser.add_argument('-submit_queue', type=str, nargs=1,
                           help='The queue to submit to on the cluster.')
        
        parser.add_argument('-work_dir', type=str, nargs=1,
                           help='Path to the directory where SIMBAD will run (will be created if it doesn\'t exist)')

        parser.add_argument('-url', type=str, nargs=1,
                           help='URL for use on SIMBAD server')

        args = vars(parser.parse_args())
        return args
    

class SIMBAD(object):
    """
    Class identify candidate structures for use in molecular replacement
    independent of sequence
    """
    
    def __init__(self):
        self.amore_installation = 'amore'
        self.mtz_in_file = None
        self.nproc = 1
        self.pdb_database = None
        self.working_dir = None
        self.anomalous_signal = False
        
        # matt coef options
        self.space_group = None
        self.resolution = None
        self.cell_parameters = None
        
        # AMORE options
        self.SHRES=3.0
        self.PKLIM=0.5
        self.NPIC=50
        self.ROTASTEP=1.0

        
        # Cluster options
        self.submit_array = None
        self.submit_cluster = None
        self.submit_qtype = None
        self.submit_queue = None
        
        self.simbad_log = None
        return
    
    def process_command_line(self, args):
        """
        Process command line arguments and
        store dictionary entries as class variables

        :param args:
        :return:
        """

        optd = ArgParser().main(args)
        
        if 'mtz' in optd.keys() and optd['mtz'] and os.path.exists(optd['mtz']):
            self.mtz_in_file = os.path.abspath(optd['mtz'])
        else:
            raise RuntimeError("MTZ not defined")
        
        if 'work_dir' in optd.keys() and optd['work_dir'] and os.path.exists(optd['work_dir']):
            self.working_dir = os.path.relpath(optd['work_dir'])
        else:
            try:
                os.mkdir(optd['work_dir'])
            except:
                msg = "Cannot create work_dir {0}".format(optd['work_dir'])
                raise RuntimeError(msg)
            self.working_dir = os.path.relpath(optd['work_dir'])
        
        mtz_util.processReflectionFile(optd)
        if optd['FP'] != None:
            self.fp = optd['FP']
        if optd['SIGF'] != None:
            self.sigf = optd['SIGF']
        if optd['DANO'] != None and optd['SIGDANO'] != None:
            self.anomalous_signal = True
        
        self.space_group, self.resolution, self.cell_parameters = mtz_util.set_crystal_data(self.mtz_in_file)
        
        return
    
    def main(self, args=None):
        """
        Main SIMBAD routine
        :param args:
        :return:
        """

        # Set command line options
        argso = self.process_command_line(args=args)
        self.amopt = amopt = config_util.SIMBADConfigOptions()
        amopt.populate(argso)
        
        # Set up things such as logging
        self.setup()
        
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
        #Run the jobs
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
        command_line = os.linesep.join(map(str,cmd))
        key="""TITLE   ** spmi  packing h k l F for crystal**
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
        command_line = os.linesep.join(map(str,cmd))
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
    
    def setup(self):
        
        # Set up logging
        simbad_log = os.path.join(self.working_dir, 'SIMBAD.log')
        debug_log = os.path.join(self.working_dir, 'debug.log')
        self.simbad_log = simbad_log
        
        setup_file_logging(simbad_log, debug_log)
    
    def setup_ccp4(self, amoptd):
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
            #exit_util.exit_error(msg)
    
        # Record the CCP4 version we're running with  - also required in pyrvapi_results
        amoptd['ccp4_version'] = simbad_util.ccp4_version()
        
        return os.environ['CCP4']
    
if __name__ == "__main__":
    try:
        SIMBAD().main()
    except Exception as e:
        msg = "Error running main SIMBAD program: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
        
    
    