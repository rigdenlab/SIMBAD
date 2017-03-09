"""Class to run MR on SIMBAD results using code from MrBump"""

__author__ = "Adam Simpkin"
__date__ = "09 Mar 2017"
__version__ = "0.1"

import os
import sys
import logging
import types
import copy_reg
from contextlib import contextmanager
import multiprocessing

from simbad.util import simbad_util
from simbad.constants import SIMBAD_EGG_ROOT

# Set up MrBUMP imports
if os.path.isdir(os.path.join(os.environ["CCP4"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
mrbump_incl = os.path.join(mrbump, "include")

sys.path.append(os.path.join(mrbump_incl, 'dev'))
sys.path.append(os.path.join(mrbump_incl, 'file_info'))
sys.path.append(os.path.join(mrbump_incl, 'modelling'))
sys.path.append(os.path.join(mrbump_incl, 'mr'))
sys.path.append(os.path.join(mrbump_incl, 'output'))
sys.path.append(os.path.join(mrbump_incl, 'seq_align'))
sys.path.append(os.path.join(mrbump_incl, 'ccp4'))
sys.path.append(os.path.join(mrbump_incl, 'structures'))
sys.path.append(os.path.join(mrbump_incl, 'tools'))
sys.path.append(os.path.join(mrbump_incl, 'initialisation'))
sys.path.append(os.path.join(mrbump_incl, 'building'))
sys.path.append(os.path.join(mrbump_incl, 'cluster'))
sys.path.append(os.path.join(mrbump_incl, 'dispatchers'))
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
sys.path.append(os.path.join(mrbump_incl, 'phasing'))

import cluster_run

_logger = logging.getLogger(__name__)

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

class MrSubmit(object):
    """Class to run MR on a defined set of models"""

    def __init__(self, mr_program, refine_program, work_dir, space_group, fp, sigf, free, resolution):
        self.input_file = None
        self.mr_program = None
        self.refine_program = None
        self.work_dir = None

        self.space_group = space_group
        self.fp = fp
        self.sigf = sigf
        self.free = free
        self.resolution = resolution


    def parse_options(self, model):
        """Function to set up input options for the job"""

        self.input_file = os.path.join(self.optd.d['work_dir'], self.optd.d['mode'], model.pdb_code, 'input.txt')

        # Set path to MR keyword file
        self.mr_program = self.optd.d['MR_program']
        if self.mr_program.upper() == 'MOLREP' and not self.optd.d['mr_keywords']:
            self.mr_keyfile = os.path.join(SIMBAD_EGG_ROOT, 'static', 'default_molrep_keywords.txt')
        elif self.mr_program.upper() == 'PHASER' and not self.optd.d['mr_keywords']:
            self.mr_keyfile = os.path.join(SIMBAD_EGG_ROOT, 'static', 'default_phaser_keywords.txt')
        elif self.optd.d['mr_keywords']:
            self.mr_keyfile = self.optd.d['mr_keywords']
        else:
            msg = "Unable to find MR keyword file for MR program: {0}".format(self.optd.d['MR_program'])
            _logger.debug(msg)
            raise RuntimeError(msg)

        # Set path to refinement keyword file
        self.refine_program = self.optd.d['refine_program']
        if self.refine_program.upper() == "REFMAC5" and not self.optd.d['refine_keywords']:
            self.refine_keyfile = os.path.join(SIMBAD_EGG_ROOT, 'static', 'default_refmac_keywords.txt')
        elif self.optd.d['refine_keywords']:
            self.refine_keyfile = self.optd.d['refine_keywords']
        else:
            msg = "Unable to find refinement keyword file for refinement program: {0}".format(self.optd.d['refine_program'])
            _logger.debug(msg)
            raise RuntimeError(msg)
        return

    def run_job(self, model):
        _logger.info("Running MR and refinement on {0}".format(model.pdb_code))

        # parse options for the model
        self.parse_options(model)

        # Generate MR input files
        self.MR_setup(model)

        # Run job
        cljob = cluster_run.ClusterJob()

        if self.input_file:
            # Parse input file
            cljob.parse_input(self.input_file)

        if self.mr_keyfile and self.mr_program and self.refine_keyfile and self.refine_program:
            with suppress_stdout():
                cljob.run(self.mr_program, self.refine_program, self.mr_keyfile, refine_keyfile=self.refine_keyfile)

        if self.optd.d['early_term']:
            terminate = self.solution_found(model)
            return terminate

        return False

    def MR_setup(self, model, output_dir, work_dir):
        """Code to generate directories for MR and refinement for a model
        and generate an input file containing information about the MR/refine run"""

        # Create individual directories for every results
        if self.optd.d['MR_program'].upper() == "MOLREP":
            os.chdir(work_dir)
            os.mkdir(os.path.join(output_dir, model.pdb_code))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'molrep'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'molrep', 'refine'))
        elif self.optd.d['MR_program'].upper() == "PHASER":
            os.chdir(work_dir)
            os.mkdir(os.path.join(output_dir, model.pdb_code))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'phaser'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'phaser', 'refine'))

        simbad_util.generate_mr_input_file(self.optd, model.pdb_code, 'lattice')

        return

    def multiprocessing(self, results):
        """Code to run MR and refinement in parallel"""

        def run(job_queue):
            """processes element of job queue if queue not empty"""
            TIME_OUT_IN_SECONDS = 60

            while not job_queue.empty():
                model = job_queue.get(timeout=TIME_OUT_IN_SECONDS)
                terminate = self.run_job(model)

                if terminate:
                    print "MR with {0} was successful so removing remaining jobs from inqueue".format(model.pdb_code)
                    while not job_queue.empty():
                        job = job_queue.get()
                        _logger.debug("Removed job [{0}] from inqueue".format(job))



        # Create job queue
        job_queue = multiprocessing.Queue()

        # Add each result from results to the job queue
        for result in results:
            job_queue.put(result)

        processes = []
        # Set up processes equal to the number of processors input
        for i in range(self.optd.d['nproc']):
            process = multiprocessing.Process(target=run, args=(job_queue,))
            process.start()
            processes.append(process)

        # block the calling thread
        for p in processes:
            p.join()


    def solution_found(self, result_dir, model):
        """Function to check is a solution has been found"""
        # Set default values if no results found
        final_r_fact = 1
        final_r_free = 1

        with open(os.path.join(result_dir, model.pdb_code, 'mr', 'molrep',
                               'refine', model.pdb_code + '_ref.log'), 'r') as f:
            for line in f:
                if line.startswith('           R factor'):
                    final_r_fact = float(line.split()[2])
                elif line.startswith('             R free'):
                    final_r_free = float(line.split()[2])

        if final_r_fact < 0.45 and final_r_free < 0.45:
            return True
        else:
            return False

    def generate_mr_input_file(self, model, output_dir, work_dir, mtz, solvent, enam=False):
        """Create an input file for MR"""

        # create input file path
        input_file = os.path.join(work_dir, output_dir, model, 'input.txt')
        dire = os.path.join(work_dir, output_dir, model)
        pdbi = os.path.join(self.model_dir, '{0}.pdb'.format(model))

        # Assign variables
        hklr = os.path.join(dire, 'mr', self.mr_program, '{0}_refinement_input.mtz'.format(model))
        hklo = os.path.join(dire, 'mr', self.mr_program, '{0}_phaser_output.mtz'.format(model))
        pdbo = os.path.join(dire, 'mr', self.mr_program, '{0}_mr_output.pdb'.format(model))
        mrlo = os.path.join(dire, 'mr', self.mr_program, '{0}_mr.log'.format(model))
        refh = os.path.join(dire, 'mr', self.mr_program, 'refine', '{0}_refinement_output.mtz'.format(model))
        refp = os.path.join(dire, 'mr', self.mr_program, 'refine', '{0}_refinement_output.pdb'.format(model))
        refl = os.path.join(dire, 'mr', self.mr_program, 'refine', '{0}_ref.log'.format(model))

        # Write input file
        with open(input_file, 'w') as f:
            f.write("DIRE {0}\n".format(dire))
            f.write("SGIN {0}\n".format(self.space_group))
            f.write("HKL1 {0}\n".format(mtz))
            f.write("HKLR {0}\n".format(hklr))
            f.write("PDBO {0}\n".format(pdbo))
            f.write("MRLO {0}\n".format(mrlo))
            f.write("REFH {0}\n".format(refh))
            f.write("REFP {0}\n".format(refp))
            f.write("REFL {0}\n".format(refl))
            f.write("ENAN {0}\n".format(enam))
            f.write("FPIN {0}\n".format(self.fp))
            f.write("SIGF {0}\n".format(self.sigf))
            f.write("FREE {0}\n".format(self.free))
            f.write("SOLV {0}\n".format(solvent))
            f.write("RESO {0}\n".format(self.resolution))
            f.write("PDBI {0}\n".format(pdbi))
            if self.mr_program == "phaser":
                f.write("HKLO {0}\n".format(hklo))

        return

    def matthews_coef(self, mtz):
        """Function to run matthews coefficient to decide if the model can fit in the unit cell

        Parameters
        ----------
        model : str
            Path to input model
        min_solvent_content : int float
            Minimum solvent content [default: 30]

        Returns
        -------
        bool
            Can the model fit in the unit cell with a solvent content higher than the min_solvent_content

        """

        cmd = ["matthews_coef"]
        key = """CELL {0}
        symm {1}
        molweight {2}
        auto""".format(self.optd.d['cell_paramaters'],
                       self.optd.d['space_group'],
                       molecular_weight)
        name = os.basename(model)[0:3]
        logfile = os.path.join(self.optd.d['work_dir'], 'matt_coef_{0}.log'.format(name))
        simbad_util.run_job(cmd, logfile, key)

        # Determine if the model can fit in the unit cell
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith('  1'):
                    solvent_content = float(line.split()[2])
                    if solvent_content >= min_solvent_content:
                        result = True
                    else:
                        result = False

        # Clean up
        os.remove(logfile)

        return result


@contextmanager
def suppress_stdout():
    """Code to suppress the stdout of a function, used here to silence stdout of MrBUMP modules"""
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

