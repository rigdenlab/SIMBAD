"""Class to run MR on SIMBAD results using code from MrBump"""

__author__ = "Adam Simpkin"
__date__ = "09 Mar 2017"
__version__ = "0.1"

import os
import sys
import logging
import tempfile
import types
import copy_reg
from contextlib import contextmanager
import multiprocessing

from simbad.util import simbad_util
from simbad.util import mtz_util

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

    def __init__(self, mtz, mr_program, refine_program, model_dir, work_dir, early_term=True):
        self._early_term = None
        self._model_dir = None
        self._mtz = None
        self._mr_program = None
        self._refine_program = None
        self._work_dir = None

        self.early_term = early_term
        self.model_dir = model_dir
        self.mtz = mtz
        self.mr_program = mr_program
        self.refine_program = refine_program
        self.work_dir = work_dir

        # Set by the program for now, may change to use properties
        self.input_file = None

        # options derived from the input mtz
        self.cell_parameters = None
        self.resolution = None
        self.solvent = None
        self.space_group = None
        self.f = None
        self.sigf = None
        self.free = None

    @property
    def early_term(self):
        """Flag to decide if the program should terminate early"""
        return self._early_term

    @early_term.setter
    def early_term(self, early_term):
        """Set the early term flag to true or false"""
        self._early_term = early_term

    @property
    def model_dir(self):
        """The directory containing the input models"""
        return self._model_dir

    @model_dir.setter
    def model_dir(self, model_dir):
        """Define the path to the directory containing the input models"""
        self._model_dir = model_dir

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = mtz

    @property
    def mr_program(self):
        """The molecular replacement program to use"""
        return self._mr_program

    @mr_program.setter
    def mr_program(self, mr_program):
        """Define the molecular replacement program to use"""
        self._mr_program = mr_program

    @property
    def refine_program(self):
        """The refinement program to use"""
        return self._refine_program

    @refine_program.setter
    def refine_program(self, refine_program):
        """Define the refinement program to use"""
        self._refine_program = refine_program

    @property
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    @contextmanager
    def suppress_stdout(self):
        """Code to suppress the stdout of a function, used here to silence stdout of MrBUMP modules"""
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                yield
            finally:
                sys.stdout = old_stdout

    def run_job(self, mtz, model, output_dir, enam=False, early_term=True):
        _logger.info("Running MR and refinement on {0}".format(model.pdb_code))

        # parse options for the model
        self.parse_options(model)

        # Get information about the input mtz
        self.get_mtz_info(mtz)

        # Generate MR input file
        self.MR_setup(model, output_dir, mtz, enam)

        # Run job
        cljob = cluster_run.ClusterJob()

        if self.input_file:
            # Parse input file
            cljob.parse_input(self.input_file)

        # Create temp files for the mr and refine keyfiles
        mr_key = tempfile.NamedTemporaryFile(delete=False)
        ref_key = tempfile.NamedTemporaryFile(delete=False)

        if self.mr_program and self.refine_program:
            with self.suppress_stdout():
                cljob.run(self.mr_program, self.refine_program, mr_key.name, refine_keyfile=ref_key.name)

        # Remove temp files
        os.unlink(mr_key.name)
        os.unlink(ref_key.name)

        if early_term:
            terminate = self.solution_found(model)
            return terminate

        return False

    def get_mtz_info(self, mtz):
        # Extract crystal data from input mtz
        self.cell_parameters, self.resolution, self.space_group = mtz_util.crystal_data(mtz)

        # Extract column labels from input mtz
        self.f, self.sigf, _, _, self.free = mtz_util.get_labels(mtz)

        # Get solvent content
        self.solvent = self.matthews_coef(self.cell_parameters, self.space_group)
        return

    def MR_setup(self, model, output_dir, mtz, enam):
        """Code to generate directories for MR and refinement for a model
        and generate an input file containing information about the MR/refine run"""

        # Create individual directories for every results
        if self.mr_program.upper() == "MOLREP":
            os.chdir(self.work_dir)
            os.mkdir(os.path.join(output_dir, model.pdb_code))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'molrep'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'molrep', 'refine'))
        elif self.mr_program.upper() == "PHASER":
            os.chdir(self.work_dir)
            os.mkdir(os.path.join(output_dir, model.pdb_code))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'phaser'))
            os.mkdir(os.path.join(output_dir, model.pdb_code, 'mr', 'phaser', 'refine'))

        self.generate_mr_input_file(model.pdb_code, output_dir, mtz, enam)

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

    def generate_mr_input_file(self, model, output_dir, mtz, enam=False):
        """Create an input file for MR"""

        # create input file path
        input_file = os.path.join(self.work_dir, output_dir, model, 'input.txt')
        dire = os.path.join(self.work_dir, output_dir, model)
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
            f.write("SOLV {0}\n".format(self.solvent))
            f.write("RESO {0}\n".format(self.resolution))
            f.write("PDBI {0}\n".format(pdbi))
            if self.mr_program == "phaser":
                f.write("HKLO {0}\n".format(hklo))

        self.input_file = input_file

        return

    def matthews_coef(self, cell_parameters, space_group):
        """Function to run matthews coefficient to decide if the model can fit in the unit cell

        Parameters
        ----------
        cell_parameters
            The parameters of the unit cell
        space_group
            The space group of the crystal

        Returns
        -------
        float
            solvent content of the protein

        """

        cmd = ["matthews_coef"]
        key = """CELL {0}
symm {1}
auto""".format(cell_parameters,
               space_group)

        logfile = os.path.join(self.work_dir, 'matt_coef.log')
        simbad_util.run_job(cmd, logfile, key)

        # Determine if the model can fit in the unit cell
        solvent_content = 0.5
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith('  1'):
                    solvent_content = (float(line.split()[2]) / 100)

        # Clean up
        os.remove(logfile)

        return solvent_content