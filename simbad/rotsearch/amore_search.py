"""Module to run the AMORE rotation search"""

from __future__ import print_function

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "07 Mar 2017"
__version__ = "0.1"

import copy_reg
import iotbx.pdb
import iotbx.pdb.mining
import logging
import multiprocessing
import numpy
import os
import pandas
import types
import warnings

from simbad.util import simbad_util
from simbad.util import mtz_util

logger = logging.getLogger(__name__)

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

class _AmoreRotationScore(object):
    """A amore rotation scoring class"""

    __slots__ = ("pdb_code", "ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
                 "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak")

    def __init__(self, pdb_code, ALPHA, BETA, GAMMA, CC_F, RF_F, CC_I, CC_P, Icp,
                 CC_F_Z_score, CC_P_Z_score, Number_of_rotation_searches_producing_peak):
        self.pdb_code = pdb_code
        self.ALPHA = ALPHA
        self.BETA = BETA
        self.GAMMA = GAMMA
        self.CC_F = CC_F
        self.RF_F = RF_F
        self.CC_I = CC_I
        self.CC_P = CC_P
        self.Icp = Icp
        self.CC_F_Z_score = CC_F_Z_score
        self.CC_P_Z_score = CC_P_Z_score
        self.Number_of_rotation_searches_producing_peak = Number_of_rotation_searches_producing_peak

    def __repr__(self):
        return "{0}(pdb_code={1} ALPHA={2} BETA={3} GAMMA={4} CC_F={5} RF_F={6} CC_I={7} CC_P={8} Icp={9} " \
               "CC_F_Z_score={10} CC_P_Z_score={11} Number_of_rotation_searches_producing_peak={12}".format(
            self.__class__.__name__, self.pdb_code, self.ALPHA, self.BETA, self.GAMMA, self.CC_F, self.RF_F,
            self.CC_I, self.CC_P, self.Icp, self.CC_F_Z_score, self.CC_P_Z_score,
            self.Number_of_rotation_searches_producing_peak
        )

    def _as_dict(self):
        """Convert the :obj:`_AmoreRotationScore <simbad.rotsearch.amore_search._AmoreRotationScore>`
        object to a dictionary"""
        dict = {}
        for k in self.__slots__:
            dict[k] = getattr(self, k)
        return dict

class AmoreRotationSearch(object):
    """A class to perform the amore rotation search

    Attributes
    ----------
    amore_exe : str
        The path to the amore executable
    mtz : str
        The path to the input MTZ
    work_dir : str
        The path to the working directory
    max_to_keep : int
        The maximum number of results to keep [default: 20]
    models_dir : str
        The directory containing the models to run the rotation search on
    logs_dir : str
        The directory where logs from the job will be placed
    nproc : int
        The number of processors to run the job on
    shres : int float
        Spherical harmonic resolution [default 3.0]
    pklim : int float
        Peak limit, output all peaks above <float> [default: 0.5]
    npic : int float
        Number of peaks to output from the translation function map for each orientation [default: 50]
    rotastep : int float
        Size of rotation step [default : 1.0]

    Examples
    --------
    >>> from simbad.rotsearch.amore_search import AmoreRotationSearch
    >>> rotation_search = AmoreRotationSearch('<amore_exe>', '<mtz>', '<work_dir>', '<max_to_keep>')
    >>> rotation_search.sortfun()
    >>> rotation_search.amore_run('<models_dir>', '<logs_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>')
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results

    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from amore.

    """

    def __init__(self, amore_exe, mtz, work_dir, max_to_keep=20):
        """Initialise a new amore rotation search class

        Parameters
        ----------
        amore_exe : str
            The path to the amore executable
        mtz : str
            The path to the input MTZ
        work_dir : str
            The path to the working directory
        max_to_keep : int
            The maximum number of results to keep [default: 20]
        """

        self._amore_exe = None
        self._max_to_keep = 0
        self._mtz = None
        self._search_results = None
        self._work_dir = None

        self.amore_exe = amore_exe
        self.max_to_keep = max_to_keep
        self.mtz = mtz
        self.work_dir = work_dir

    @property
    def amore_exe(self):
        """The amore exectutable"""
        return self._amore_exe

    @amore_exe.setter
    def amore_exe(self, amore_exe):
        """Define the amore executable"""
        self._amore_exe = amore_exe

    @property
    def max_to_keep(self):
        """The maximum number of results to keep"""
        return self._max_to_keep

    @max_to_keep.setter
    def max_to_keep(self, max_to_keep):
        """Define the maximum number of results to keep"""
        self._max_to_keep = max_to_keep

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = mtz

    @property
    def search_results(self):
        """The results from the amore rotation search"""
        return sorted(self._search_results, key=lambda x: float(x.CC_F_Z_score), reverse=True)[:self._max_to_keep]

    @property
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    @staticmethod
    def cleanup(logfile):
        """Simple function to clean up log files after a run"""
        os.remove(logfile)

    @staticmethod
    def calculate_integration_box(model):
        """Function to calculate the integration radius or minimal box for an input PDB

        Parameters
        ----------
        model : str
            Path to input model

        Returns
        -------
        float
            The X coordinate
        float
            The Y coordinate
        float
            The Z coordinate
        float
            The integration radius for spherical structure

        """
        pdb_input = iotbx.pdb.pdb_input(file_name=model)
        hierarchy = pdb_input.construct_hierarchy()

        # Get resolution
        x = pdb_input.extract_remark_iii_records(2)
        resolution = iotbx.pdb.mining.extract_best_resolution(x)

        # Set a default resolution if mining fails
        if resolution == None:
            resolution = 2.0

        # Get a list of all xyz coordinates
        xyz = numpy.zeros((0, 3))
        for residue_group in hierarchy.models()[0].chains()[0].residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():
                    xyz = numpy.vstack([xyz, atom.xyz])
        
        # Get the smallest box containing the model
        #   numpy.ptp() ==> "Range of values (maximum - minimum) along an axis"
        diffs = numpy.asarray([
            numpy.ptp(xyz[:, 0]),
            numpy.ptp(xyz[:, 1]),
            numpy.ptp(xyz[:, 2])
        ])
        # Get integration radius (note, for spherical structure)
        intrad = diffs.min() * 0.75
        
        # Add together for each coordinate
        x, y, z = diffs + intrad + resolution
        
        return x.item(), y.item(), z.item(), intrad.item()

    def amore_run(self, models_dir, logs_dir, nproc=2, shres=3.0, pklim=0.5, npic=50, rotastep=1.0):
        """Run amore rotation function on a directory of models

        Parameters
        ----------
        models_dir : str
            The directory containing the models to run the rotation search on
        logs_dir : str
            The directory where logs from the job will be placed
        nproc : int
            The number of processors to run the job on
        shres : int float
            Spherical harmonic resolution [default 3.0]
        pklim : int float
            Peak limit, output all peaks above <float> [default: 0.5]
        npic : int
            Number of peaks to output from the translation function map for each orientation [default: 50]
        rotastep : int float
            Size of rotation step [default : 1.0]

        Returns
        -------
        file
            log file for each model in the models_dir

        """

        # make logs directory if it hasn't already been made
        if not os.path.isdir(logs_dir):
            os.mkdir(logs_dir)

        job_queue = multiprocessing.Queue()

        def run(job_queue, timeout=60):
            """processes element of job queue if queue not empty"""
            while not job_queue.empty():
                model = job_queue.get(timeout=timeout)
                self._amore_run(model, logs_dir, shres, pklim, npic, rotastep)

        for e in os.walk(models_dir):
            for model in e[2]:
                relpath = os.path.relpath(models_dir)
                job_queue.put(os.path.join(relpath, model))

        processes = []
        for i in range(nproc):
            process = multiprocessing.Process(target=run, args=(job_queue,))
            process.start()
            processes.append(process)

        for process in processes:
            process.join()

        if job_queue.empty():
            self.return_z_score_results(logs_dir)

    def _amore_run(self, model, logs_dir, shres, pklim, npic, rotastep):
        """Function to run tabfun and rotfun sequentially"""

        self.name = os.path.splitext(os.path.basename(model)[0:6])[0]

        # Make output directory if it doesn't exist
        output_dir = os.path.join(self.work_dir, 'output')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        logger.info("Running AMORE rotation function on {0}".format(self.name))

        # Set up variables for the run
        x, y, z, intrad = AmoreRotationSearch.calculate_integration_box(model)

        # Run tabfun
        self.tabfun(model, x, y, z)

        # Run rotfun
        self.rotfun(logs_dir, shres, intrad, pklim, npic, rotastep)

        return

    def calculate_intr_box(self, model):
        """Function to calculate the integration radius or minimal box for an input PDB

        Parameters
        ----------
        model : str
            Path to input model

        See Also
        --------
        calculate_integration_box

        """
        msg = "Deprecated - This function will be removed in a future update. Use calculate_integration_box instead"
        warnings.warn(msg)
        return AmoreRotationSearch.calculate_integration_box(model)

    def matthews_coef(self, model, min_solvent_content=30):
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

        # Get the molecular weight of the input model
        molecular_weight = self.rwcontents(model)

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

    def return_z_score_results(self, log_dir):
        """Function to return the z_score results

        Parameters
        ----------
        log_dir : str
            Path to the directory containing the logs from the amore run

        Returns
        -------
        class object
            Class object containing each model and its associated scores from AMORE
        """
        results = []
        for e in os.walk(log_dir):
            for log in e[2]:
                for line in open(os.path.join(log_dir, log)):
                    if line.startswith(" SOLUTIONRCD "):
                        fields = line.split()

                        if float(fields[-3]) > 0:
                            try:
                                ALPHA = float(fields[2])
                                BETA = float(fields[3])
                                GAMMA = float(fields[4])
                                CC_F = float(fields[8])
                                RF_F = float(fields[9])
                                CC_I = float(fields[10])
                                CC_P = float(fields[11])
                                Icp = float(fields[12])
                                CC_F_Z_score = float(fields[-3])
                                CC_P_Z_score = float(fields[-2])
                                Num_of_rot = float(fields[-1])

                            except ValueError:
                                ALPHA = float(fields[2])
                                BETA = float(fields[3])
                                GAMMA = float(fields[4])
                                CC_F = 'N/A'
                                RF_F = 'N/A'
                                CC_I = 'N/A'
                                CC_P = 'N/A'
                                Icp = 'N/A'
                                CC_F_Z_score = float(fields[-3])
                                CC_P_Z_score = float(fields[-2])
                                Num_of_rot = float(fields[-1])

                            break
                if 'clogs' in log_dir:
                    pdb_code = log[0:6]
                else:
                    pdb_code = log[0:7]

                score = _AmoreRotationScore(pdb_code, ALPHA, BETA, GAMMA, CC_F, RF_F, CC_I, CC_P, Icp, CC_F_Z_score,
                                            CC_P_Z_score, Num_of_rot)
                results.append(score)

        self._search_results = results

        return

    def rotfun(self, logs_dir, shres, intrad, pklim, npic, rotastep):
        """Function to perform first amore rotation function,

        Parameters
        ----------
        logs_dir : str
            The directory where logs from the job will be placed
        nproc : int
            The number of processors to run the job on
        shres : int float
            Spherical harmonic resolution [default 3.0]
        pklim : int float
            Peak limit, output all peaks above <float> [default: 0.5]
        npic : int float
            Number of peaks to output from the translation function map for each orientation [default: 50]
        rotastep : int float
            Size of rotation step [default : 1.0]

        Returns
        -------
        file
            A log file in the logs_dir containing information about the run
        """

        cmd = [self.amore_exe,
               'table1', os.path.join(self.work_dir, 'output', '{0}_sfs.tab'.format(self.name)),
               'HKLPCK1', os.path.join(self.work_dir, 'output', '{0}.hkl'.format(self.name)),
               'hklpck0', os.path.join(self.work_dir, 'spmipch.hkl'),
               'clmn1', os.path.join(self.work_dir, 'output', '{0}.clmn'.format(self.name)),
               'clmn0', os.path.join(self.work_dir, 'output', '{0}_spmipch.clmn'.format(self.name)),
               'MAPOUT', os.path.join(self.work_dir, 'output', 'amore_cross.map')]
        command_line = os.linesep.join(map(str, cmd))

        key = """ROTFUN
TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1
GENE 1   RESO 100.0 {0}  CELL_MODEL 80 75 65
CLMN CRYSTAL ORTH  1 RESO  20.0  {0}  SPHERE   {1}
CLMN MODEL 1     RESO  20.0  {0} SPHERE   {1}
ROTA  CROSS  MODEL 1  PKLIM {2}  NPIC {3} STEP {4}""".format(shres,
                                                             intrad,
                                                             pklim,
                                                             npic,
                                                             rotastep)

        logfile = os.path.join(self.work_dir, logs_dir, '{0}.log'.format(self.name))

        simbad_util.run_job(command_line, logfile, key)

        return

    def rwcontents(self, model):
        """Function to run rwcontents to get the molecular weight of a model

        Parameters
        ----------
        model : str
            Path to input model

        Returns
        -------
        float
            Molecular weight of input model
        """

        cmd = ['rwcontents',
               'xyzin', model]

        name = os.basename(model)[0:3]
        logfile = 'rwcontents_{0}.log'.format(name)
        simbad_util.run_job(cmd, logfile)

        # Exctract molecular weight from log file
        molecular_weight = None
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith(" Molecular Weight of protein"):
                    molecular_weight = float(line.split()[-1])
        if not molecular_weight:
            msg = "Cannot find Molecular weight in logfile {0]".format(logfile)
            logger.debug(msg)
            raise RuntimeError(msg)

        # Clean up
        os.remove(logfile)

        return molecular_weight

    def sortfun(self):
        """A function to prepare files for amore rotation function

        Parameters
        ----------
        self.mtz : str
            mtz file input to AmoreRotationSearch()
        self.work_dir : str
            working directory input to AmoreRotationSearch()

        Returns
        -------
        file
            spmipch.hkl file needed for rotfun and tabfun

        """

        logger.info("Preparing files for AMORE rotation function")

        # Get column labels for f and sigf
        f,sigf,_,_,_ = mtz_util.get_labels(self.mtz)

        cmd = [self.amore_exe,
               'hklin', self.mtz,
               'hklpck0', os.path.join(self.work_dir, 'spmipch.hkl')]

        command_line = os.linesep.join(map(str, cmd))

        key = """TITLE   ** spmi  packing h k l F for crystal**
SORTFUN RESOL 100.  2.5
LABI FP={0}  SIGFP={1}""".format(f, sigf)

        logfile = os.path.join(self.work_dir, 'SORTFUN.log')
        simbad_util.run_job(command_line, logfile, key)

        self.cleanup(logfile)

        return

    def summarize(self, csv_file='amore.csv'):
        """Summarize the search results

        Parameters
        ----------
        csv_file : str
           The path for a backup CSV file

        Raises
        ------
            No results found
        """

        search_results = self.search_results
        if not search_results:
            msg = "No results found"
            raise RuntimeError(msg)

        df = pandas.DataFrame(
            [r._as_dict() for r in search_results],
            index=[r.pdb_code for r in search_results],
            columns=["ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
                     "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak"],
        )
        # Create a CSV for reading later
        df.to_csv(csv_file)
        # Display table in stdout
        summary_table = """
The AMORE rotation search found the following structures:

%s
"""
        logger.info(summary_table % df.to_string())

    def tabfun(self, model, x=200, y=200, z=200):
        """Function to perform amore table function

        Parameters
        ----------
        model : str
            Path to input model
        x : int float
            x value for minimal box [default: 200]
        y : int float
            y value for minimal box [default: 200]
        z : int float
            z value for minimal box [default: 200]


        Returns
        -------
        file
            Output PDB for use in rotfun
        file
            Output table file for use in rotfun

        """

        cmd = [self.amore_exe,
               'xyzin1', model,
               'xyzout1', os.path.join(self.work_dir, 'output', '{0}_rot.pdb'.format(self.name)),
               'table1', os.path.join(self.work_dir, 'output', '{0}_sfs.tab'.format(self.name))]
        command_line = os.linesep.join(map(str, cmd))

        key = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL {0} {1} {2} 90 90 120 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0""".format(x, y, z)

        logfile = os.path.join(self.work_dir, '{0}_tabfun.log'.format(self.name))

        simbad_util.run_job(command_line, logfile, key)

        self.cleanup(logfile)

        return
