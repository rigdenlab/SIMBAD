"""Module to run the AMORE rotation search"""

from __future__ import print_function

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "07 Mar 2017"
__version__ = "0.1"

import copy_reg
import logging
import numpy
import os
import pandas
import types

from simbad.parsers import rotsearch_parser
from simbad.util import mtz_util
from simbad.util import simbad_util
from simbad.util import workers_util

import iotbx.pdb
import iotbx.pdb.mining

logger = logging.getLogger(__name__)


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


class _AmoreRotationScore(object):
    """An amore rotation scoring class"""

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
        return {k: getattr(self, k) for k in self.__slots__}


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
    min_solvent_content : int float
        The minimum solvent content present in the unit cell with the input model [default: 20]

    Examples
    --------
    >>> from simbad.rotsearch.amore_search import AmoreRotationSearch
    >>> rotation_search = AmoreRotationSearch('<amore_exe>', '<mtz>', '<work_dir>', '<max_to_keep>')
    >>> rotation_search.sortfun()
    >>> rotation_search.amore_run('<models_dir>', '<logs_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>',
    ...                           '<min_solvent_content>')
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
        if resolution is None:
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

    def amore_run(self, models_dir, logs_dir, nproc=2, shres=3.0, pklim=0.5, npic=50, 
                  rotastep=1.0, min_solvent_content=20, submit_cluster=False, submit_qtype=None, 
                  submit_queue=False, submit_array=None, submit_max_array=None, monitor=None):
        """Run amore rotation function on a directory of models

        Parameters
        ----------
        models_dir : str
            The directory containing the models to run the rotation search on
        logs_dir : str
            The directory where logs from the job will be placed
        nproc : int, optional
            The number of processors to run the job on
        shres : int, float, optional
            Spherical harmonic resolution [default 3.0]
        pklim : int, float, optional
            Peak limit, output all peaks above <float> [default: 0.5]
        npic : int, optional
            Number of peaks to output from the translation function map for each orientation [default: 50]
        rotastep : int, float, optional
            Size of rotation step [default : 1.0]
        min_solvent_content : int, float, optional
            The minimum solvent content present in the unit cell with the input model [default: 30]

        Returns
        -------
        file
            log file for each model in the models_dir

        """
        
        # make logs directory if it hasn't already been made
        if not os.path.isdir(logs_dir):
            os.mkdir(logs_dir)

        # Get the space group and cell parameters for the input mtz
        space_group, _, cell_parameters = mtz_util.crystal_data(self.mtz)
        
        log_files = []
        job_scripts = []
        for e in os.walk(models_dir):
            for model in e[2]:
                relpath = os.path.relpath(models_dir)
                input_model = os.path.join(relpath, model)
                self.name = os.path.basename(model).split('.')[0]
                
                # Ignore models below minimum solvent content
                if self.matthews_coef(input_model, cell_parameters, space_group, min_solvent_content):
                    output_dir = os.path.join(self.work_dir, 'output')
                    if not os.path.exists(output_dir):
                        os.mkdir(output_dir)
                    
                    logger.debug("Generating script to perform AMORE rotation function on %s", self.name)
                        
                    # Set up variables for the run
                    x, y, z, intrad = AmoreRotationSearch.calculate_integration_box(input_model)
                    tab_cmd, tab_key = self.tabfun(input_model, x, y, z)
                    rot_cmd, rot_key = self.rotfun(logs_dir, shres, intrad, pklim, npic, rotastep)
                    
                    script = simbad_util.tmp_file_name(delete=False, suffix=simbad_util.SCRIPT_EXT)
                    logfile = os.path.join(self.work_dir, logs_dir, '{0}.log'.format(self.name))
                    with open(script, 'w') as f_out:
                        f_out.write(simbad_util.SCRIPT_HEADER + os.linesep * 2)
                        f_out.write(" ".join(map(str, tab_cmd)) + " << eof" + os.linesep)
                        f_out.write(tab_key + os.linesep + "eof" + os.linesep * 2)
                        f_out.write(" ".join(map(str, rot_cmd)) + " << eof > " + logfile + os.linesep)
                        f_out.write(rot_key + os.linesep + "eof" + os.linesep * 2)
                    
                    os.chmod(script, 0o777)
                    job_scripts.append(script)
                    log_files.append(logfile)
                else:
                    msg = "Skipping {0}: solvent content is predicted to be less than {1}".format(self.name, min_solvent_content)
                    logger.debug(msg)
                                        
        # Execute the scripts
        success = workers_util.run_scripts(
            job_scripts=job_scripts,
            monitor=monitor,
            check_success=None,
            early_terminate=False,
            nproc=nproc,
            job_time=7200,
            job_name='simbad_rot',
            submit_cluster=submit_cluster,
            submit_qtype=submit_qtype,
            submit_queue=submit_queue,
            submit_array=submit_array,
            submit_max_array=submit_max_array,
        )

        if not success == None:
            results = []
            for logfile in log_files:
                RP = rotsearch_parser.RotsearchParser(logfile)
               
                pdb_code = os.path.basename(logfile).split('.')[0]
                
                score = _AmoreRotationScore(pdb_code, RP.alpha, RP.beta, RP.gamma, RP.cc_f, RP.rf_f, RP.cc_i, 
                                            RP.cc_p, RP.icp, RP.cc_f_z_score, RP.cc_p_z_score, RP.num_of_rot)
                
                # Ignore results for searches which didn't work
                if not RP.cc_f_z_score == None:
                    results.append(score)

            self._search_results = results
        return

    def matthews_coef(self, model, cell_parameters, space_group, min_solvent_content=30):
        """Function to run matthews coefficient to decide if the model can fit in the unit cell

        Parameters
        ----------
        model : str
            Path to input model
        cell_parameters : str
            The parameters describing the unit cell of the crystal
        space_group : str
            The space group of the crystal
        min_solvent_content : int, float, optional
            Minimum solvent content [default: 30]

        Returns
        -------
        bool
            Can the model fit in the unit cell with a solvent content higher than the min_solvent_content

        """
        # Get the molecular weight of the input model
        molecular_weight = simbad_util.molecular_weight(model)

        cmd = ["matthews_coef"]
        key = """CELL {0}
        symm {1}
        molweight {2}
        auto""".format(cell_parameters,
                       space_group,
                       molecular_weight)

        logfile = os.path.join(self.work_dir, 'matt_coef_{0}.log'.format(self.name))
        simbad_util.run_job(cmd, logfile=logfile, stdin=key)

        # Determine if the model can fit in the unit cell
        solvent_content = 0
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

    def rotfun(self, logs_dir, shres, intrad, pklim, npic, rotastep):
        """Function to perform first amore rotation function,

        Parameters
        ----------
        logs_dir : str
            The directory where logs from the job will be placed
        nproc : int
            The number of processors to run the job on
        intrad : int, float
            The tolerance
        shres : int, float
            Spherical harmonic resolution
        pklim : int, float
            Peak limit, output all peaks above :obj:`float`
        npic : int, float
            Number of peaks to output from the translation function map for each orientation
        rotastep : int, float
            Size of rotation step

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

        return cmd, key

    def sortfun(self):
        """A function to prepare files for amore rotation function

        Parameters
        ----------
        self.mtz : str
            mtz file input to :obj:`AmoreRotationSearch`
        self.work_dir : str
            working directory input to :obj:`AmoreRotationSearch`

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

        key = """TITLE   ** spmi  packing h k l F for crystal**
SORTFUN RESOL 100.  2.5
LABI FP={0}  SIGFP={1}""".format(f, sigf)

        logfile = os.path.join(self.work_dir, 'SORTFUN.log')
        simbad_util.run_job(cmd, logfile=logfile, stdin=key)
        self.cleanup(logfile)

    def summarize(self, csv_file):
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
        df.to_csv(os.path.join(self.work_dir, csv_file))
        # Display table in stdout
        summary_table = """
The AMORE rotation search found the following structures:

%s
"""
        logger.info(summary_table, df.to_string())

    def tabfun(self, model, x=200, y=200, z=200):
        """Function to perform amore table function

        Parameters
        ----------
        model : str
            Path to input model
        x : int, float, optional
            x value for minimal box [default: 200]
        y : int, float, optional
            y value for minimal box [default: 200]
        z : int, float, optional
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
               'xyzout1', os.path.join(self.work_dir, 'output', '{0}.pdb'.format(self.name)),
               'table1', os.path.join(self.work_dir, 'output', '{0}_sfs.tab'.format(self.name))]

        key = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL {0} {1} {2} 90 90 120 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0""".format(x, y, z)

        return cmd, key

