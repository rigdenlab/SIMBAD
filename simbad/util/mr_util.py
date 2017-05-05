"""Class to run MR on SIMBAD results using code from MrBump"""

import copy_reg
import logging
import multiprocessing
import os
import pandas
import types

from simbad.parsers import molrep_parser
from simbad.parsers import phaser_parser
from simbad.parsers import refmac_parser
from simbad.util import anomalous_util
from simbad.util import molrep_util
from simbad.util import mtz_util
from simbad.util import phaser_util
from simbad.util import refmac_util
from simbad.util import simbad_util

__author__ = "Adam Simpkin"
__date__ = "09 Mar 2017"
__version__ = "0.2"

logger = logging.getLogger(__name__)


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


copy_reg.pickle(types.MethodType, _pickle_method)


class _MrScore(object):
    """A molecular replacement scoring class"""

    __slots__ = ("pdb_code", "final_r_fact", "final_r_free", "molrep_score", "molrep_tfscore",
                 "phaser_tfz", "phaser_llg", "phaser_rfz", "peaks_over_6_rms", "peaks_over_6_rms_within_2A_of_model",
                 "peaks_over_12_rms", "peaks_over_12_rms_within_2A_of_model")

    def __init__(self, pdb_code, final_r_fact, final_r_free, molrep_score=None, molrep_tfscore=None, phaser_tfz=None,
                 phaser_llg=None, phaser_rfz=None, peaks_over_6_rms=None, peaks_over_6_rms_within_2A_of_model=None,
                 peaks_over_12_rms=None, peaks_over_12_rms_within_2A_of_model=None):
        self.pdb_code = pdb_code
        self.molrep_score = molrep_score
        self.molrep_tfscore = molrep_tfscore
        self.phaser_tfz = phaser_tfz
        self.phaser_llg = phaser_llg
        self.phaser_rfz = phaser_rfz
        self.final_r_fact = final_r_fact
        self.final_r_free = final_r_free
        self.peaks_over_6_rms = peaks_over_6_rms
        self.peaks_over_6_rms_within_2A_of_model = peaks_over_6_rms_within_2A_of_model
        self.peaks_over_12_rms = peaks_over_12_rms
        self.peaks_over_12_rms_within_2A_of_model = peaks_over_12_rms_within_2A_of_model

    def __repr__(self):
        return "{0}(pdb_code={1}  final_r_fact={2} final_r_free={3} molrep_score={4} molrep_tfscore={5} " \
               "phaser_tfz={6}, phaser_llg={7}, phaser_rfz={8}, peaks_over_6_rms={9}, " \
               "peaks_over_6_rms_within_2A_of_model={10}, peaks_over_12_rms={11}, " \
               "peaks_over_12_rms_within_2A_of_model={12})".format(self.__class__.__name__, self.pdb_code,
                                                                   self.final_r_fact, self.final_r_free,
                                                                   self.molrep_score, self.molrep_tfscore,
                                                                   self.phaser_tfz, self.phaser_llg,
                                                                   self.phaser_rfz, self.peaks_over_6_rms,
                                                                   self.peaks_over_6_rms_within_2A_of_model,
                                                                   self.peaks_over_12_rms,
                                                                   self.peaks_over_12_rms_within_2A_of_model)

    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.util.mr_util._MrScore>`
        object to a dictionary"""
        dictionary = {}
        for k in self.__slots__:
            dictionary[k] = getattr(self, k)
        return dictionary


class MrSubmit(object):
    """Class to run MR on a defined set of models

    Attributes
    ----------
    mtz : str
        Path to the input MTZ file
    mr_program : str
        Name of the molecular replacement program to use
    refine_program : str
        Name of the refinement program to use
    model_dir : str
        Path to the directory containing input models
    output_dir : str
        Path to the directory to output results
    early_term : bool
        Terminate early if a solution is found [default: True]
    enam : bool
        Test enantimorphic space groups [default: False]
    results : class
        Results from :obj: '_LatticeParameterScore' or :obj: '_AmoreRotationScore'
    time_out : int, optional
        Number of seconds for multiprocessing job to timeout [default: 60]
    nproc : int, optional
        Number of processors to use [default: 2]

    Examples
    --------
    >>> from simbad.util.mr_util import MrSubmit
    >>> MR = MrSubmit('<mtz>', '<mr_program>', '<refine_program>', '<model_dir>', '<output_dir>', '<early_term>',
    >>>               '<enam>')
    >>> MR.multiprocessing('<results>', '<time_out>', '<nproc>')

    If a solution is found and early_term is set to True, the queued jobs will be terminated.
    """

    def __init__(self, mtz, mr_program, refine_program, model_dir, output_dir, early_term=True, enant=False):
        """Initialise MrSubmit class"""
        self.input_file = None
        self._early_term = None
        self._enant = None
        self._mtz = None
        self._mr_program = None
        self._output_dir = None
        self._refine_program = None
        self._search_results = []

        # options derived from the input mtz
        self._cell_parameters = None
        self._resolution = None
        self._solvent = None
        self._space_group = None
        self._f = None
        self._sigf = None
        self._dano = None
        self._sigdano = None
        self._free = None

        self.early_term = early_term
        self.enant = enant
        self.model_dir = os.path.abspath(model_dir)
        self.mtz = mtz
        self.mr_program = mr_program
        self.output_dir = output_dir
        self.refine_program = refine_program

    @property
    def early_term(self):
        """Flag to decide if the program should terminate early"""
        return self._early_term

    @early_term.setter
    def early_term(self, early_term):
        """Set the early term flag to true or false"""
        self._early_term = early_term

    @property
    def enant(self):
        """Flag to decide if enantiomorphic spacegroups should be trialled"""
        return self._enant

    @enant.setter
    def enant(self, enant):
        """Set the enant flag to true or false"""
        self._enant = enant

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = os.path.abspath(mtz)
        self.get_mtz_info(mtz)

    @property
    def cell_parameters(self):
        """The cell parameters of the input MTZ file"""
        return self._cell_parameters

    @property
    def resolution(self):
        """The resolution of the input MTZ file"""
        return self._resolution

    @property
    def search_results(self):
        """The results from the amore rotation search"""
        return sorted(self._search_results, key=lambda x: float(x.final_r_free), reverse=False)

    @property
    def solvent(self):
        """The predicted solvent content of the input MTZ file"""
        return self._solvent

    @property
    def space_group(self):
        """The space group of the input MTZ file"""
        return self._space_group

    @property
    def f(self):
        """The F column label of the input MTZ file"""
        return self._f

    @property
    def sigf(self):
        """The SIGF column label of the input MTZ file"""
        return self._sigf

    @property
    def free(self):
        """The FREE column label of the input MTZ file"""
        return self._free

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
    def output_dir(self):
        """The path to the output directory"""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        """Define the output directory"""
        self._output_dir = output_dir

    def _run_job(self, model):
        """Function to run MR on each model"""
        logger.info("Running MR and refinement on %s", model.pdb_code)

        try:
            os.mkdir(self.output_dir)
        except OSError:
            pass

        # Make output directories
        os.mkdir(os.path.join(self.output_dir, model.pdb_code))
        os.mkdir(os.path.join(self.output_dir, model.pdb_code, 'mr'))
        os.mkdir(os.path.join(self.output_dir, model.pdb_code, 'mr', self.mr_program))
        os.mkdir(os.path.join(self.output_dir, model.pdb_code, 'mr', self.mr_program, 'refine'))

        # Set up MR input paths
        mr_pdbin = os.path.join(self.model_dir, '{0}.pdb'.format(model.pdb_code))
        mr_workdir = os.path.join(self.output_dir, model.pdb_code, 'mr', self.mr_program)
        mr_logfile = os.path.join(mr_workdir, '{0}_mr.log'.format(model.pdb_code))
        mr_pdbout = os.path.join(mr_workdir, '{0}_mr_output.pdb'.format(model.pdb_code))

        # Set up refinement input paths
        ref_workdir = os.path.join(mr_workdir, 'refine')
        ref_hklout = os.path.join(ref_workdir, '{0}_refinement_output.mtz'.format(model.pdb_code))
        ref_logfile = os.path.join(ref_workdir, '{0}_ref.log'.format(model.pdb_code))
        ref_pdbout = os.path.join(ref_workdir, '{0}_refinement_output.pdb'.format(model.pdb_code))

        # Run job
        if self.mr_program.upper() == 'MOLREP':
            # Set up class with MOLREP input arguments
            molrep = molrep_util.Molrep(self.enant, self.mtz, mr_logfile, mr_pdbin, mr_pdbout, self.space_group,
                                        mr_workdir)
            # Run MOLREP
            molrep.run()

            # Set up class with REFMAC input arguments
            refmac = refmac_util.Refmac(self.mtz, ref_hklout, ref_logfile, mr_pdbout, ref_pdbout, ref_workdir)
            # Run REFMAC
            refmac.run()

        elif self.mr_program.upper() == 'PHASER':
            hklout = os.path.join(mr_workdir, '{0}_mr_output.mtz'.format(model.pdb_code))

            # Set up class with PHASER input arguments
            phaser = phaser_util.Phaser(self.enant, self.f, self.mtz, hklout, mr_logfile, mr_pdbin, mr_pdbout,
                                        self.sigf, self.solvent, mr_workdir)
            # Run PHASER
            phaser.run()

            # Set up class with REFMAC input arguments
            refmac = refmac_util.Refmac(hklout, ref_hklout, ref_logfile, mr_pdbout, ref_pdbout, ref_workdir)
            # Run REFMAC
            refmac.run()

        if self.early_term and self.early_term != "False":
            try:
                terminate = self.solution_found(model)
                return terminate
            except:
                pass

        return False

    def get_mtz_info(self, mtz):
        """Get various information from the input MTZ

         Parameters
         ----------
         mtz : str
            Path to the input MTZ

        Returns
        -------
        self._cell_parameters : list
            The parameters that descibe the unit cell
        self._resolution : float
            The resolution of the data
        self._space_group : str
            The space group of the data
        self._f : str
            The column label for F
        self._sigf : str
            The column label for SIGF
        self._free : str
            The column label for FREE
        self._solvent : float
            The predicted solvent content of the protein
        """
        # Extract crystal data from input mtz
        self._space_group, _, self._cell_parameters = mtz_util.crystal_data(mtz)

        # Extract column labels from input mtz
        self._f, self._sigf, self._dano, self._sigdano, self._free = mtz_util.get_labels(mtz)

        # Get solvent content
        self._solvent = self.matthews_coef(self._cell_parameters, self._space_group)

    def multiprocessing(self, results, time_out=60, nproc=2):
        """Code to run MR and refinement in parallel

        Parameters
        ----------
        results : class
            Results from :obj: '_LatticeParameterScore' or :obj: '_AmoreRotationScore'
        time_out: int
            Number of seconds for multiprocessing job to timeout [default: 60]
        nproc : int
            Number of processors to use [default: 2]

        Returns
        -------
        file
            PDB from the molecular replacement job
        file
            Log file from the molecular replacement job
        file
            PDB from the refinement job
        file
            MTZ from the refinement job
        file
            Log file from the refinement job
        """

        def run(job_queue):
            """processes element of job queue if queue not empty"""
            while not job_queue.empty():
                model = job_queue.get(timeout=time_out)
                terminate = self._run_job(model)

                if terminate:
                    print "MR with {0} was successful so removing remaining jobs from inqueue".format(model.pdb_code)
                    while not job_queue.empty():
                        job = job_queue.get()
                        logger.debug("Removed job [%s] from inqueue", job.pdb_code)

        # Create job queue
        job_queue = multiprocessing.Queue()

        # Add each result from results to the job queue
        for result in results:
            job_queue.put(result)

        processes = []
        # Set up processes equal to the number of processors input
        for i in range(nproc):
            process = multiprocessing.Process(target=run, args=(job_queue,))
            process.start()
            processes.append(process)

        # block the calling thread
        for p in processes:
            p.join()

        if job_queue.empty():
            if self.mr_program.lower() == "molrep":
                for result in results:
                    try:
                        MP = molrep_parser.MolrepParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 'molrep',
                                                                     '{0}_mr.log'.format(result.pdb_code)))
                        molrep_score = MP.score
                        molrep_tfscore = MP.tfscore

                        RP = refmac_parser.RefmacParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 'molrep',
                                                                     'refine', '{0}_ref.log'.format(result.pdb_code)))
                        final_r_free = RP.final_r_free
                        final_r_fact = RP.final_r_fact

                        if self._dano is not None:
                            AS = anomalous_util.AnomSearch(self.mtz, self.output_dir, self.mr_program)
                            AS.run(result)
                            a = AS.search_results()

                            score = _MrScore(pdb_code=result.pdb_code, molrep_score=molrep_score,
                                             molrep_tfscore=molrep_tfscore, final_r_fact=final_r_fact,
                                             final_r_free=final_r_free, peaks_over_6_rms=a.peaks_over_6_rms,
                                             peaks_over_6_rms_within_2A_of_model=a.peaks_over_6_rms_within_2A_of_model,
                                             peaks_over_12_rms=a.peaks_over_12_rms,
                                             peaks_over_12_rms_within_2A_of_model=a.peaks_over_12_rms_within_2A_of_model
                                             )
                        else:
                            score = _MrScore(pdb_code=result.pdb_code, molrep_score=molrep_score,
                                             molrep_tfscore=molrep_tfscore, final_r_fact=final_r_fact,
                                             final_r_free=final_r_free)

                        self._search_results.append(score)
                    except:
                        pass
            elif self.mr_program.lower() == "phaser":
                for result in results:
                    try:
                        PP = phaser_parser.PhaserParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 'phaser',
                                                                      '{0}_mr.log'.format(result.pdb_code)))
                        phaser_tfz = PP.tfz
                        phaser_llg = PP.llg
                        phaser_rfz = PP.rfz

                        RP = refmac_parser.RefmacParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 'phaser',
                                                                     'refine', '{0}_ref.log'.format(result.pdb_code)))
                        final_r_free = RP.final_r_free
                        final_r_fact = RP.final_r_fact

                        if self._dano is not None:
                            AS = anomalous_util.AnomSearch(self.mtz, self.output_dir, self.mr_program)
                            AS.run(result)
                            a = AS.search_results()

                            score = _MrScore(pdb_code=result.pdb_code, phaser_tfz=phaser_tfz, phaser_llg=phaser_llg,
                                             phaser_rfz=phaser_rfz, final_r_fact=final_r_fact,
                                             final_r_free=final_r_free, peaks_over_6_rms=a.peaks_over_6_rms,
                                             peaks_over_6_rms_within_2A_of_model=a.peaks_over_6_rms_within_2A_of_model,
                                             peaks_over_12_rms=a.peaks_over_12_rms,
                                             peaks_over_12_rms_within_2A_of_model=a.peaks_over_12_rms_within_2A_of_model
                                             )
                        else:
                            score = _MrScore(pdb_code=result.pdb_code, phaser_tfz=phaser_tfz, phaser_llg=phaser_llg,
                                             phaser_rfz=phaser_rfz, final_r_fact=final_r_fact,
                                             final_r_free=final_r_free)

                        self._search_results.append(score)
                    except NameError:
                        pass

    def solution_found(self, model):
        """Function to check is a solution has been found

        Parameters
        ----------
        result_dir : str
            Path to the results directory
        model : class
            Class object containing the PDB code for the input model

        Returns
        -------
        bool
            Solution found <True|False>
        """
        RP = refmac_parser.RefmacParser(os.path.join(self.output_dir, model.pdb_code, 'mr', self.mr_program, 'refine',
                                                     '{0}_ref.log'.format(model.pdb_code)))

        final_r_free = RP.final_r_free
        final_r_fact = RP.final_r_fact

        if final_r_fact < 0.45 and final_r_free < 0.45:
            return True
        else:
            return False

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

        logfile = 'matt_coef.log'
        ret = simbad_util.run_job(cmd, logfile=logfile, stdin=key)
        if ret and ret != 0:
            msg = "matthews_coef exited with non-zero return code ({0}). Log is {1}".format(ret, logfile)
            logger.critical(msg)
            raise RuntimeError(msg)

        # Determine if the model can fit in the unit cell
        solvent_content = 0.5
        for line in open(logfile, 'r'):
            if line.startswith('  1'):
                solvent_content = (float(line.split()[2]) / 100)
                break

        os.remove(logfile)
        return solvent_content

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
        if search_results is None:
            msg = "No results found"
            raise RuntimeError(msg)

        # Decide on which column labels to output
        if self.mr_program.lower() == "molrep":
            if self._dano:
                columns = ["molrep_score", "molrep_tfscore", "final_r_fact", "final_r_free", "peaks_over_6_rms",
                           "peaks_over_6_rms_within_2A_of_model", "peaks_over_12_rms",
                           "peaks_over_12_rms_within_2A_of_model"]
            else:
                columns = ["molrep_score", "molrep_tfscore", "final_r_fact", "final_r_free"]
        elif self.mr_program.lower() == "phaser":
            if self._dano:
                columns = ["phaser_tfz", "phaser_llg", "phaser_rfz", "final_r_fact", "final_r_free", "peaks_over_6_rms",
                           "peaks_over_6_rms_within_2A_of_model", "peaks_over_12_rms",
                           "peaks_over_12_rms_within_2A_of_model"]
            else:
                columns = ["phaser_tfz", "phaser_llg", "phaser_rfz", "final_r_fact", "final_r_free"]

        df = pandas.DataFrame(
            [r._as_dict() for r in search_results],
            index=[r.pdb_code for r in search_results],
            columns=columns,
        )
        # Create a CSV for reading later
        df.to_csv(os.path.join(self.output_dir, csv_file))
        # Display table in stdout
        summary_table = """
MR/refinement gave the following results:

%s
"""
        logger.info(summary_table % df.to_string())
