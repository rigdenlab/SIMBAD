"""Class to run MR on SIMBAD results using code from MrBump"""

from __future__ import print_function

__author__ = "Adam Simpkin"
__date__ = "09 Mar 2017"
__version__ = "1.0"

import copy_reg
import logging
import multiprocessing
import os
import pandas
import simbad.mr
import tempfile
import types

from simbad.mr import anomalous_util
from simbad.parsers import molrep_parser
from simbad.parsers import phaser_parser
from simbad.parsers import refmac_parser
from simbad.util import mtz_util
from simbad.util import simbad_util
from simbad.util import workers_util

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
        """Convert the :obj:`_MrScore <simbad.mr._MrScore>`
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
    enant : bool
        Test enantimorphic space groups [default: False]
    results : class
        Results from :obj: '_LatticeParameterScore' or :obj: '_AmoreRotationScore'
    time_out : int, optional
        Number of seconds for multiprocessing job to timeout [default: 60]
    nproc : int, optional
        Number of processors to use [default: 2]
    submit_cluster : bool
        specify whether to submit to a cluster [default: False]
    
    TODO: Put these here for now but need to make a submit obj    
    
    submit_qtype : str
    submit_queue : bool
    submit_array : str
    submit_max_array : str
    monitor : str

    Examples
    --------
    >>> from simbad.mr import MrSubmit
    >>> MR = MrSubmit('<mtz>', '<mr_program>', '<refine_program>', '<model_dir>', '<output_dir>', '<early_term>',
    >>>               '<enam>')
    >>> MR.submit_jobs('<results>', '<time_out>', '<nproc>', '<submit_cluster>', '<submit_qtype>', 
    >>>                '<submit_queue>', '<submit_array>', '<submit_max_array>', '<early_terminate>',
    >>>                '<monitor>')

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
        # Catch arguments added in a strings
        if early_term == 'False':
            early_term = False
        elif early_term == 'True':
            early_term = True
            
        self._early_term = early_term

    @property
    def enant(self):
        """Flag to decide if enantiomorphic spacegroups should be trialled"""
        return self._enant

    @enant.setter
    def enant(self, enant):
        """Set the enant flag to true or false"""
        # Catch arguments added in a strings
        if enant == 'False':
            enant = False
        elif enant == 'True':
            enant = True
        
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

    def submit_jobs(self, results, time_out=7200, nproc=1, submit_cluster=False, submit_qtype=None, 
                    submit_queue=False, submit_array=None, submit_max_array=None, monitor=None):

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

        # Set up path to MR executable files
        exe_path = os.path.dirname(simbad.mr.__file__)

        if self.mr_program.upper() == 'MOLREP':
            mr_exectutable = os.path.join(exe_path, 'molrep_mr.py')
        elif self.mr_program.upper() == 'PHASER':
            mr_exectutable = os.path.join(exe_path, 'phaser_mr.py')

        if self.refine_program.upper() == "REFMAC5":
            ref_exectutable = os.path.join(exe_path, 'refmac_refine.py')

        job_scripts = []
        log_files = []
        for result in results:
            mr_pdbin = os.path.join(self.model_dir, '{0}.pdb'.format(result.pdb_code))
            mr_workdir = os.path.join(self.output_dir, result.pdb_code, 'mr', self.mr_program)
            mr_logfile = os.path.join(mr_workdir, '{0}_mr.log'.format(result.pdb_code))
            mr_pdbout = os.path.join(mr_workdir, '{0}_mr_output.pdb'.format(result.pdb_code))
     
            ref_workdir = os.path.join(mr_workdir, 'refine')
            ref_hklout = os.path.join(ref_workdir, '{0}_refinement_output.mtz'.format(result.pdb_code))
            ref_logfile = os.path.join(ref_workdir, '{0}_ref.log'.format(result.pdb_code))
            ref_pdbout = os.path.join(ref_workdir, '{0}_refinement_output.pdb'.format(result.pdb_code))

            # TODO: don't just skip, remove the entry from the actual search_results
            if not os.path.isfile(mr_pdbin):
                continue
             
            # Common MR keywords
            mr_cmd = [mr_exectutable, "-enant", self.enant, "-hklin", self.mtz, "-pdbin", mr_pdbin,
                      "-pdbout", mr_pdbout, "-logfile", mr_logfile, "-work_dir", mr_workdir]

            # Common refine keywords
            ref_cmd = [ref_exectutable, "-hklout", ref_hklout, "-pdbin", mr_pdbout,
                       "-pdbout", ref_pdbout, "-logfile", ref_logfile, "-work_dir", ref_workdir]

            # Extend commands with program-specific options
            if self.mr_program.upper() == 'MOLREP':
                mr_cmd += ["-space_group", self.space_group]
                ref_cmd += ["-hklin", self.mtz]
     
            elif self.mr_program.upper() == 'PHASER':
                hklout = os.path.join(mr_workdir, '{0}_mr_output.mtz'.format(model.pdb_code))
                mr_cmd += [
                    "-f", self.f,
                    "-hklout", hklout,
                    "-sigf", self.sigf,
                    "-solvent", self.solvent,
                ]
                ref_cmd += ["-hklin", hklout]
            
            # Create a run script
            prefix = result.pdb_code + '_simbad'
            log = os.path.join(self.output_dir, prefix + '.log')
            script = os.path.join(self.output_dir, prefix + simbad_util.SCRIPT_EXT)
            with open(script, 'w') as f_out:
                f_out.write(simbad_util.SCRIPT_HEADER + os.linesep * 2)
                f_out.write(" ".join(map(str, mr_cmd)) + os.linesep * 2)
                f_out.write(" ".join(map(str, ref_cmd)) + os.linesep)
            os.chmod(script, 0o777)
            job_scripts.append(script)
            log_files.append(log)
                
        # Execute the scripts
        success = workers_util.run_scripts(
            job_scripts=job_scripts,
            monitor=monitor,
            check_success=mr_job_succeeded,
            early_terminate=self.early_term,
            nproc=nproc,
            job_time=time_out,
            job_name='simbad_mr',
            submit_cluster=submit_cluster,
            submit_qtype=submit_qtype,
            submit_queue=submit_queue,
            submit_array=submit_array,
            submit_max_array=submit_max_array,
        )
        
        if success:
            for result in results:
                # Set default values
                molrep_score = None
                molrep_tfscore = None
                phaser_tfz = None
                phaser_llg = None
                phaser_rfz = None
                final_r_fact = None
                final_r_free = None
                peaks_over_6_rms = None
                peaks_over_6_rms_within_2A_of_model = None
                peaks_over_12_rms = None
                peaks_over_12_rms_within_2A_of_model = None
                
                if self.mr_program.lower() == "molrep":
                    MP = molrep_parser.MolrepParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 
                                                                 'molrep', '{0}_mr.log'.format(result.pdb_code)))
                    molrep_score = MP.score
                    molrep_tfscore = MP.tfscore
                    
                elif self.mr_program.lower() == "phaser":
                    PP = phaser_parser.PhaserParser(os.path.join(self.output_dir, result.pdb_code, 'mr', 
                                                                 'phaser', '{0}_mr.log'.format(result.pdb_code)))
                    phaser_tfz = PP.tfz
                    phaser_llg = PP.llg
                    phaser_rfz = PP.rfz
                
                if self._dano is not None:
                    AS = anomalous_util.AnomSearch(self.mtz, self.output_dir, self.mr_program)
                    AS.run(result)
                    a = AS.search_results()
                    
                    peaks_over_6_rms = a.peaks_over_6_rms
                    peaks_over_6_rms_within_2A_of_model = a.peaks_over_6_rms_within_2A_of_model
                    peaks_over_12_rms = a.peaks_over_12_rms
                    peaks_over_12_rms_within_2A_of_model = a.peaks_over_12_rms_within_2A_of_model
                
                
                RP = refmac_parser.RefmacParser(os.path.join(self.output_dir, result.pdb_code, 'mr', self.mr_program, 
                                                             'refine', '{0}_ref.log'.format(result.pdb_code)))
                final_r_free = RP.final_r_free
                final_r_fact = RP.final_r_fact
            
                score = _MrScore(pdb_code=result.pdb_code, 
                                 molrep_score=molrep_score,
                                 molrep_tfscore=molrep_tfscore, 
                                 phaser_tfz=phaser_tfz, 
                                 phaser_llg=phaser_llg, 
                                 phaser_rfz=phaser_rfz, 
                                 final_r_fact=final_r_fact,
                                 final_r_free=final_r_free, 
                                 peaks_over_6_rms=peaks_over_6_rms,
                                 peaks_over_6_rms_within_2A_of_model=peaks_over_6_rms_within_2A_of_model,
                                 peaks_over_12_rms=peaks_over_12_rms,
                                 peaks_over_12_rms_within_2A_of_model=peaks_over_12_rms_within_2A_of_model
                                 )
                self._search_results.append(score)
                    
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

        logfile = 'matt_coef.log'
        ret = simbad_util.run_job(cmd, logfile=logfile, stdin=key)
        if ret != 0:
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
        logger.info(summary_table, df.to_string())


def mr_job_succeeded(script):
    """Check a Molecular Replacement job for it's success
    
    Parameters
    ----------
    script : str
       The path to the execution script

    Returns
    -------
    bool
       Success status of the MR run

    """
    for line in open(script, 'r'):
        if 'refmac_refine.py' in line:
            line = line.strip().split()
            logfile = line[line.index('-logfile') + 1]
    if os.path.isfile(logfile):
        RP = refmac_parser.RefmacParser(logfile)
        if RP.final_r_fact < 0.45 and RP.final_r_free < 0.45:
            return True
        return False
    else:
        logger.critical("Cannot find logfile: %s", logfile)
        return False

