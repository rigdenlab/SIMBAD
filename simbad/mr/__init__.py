"""Class to run MR on SIMBAD results"""

from __future__ import division

__author__ = "Adam Simpkin"
__date__ = "09 Mar 2017"
__version__ = "1.0"

import logging
import numpy
import os

from pyjob import Job, cexec
from pyjob.misc import make_script, tmp_file

from simbad.mr import anomalous_util
from simbad.parsers import molrep_parser
from simbad.parsers import phaser_parser
from simbad.parsers import refmac_parser
from simbad.util import pdb_edit
from simbad.util import matthews_coef
from simbad.util import mtz_util

logger = logging.getLogger(__name__)

# Make clear which binaries we currently support
KNOWN_MR_PROGRAMS = ["molrep", "phaser"]
KNOWN_REF_PROGRAMS = ["refmac5"]


class _MrScore(object):
    """A molecular replacement scoring class"""

    __slots__ = ("pdb_code", "final_r_fact", "final_r_free", "molrep_score", "molrep_tfscore",
                 "phaser_tfz", "phaser_llg", "phaser_rfz", "peaks_over_6_rms", "peaks_over_6_rms_within_4a_of_model",
                 "peaks_over_9_rms", "peaks_over_9_rms_within_4a_of_model")

    def __init__(self, pdb_code):
        self.pdb_code = pdb_code
        self.molrep_score = None
        self.molrep_tfscore = None
        self.phaser_tfz = None
        self.phaser_llg = None
        self.phaser_rfz = None
        self.final_r_fact = 1.0
        self.final_r_free = 1.0
        self.peaks_over_6_rms = None
        self.peaks_over_6_rms_within_4a_of_model = None
        self.peaks_over_9_rms = None
        self.peaks_over_9_rms_within_4a_of_model = None

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code}  final_r_fact={final_r_fact} final_r_free={final_r_free}"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.mr._MrScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}


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
    output_dir : str
        Path to the directory to output results
    enant : bool
        Test enantimorphic space groups [default: False]
    results : obj
        Results from :obj: '_LatticeParameterScore' or :obj: '_AmoreRotationScore'

    Examples
    --------
    >>> from simbad.mr import MrSubmit
    >>> MR = MrSubmit('<mtz>', '<mr_program>', '<refine_program>', '<output_dir>', '<enam>')
    >>> MR.submit_jobs('<results>', '<nproc>', '<submit_cluster>', '<submit_qtype>', '<submit_queue>',
    ...                '<submit_array>', '<submit_max_array>', '<timeout>', '<process_all>', '<monitor>')

    If a solution is found and process_all is not set, the queued jobs will be terminated.
    """

    def __init__(self, mtz, mr_program, refine_program, output_dir, timeout, enant=False):
        """Initialise MrSubmit class"""
        self.input_file = None
        self._process_all = None
        self._enant = None
        self._mtz = None
        self._mr_program = None
        self._output_dir = None
        self._refine_program = None
        self._search_results = []
        self._timeout = None

        # options derived from the input mtz
        self._cell_parameters = None
        self._resolution = None
        self._space_group = None
        self._f = None
        self._sigf = None
        self._i = None
        self._sigi = None
        self._dano = None
        self._sigdano = None
        self._free = None

        self.enant = enant
        self.mtz = mtz
        self.mr_program = mr_program
        self.output_dir = output_dir
        self.refine_program = refine_program
        self.timeout = timeout

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
    def i(self):
        """The I column label of the input MTZ file"""
        return self._i

    @property
    def sigi(self):
        """The SIGI column label of the input MTZ file"""
        return self._sigi

    @property
    def free(self):
        """The FREE column label of the input MTZ file"""
        return self._free

    @property
    def mr_python_module(self):
        """The MR python module"""
        if self.mr_program == "molrep":
            return "simbad.mr.molrep_mr"
        elif self.mr_program == "phaser":
            return "simbad.mr.phaser_mr"

    @property
    def mr_program(self):
        """The molecular replacement program to use"""
        return self._mr_program

    @mr_program.setter
    def mr_program(self, mr_program):
        """Define the molecular replacement program to use"""
        if mr_program.lower() in KNOWN_MR_PROGRAMS:
            self._mr_program = mr_program.lower()
        else:
            msg = "Unknown MR program!"
            raise RuntimeError(msg)

    @property
    def refine_python_module(self):
        """The Refinement python module"""
        if self.refine_program == "refmac5":
            return "simbad.mr.refmac_refine"

    @property
    def refine_program(self):
        """The refinement program to use"""
        return self._refine_program

    @refine_program.setter
    def refine_program(self, refine_program):
        """Define the refinement program to use"""
        if refine_program.lower() in KNOWN_REF_PROGRAMS:
            self._refine_program = refine_program
        else:
            msg = "Unknown Refinement program!"
            raise RuntimeError(msg)

    @property
    def output_dir(self):
        """The path to the output directory"""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        """Define the output directory"""
        self._output_dir = output_dir

    @property
    def timeout(self):
        """The time in minutes before phaser is killed"""
        return self._timeout

    @timeout.setter
    def timeout(self, timeout):
        """Define the time in minutes before phaser should be killed"""
        self._timeout = timeout

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
        self._space_group, _, cell_parameters = mtz_util.crystal_data(mtz)
        self._cell_parameters = " ".join(map(str, cell_parameters))

        # Extract column labels from input mtz
        self._f, self._sigf, self._i, self._sigi, self._dano, self._sigdano, self._free = mtz_util.get_labels(
            mtz)

    def submit_jobs(self, results, nproc=1, process_all=False, submit_qtype=None, submit_queue=False, monitor=None):
        """Submit jobs to run in serial or on a cluster

        Parameters
        ----------
        results : class
            Results from :obj: '_LatticeParameterScore' or :obj: '_AmoreRotationScore'
        time_out : int, optional
            Number of seconds for job to timeout [default: 60]
        nproc : int, optional
            Number of processors to use [default: 2]
        process_all : bool, optional
            Terminate MR after a success [default: True]
        submit_qtype : str
            The cluster submission queue type - currently support SGE and LSF
        submit_queue : str
            The queue to submit to on the cluster
        monitor : str

        Returns
        -------
        file
            Output pdb from mr
        file
            Output hkl from mr - if using phaser
        file
            Output log file from mr program
        file
            Output pdb from refinement
        file
            Output hkl from refinement
        file
            Output log file from refinement program

        """
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

        run_files = []
        for result in results:
            mr_pdbin = result.pdb_path
            mr_workdir = os.path.join(
                self.output_dir, result.pdb_code, 'mr', self.mr_program)
            mr_logfile = os.path.join(
                mr_workdir, '{0}_mr.log'.format(result.pdb_code))
            mr_pdbout = os.path.join(
                mr_workdir, '{0}_mr_output.pdb'.format(result.pdb_code))

            ref_workdir = os.path.join(mr_workdir, 'refine')
            ref_hklout = os.path.join(
                ref_workdir, '{0}_refinement_output.mtz'.format(result.pdb_code))
            ref_logfile = os.path.join(
                ref_workdir, '{0}_ref.log'.format(result.pdb_code))
            ref_pdbout = os.path.join(
                ref_workdir, '{0}_refinement_output.pdb'.format(result.pdb_code))

            diff_mapout1 = os.path.join(
                ref_workdir, '{0}_refmac_2fofcwt.map'.format(result.pdb_code))
            diff_mapout2 = os.path.join(
                ref_workdir, '{0}_refmac_fofcwt.map'.format(result.pdb_code))

            # See if the input model can fit in the unit cell, otherwise move to a single chain
            solvent = matthews_coef.solvent_content(mr_pdbin, self.cell_parameters, self.space_group)
            if solvent > 30:
                n_copies = 1
            else:
                pdb_edit.to_single_chain(mr_pdbin, mr_pdbin)
                nres = pdb_edit.number_of_residues(mr_pdbin)
                solvent, n_copies = matthews_coef.matthews_coef(
                    self.cell_parameters, self.space_group, nres)
                msg = "{0} is predicted to be too large to fit in the unit cell with a solvent content " \
                      "of at least 30%, therefore MR will use only the first chain".format(
                          result.pdb_code)
                logger.debug(msg)

            # Common MR keywords
            mr_cmd = ["ccp4-python", "-m", self.mr_python_module, "-enant", self.enant, "-hklin", self.mtz,
                      "-pdbin", mr_pdbin, "-pdbout", mr_pdbout, "-logfile", mr_logfile, "-work_dir", mr_workdir,
                      "-nmol", n_copies]

            # Common refine keywords
            ref_cmd = ["ccp4-python", "-m", self.refine_python_module, "-hklout", ref_hklout, "-pdbin", mr_pdbout,
                       "-pdbout", ref_pdbout, "-logfile", ref_logfile, "-work_dir", ref_workdir]

            # Extend commands with program-specific options
            if self.mr_program == "molrep":
                mr_cmd += ["-space_group", self.space_group]
                ref_cmd += ["-hklin", self.mtz]

            elif self.mr_program == "phaser":
                hklout = os.path.join(
                    mr_workdir, '{0}_mr_output.mtz'.format(result.pdb_code))
                mr_cmd += [
                    "-i", self.i,
                    "-hklout", hklout,
                    "-sigi", self.sigi,
                    "-solvent", solvent,
                    "-timeout", self.timeout,
                ]
                ref_cmd += ["-hklin", hklout]

            # ====
            # Create a run script - prefix __needs__ to contain mr_program so we can find log
            # Leave order of this as SGE does not like scripts with numbers as first char
            # ====
            prefix, stem = self.mr_program + "_", result.pdb_code

            # Set stdin 1
            fft_cmd1, fft_stdin1 = self.fft(
                ref_hklout, diff_mapout1, "2mfo-dfc")
            run_stdin_1 = tmp_file(
                directory=self.output_dir, prefix=prefix, stem=stem, suffix="_1.stdin")
            with open(run_stdin_1, 'w') as f_out:
                f_out.write(fft_stdin1)

            # Set up stdin 2
            fft_cmd2, fft_stdin2 = self.fft(
                ref_hklout, diff_mapout2, "mfo-dfc")
            run_stdin_2 = tmp_file(
                directory=self.output_dir, prefix=prefix, stem=stem, suffix="_2.stdin")
            with open(run_stdin_2, 'w') as f_out:
                f_out.write(fft_stdin2)

            # Set up script and log
            run_script = make_script(
                [mr_cmd,
                 ref_cmd,
                 fft_cmd1 + ["<", run_stdin_1],
                 fft_cmd2 + ["<", run_stdin_2]],
                directory=self.output_dir, prefix=prefix, stem=stem
            )
            run_log = run_script.rsplit(".", 1)[0] + '.log'

            # Save a copy of the files we need to run
            run_files += [(run_script, run_stdin_1, run_stdin_2,
                           run_log, mr_pdbout, mr_logfile, ref_logfile)]

        logger.info("Running %s Molecular Replacement", self.mr_program)
        # Extract relevant data
        run_scripts, _, _, _, mr_pdbouts, mr_logfiles, ref_logfiles = zip(
            *run_files)

        # Execute the scripts
        j = Job(submit_qtype)
        j.submit(run_scripts, directory=self.output_dir, nproc=nproc, name='simbad_mr',
                 submit_queue=submit_queue, permit_nonzero=True)
        # This way we can accept booleans and strings
        interval = int(numpy.log(len(run_scripts)) / 3)
        interval_in_seconds = interval if interval >= 5 else 5
        if process_all:
            j.wait(interval=interval_in_seconds, monitor=monitor)
        else:
            j.wait(interval=interval_in_seconds, monitor=monitor,
                   check_success=mr_succeeded_log)

        # Go through the result and see what has worked
        mr_results = []
        for result, mr_logfile, mr_pdbout, ref_logfile in zip(results, mr_logfiles, mr_pdbouts, ref_logfiles):

            # Quick check to see it's run
            if not os.path.isfile(mr_logfile):
                logger.debug("Cannot find %s MR log file: %s",
                             self.mr_program, mr_logfile)
                continue
            elif not os.path.isfile(ref_logfile):
                logger.debug("Cannot find %s refine log file: %s",
                             self.mr_program, ref_logfile)
                continue
            elif not os.path.isfile(mr_pdbout):
                logger.debug("Cannot find %s output file: %s",
                             self.mr_program, mr_pdbout)
                continue

            # Define a new store container
            score = _MrScore(pdb_code=result.pdb_code)

            # Decide which score to take
            if self.mr_program == "molrep":
                mp = molrep_parser.MolrepParser(mr_logfile)
                score.molrep_score = mp.score
                score.molrep_tfscore = mp.tfscore
            elif self.mr_program == "phaser":
                pp = phaser_parser.PhaserParser(mr_logfile)
                score.phaser_tfz = pp.tfz
                score.phaser_llg = pp.llg
                score.phaser_rfz = pp.rfz

            if self._dano is not None:
                try:
                    anoms = anomalous_util.AnomSearch(
                        self.mtz, self.output_dir, self.mr_program)
                    anoms.run(result)
                    a = anoms.search_results()
                    score.peaks_over_6_rms = a.peaks_over_6_rms
                    score.peaks_over_6_rms_within_4a_of_model = a.peaks_over_6_rms_within_4a_of_model
                    score.peaks_over_9_rms = a.peaks_over_9_rms
                    score.peaks_over_9_rms_within_4a_of_model = a.peaks_over_9_rms_within_4a_of_model
                except RuntimeError:
                    logger.debug(
                        "RuntimeError: Unable to create phased anomalous fourier map for: %s", result.pdb_code)
                except IOError:
                    logger.debug(
                        "IOError: Unable to create phased anomalous fourier map for: %s", result.pdb_code)

            # Analyse the refinement log file
            if os.path.isfile(ref_logfile):
                rp = refmac_parser.RefmacParser(ref_logfile)
                score.final_r_free = rp.final_r_free
                score.final_r_fact = rp.final_r_fact
            else:
                logger.debug("Cannot find %s log file: %s",
                             self.refine_program, ref_logfile)

            # Store this score container information
            mr_results += [score]

        # Save the results
        self._search_results = mr_results

        return

    @staticmethod
    def fft(hklin, mapout, map_type):
        """Function to run fft to generate difference maps for uglymol

        Parameters
        ----------
        hklin : str
           Path to input HKL file
        mapout : str
           Path to output MAP file
        map_type : str
           Define type of run, either mfo-dfc or 2mfo-dfc

        Returns
        -------
        list
           cmd
        str
           stdin

        Raises
        ------
        ValueError
           Unknown map type

        """

        cmd = ["fft", "hklin", hklin, "mapout", mapout]
        if map_type == "2mfo-dfc":
            stdin = "title Sigmaa style 2mfo-dfc map calculated with refmac coefficients" + os.linesep \
                    + "labi F1=FWT PHI=PHWT" + os.linesep \
                    + "end" + os.linesep
        elif map_type == "mfo-dfc":
            stdin = "title Sigmaa style mfo-dfc map calculated with refmac coefficients" + os.linesep \
                    + "labi F1=DELFWT PHI=PHDELWT" + os.linesep \
                    + "end" + os.linesep
        else:
            msg = "Unknown map type!"
            raise ValueError(msg)
        return cmd, stdin

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
        from simbad.util import summarize_result

        columns = []
        if self.mr_program == "molrep":
            columns += ["molrep_score", "molrep_tfscore"]

        elif self.mr_program == "phaser":
            columns += ["phaser_tfz", "phaser_llg", "phaser_rfz"]

        columns += ["final_r_fact", "final_r_free"]

        if self._dano:
            columns += ["peaks_over_6_rms", "peaks_over_6_rms_within_4a_of_model",
                        "peaks_over_9_rms", "peaks_over_9_rms_within_4a_of_model"]

        summarize_result(self.search_results,
                         csv_file=csv_file, columns=columns)


def _mr_job_succeeded(r_fact, r_free):
    """Check values for job success"""
    return r_fact < 0.45 and r_free < 0.45


def mr_succeeded_log(log):
    """Check a Molecular Replacement job for it's success
    
    Parameters
    ----------
    log : str
       The path to a log file

    Returns
    -------
    bool
       Success status of the MR run

    """
    mr_prog, pdb = os.path.basename(log).replace('.log', '').split('_', 1)
    refmac_log = os.path.join(os.path.dirname(
        log), pdb, "mr", mr_prog, "refine", pdb + "_ref.log")
    if os.path.isfile(refmac_log):
        rp = refmac_parser.RefmacParser(refmac_log)
        return _mr_job_succeeded(rp.final_r_fact, rp.final_r_free)
    return False


def mr_succeeded_csvfile(f):
    """Check a Molecular Replacement job for it's success
    
    Parameters
    ----------
    f : str
        The path to f
        
    Returns
    -------
    bool
       Success status of the MR run

    """
    import pandas as pd
    df = pd.read_csv(f)
    data = zip(df.final_r_fact.tolist(), df.final_r_free.tolist())
    return any(_mr_job_succeeded(rfact, rfree) for rfact, rfree in data)
