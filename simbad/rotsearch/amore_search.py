"""Module to run the amore rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "15 April 2018"
__version__ = "0.4"

import glob
import logging
import os
import shutil
import uuid

logger = logging.getLogger(__name__)

import pyjob
from pyjob.script import ScriptCollector, Script

import simbad.db
import simbad.mr
import simbad.rotsearch
import simbad.core.amore_score
import simbad.core.dat_score
import simbad.parsers.phaser_parser
import simbad.parsers.refmac_parser
import simbad.parsers.rotsearch_parser
import simbad.util
import simbad.util.pdb_util
import simbad.util.mtz_util
import simbad.util.matthews_prob

from simbad.util import EXPORT, CMD_PREFIX, CCP4_SOURCE, CCP4_SCRATCH, MKDIR_CMD, RM_CMD


class AmoreRotationSearch(simbad.rotsearch._RotationSearch):
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
    search_results : list
        The search results

    Examples
    --------
    >>> from simbad.rotsearch.amore_search import AmoreRotationSearch
    >>> rotation_search = AmoreRotationSearch('<mtz>', '<mr_program>', '<tmp_dir>', '<work_dir>',
    ...                                       '<amore_exe>', '<max_to_keep>', '<skip_mr>', '<process_all>')
    >>> rotation_search.run(
    ...     '<models_dir>', '<output_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>',
    ...     '<min_solvent_content>', '<submit_nproc>', '<submit_qtype>', '<submit_queue>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results


    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from amore.

    """

    def __init__(self, mtz, mr_program, tmp_dir, work_dir, amore_exe=None, max_to_keep=20, skip_mr=False, process_all=False, **kwargs):

        super(AmoreRotationSearch, self).__init__(mtz, mr_program, tmp_dir, work_dir,
                                                  max_to_keep=max_to_keep, skip_mr=skip_mr, process_all=process_all)

        self.amore_exe = amore_exe

        self.hklpck0 = None
        self.shres = None
        self.pklim = None
        self.progress = -5
        self.npic = None
        self.rotastep = None
        self.ccp4_scr = None
        self.script_log_dir = None

        self.columns = [
            "ALPHA",
            "BETA",
            "GAMMA",
            "CC_F",
            "RF_F",
            "CC_I",
            "CC_P",
            "Icp",
            "CC_F_Z_score",
            "CC_P_Z_score",
            "Number_of_rotation_searches_producing_peak",
        ]
        self.score_column = "CC_F_Z_score"

        self.template_hklpck1 = os.path.join(CCP4_SCRATCH, "{0}.hkl")
        self.template_clmn0 = os.path.join(CCP4_SCRATCH, "{0}_spmipch.clmn")
        self.template_clmn1 = os.path.join(CCP4_SCRATCH, "{0}.clmn")
        self.template_mapout = os.path.join(CCP4_SCRATCH, "{0}_amore_cross.map")
        self.template_table1 = os.path.join(CCP4_SCRATCH, "{0}_sfs.tab")
        self.template_model = os.path.join(CCP4_SCRATCH, "{0}.pdb")
        self.template_rot_log = os.path.join(CCP4_SCRATCH, "{0}_rot.log")
        self.template_tmp_dir = None

    def run(
        self,
        models_dir,
        nproc=2,
        shres=3.0,
        pklim=0.5,
        npic=50,
        rotastep=1.0,
        min_solvent_content=20,
        submit_qtype=None,
        submit_queue=None,
        chunk_size=0,
        **kwargs
    ):
        """Run amore rotation function on a directory of models

        Parameters
        ----------
        models_dir : str
            The directory containing the models to run the rotation search on
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
        submit_qtype : str
            The cluster submission queue type - currently support SGE and LSF
        submit_queue : str
            The queue to submit to on the cluster
        chunk_size : int, optional
            The number of jobs to submit at the same time

        Returns
        -------
        file
            log file for each model in the models_dir

        """
        from phaser import InputMR_DAT, runMR_DAT, InputCCA, runCCA

        self.shres = shres
        self.pklim = pklim
        self.npic = npic
        self.rotastep = rotastep

        self.submit_qtype = submit_qtype
        self.submit_queue = submit_queue

        self.simbad_dat_files = simbad.db.find_simbad_dat_files(models_dir)

        i = InputMR_DAT()
        i.setHKLI(self.mtz)
        i.setMUTE(True)
        run_mr_data = runMR_DAT(i)

        sol_calc = simbad.util.matthews_prob.SolventContent(self.mtz_obj.cell.volume_per_image())

        dir_name = "simbad-tmp-" + str(uuid.uuid1())
        self.script_log_dir = os.path.join(self.work_dir, dir_name)
        os.mkdir(self.script_log_dir)

        self.hklpck0 = self._generate_hklpck0()

        self.ccp4_scr = os.environ["CCP4_SCR"]
        default_tmp_dir = os.path.join(self.work_dir, "tmp")
        if self.tmp_dir:
            self.template_tmp_dir = os.path.join(self.tmp_dir, dir_name + "-{0}")
        else:
            self.template_tmp_dir = os.path.join(default_tmp_dir, dir_name + "-{0}")

        predicted_molecular_weight = 0
        if run_mr_data.Success():
            i = InputCCA()
            i.setSPAC_HALL(run_mr_data.getSpaceGroupHall())
            i.setCELL6(run_mr_data.getUnitCell())
            i.setMUTE(True)
            run_cca = runCCA(i)

            if run_cca.Success():
                predicted_molecular_weight = run_cca.getAssemblyMW()

        dat_models = []
        for dat_model in self.simbad_dat_files:
            name = os.path.basename(dat_model.replace(".dat", ""))
            try:
                pdb_struct = simbad.util.pdb_util.PdbStructure.from_file(dat_model)
            except Exception: # Catch all issues here
                msg = "Skipping %s: Problem with dat file"
                logger.debug(msg, name)
                continue
            solvent_fraction = sol_calc.calculate_from_struct(pdb_struct)
            solvent_content = solvent_fraction * 100
            if solvent_content < min_solvent_content:
                msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                logger.debug(msg, name, min_solvent_content)
                continue
            x, y, z, intrad = pdb_struct.integration_box
            model_molecular_weight = pdb_struct.molecular_weight
            mw_diff = abs(predicted_molecular_weight - model_molecular_weight)

            info = simbad.core.dat_score.DatModelScore(name, dat_model, mw_diff, x, y, z, intrad, solvent_content, None)
            dat_models.append(info)

        sorted_dat_models = sorted(dat_models, key=lambda x: float(x.mw_diff), reverse=False)
        n_files = len(sorted_dat_models)
        chunk_size = simbad.rotsearch.get_chunk_size(n_files, chunk_size)
        total_chunk_cycles = simbad.rotsearch.get_total_chunk_cycles(n_files, chunk_size)

        results = []
        iteration_range = range(0, n_files, chunk_size)
        for cycle, i in enumerate(iteration_range):
            logger.info("Working on chunk %d out of %d", cycle + 1, total_chunk_cycles)

            if self.solution:
                logger.info("Early termination criteria met, skipping chunk %d", cycle + 1)
                continue

            collector = ScriptCollector(None)
            amore_files = []
            for dat_model in sorted_dat_models[i: i + chunk_size]:
                script, run_file = self.generate_script(dat_model)
                collector.add(script)
                amore_files.append(run_file)

            if len(collector.scripts) > 0:
                logger.info("Running AMORE tab/rot functions")
                amore_logs, dat_models = zip(*amore_files)
                simbad.util.submit_chunk(
                    collector, self.script_log_dir, nproc, "simbad_amore", submit_qtype, submit_queue, True,
                    self.progress_monitor,
                    self.rot_succeeded_log
                )

                for dat_model, amore_log in zip(dat_models, amore_logs):
                    base = os.path.basename(amore_log)
                    pdb_code = base.replace("amore_", "").replace(".log", "")
                    try:
                        rotsearch_parser = simbad.parsers.rotsearch_parser.AmoreRotsearchParser(amore_log)
                        score = simbad.core.amore_score.AmoreRotationScore(
                            pdb_code,
                            dat_model,
                            rotsearch_parser.alpha,
                            rotsearch_parser.beta,
                            rotsearch_parser.gamma,
                            rotsearch_parser.cc_f,
                            rotsearch_parser.rf_f,
                            rotsearch_parser.cc_i,
                            rotsearch_parser.cc_p,
                            rotsearch_parser.icp,
                            rotsearch_parser.cc_f_z_score,
                            rotsearch_parser.cc_p_z_score,
                            rotsearch_parser.num_of_rot,
                        )
                        if rotsearch_parser.cc_f_z_score:
                            results += [score]
                    except IOError:
                        pass

            else:
                logger.critical("No structures to be trialled")

        self._search_results = results
        shutil.rmtree(self.script_log_dir)

        if os.path.isdir(default_tmp_dir):
            shutil.rmtree(default_tmp_dir)

    def generate_script(self, dat_model):
        logger.debug("Generating script to perform AMORE rotation " + "function on %s", dat_model.pdb_code)

        pdb_model = self.template_model.format(dat_model.pdb_code)
        table1 = self.template_table1.format(dat_model.pdb_code)
        hklpck1 = self.template_hklpck1.format(dat_model.pdb_code)
        clmn0 = self.template_clmn0.format(dat_model.pdb_code)
        clmn1 = self.template_clmn1.format(dat_model.pdb_code)
        mapout = self.template_mapout.format(dat_model.pdb_code)

        conv_py = "\"from simbad.db import convert_dat_to_pdb; convert_dat_to_pdb(r'{}', r'{}')\""
        conv_py = conv_py.format(dat_model.dat_path, pdb_model)

        tab_cmd = [self.amore_exe, "xyzin1", pdb_model, "xyzout1", pdb_model, "table1", table1]
        tab_stdin = self.tabfun_stdin_template.format(x=dat_model.x, y=dat_model.y, z=dat_model.z, a=90, b=90, c=120)

        rot_cmd = [self.amore_exe, "table1", table1, "HKLPCK1", hklpck1, "hklpck0", self.hklpck0, "clmn1", clmn1, "clmn0", clmn0, "MAPOUT", mapout]
        rot_stdin = self.rotfun_stdin_template.format(shres=self.shres, intrad=dat_model.intrad, pklim=self.pklim, npic=self.npic, step=self.rotastep)
        rot_log = self.template_rot_log.format(dat_model.pdb_code)

        tmp_dir = self.template_tmp_dir.format(dat_model.pdb_code)

        source = simbad.util.source_ccp4()

        cmd = [
            [source],
            [EXPORT, "CCP4_SCR=" + tmp_dir],
            [MKDIR_CMD, CCP4_SCRATCH, os.linesep],
            [CMD_PREFIX, CCP4_SOURCE + "/bin/ccp4-python", "-c", conv_py, os.linesep],
        ]

        if os.name == "nt":
            tab_stdin_file = os.path.join(CCP4_SCRATCH, 'tab_stdin')
            rot_stdin_file = os.path.join(CCP4_SCRATCH, 'rot_stdin')

            for i, line in enumerate(tab_stdin.split('\n')):
                if i == 0:
                    cmd += [["echo", line, ">", tab_stdin_file]]
                else:
                    cmd += [["echo", line, ">>", tab_stdin_file]]

            for i, line in enumerate(rot_stdin.split('\n')):
                if i == 0:
                    cmd += [["echo", line, ">", rot_stdin_file]]
                else:
                    cmd += [["echo", line, ">>", rot_stdin_file]]

            cmd += [
                tab_cmd + ["<", tab_stdin_file, ">", os.devnull],
                rot_cmd + ["<", rot_stdin_file, ">", rot_log],
                [EXPORT, '"match="', os.linesep],
                ['for /F "delims=" %%G in', "('findstr SOLUTIONRCD", rot_log + "')",
                 'do (if not defined match set "match=%%G" & goto :found)', os.linesep],
                [":found", os.linesep],
                ["echo %match%", os.linesep]
            ]
        else:
            cmd += [
                tab_cmd + ["<< eof >", os.devnull],
                [tab_stdin],
                ["eof"],
                [os.linesep],
                rot_cmd + ["<< eof >", rot_log],
                [rot_stdin],
                ["eof"],
                [os.linesep],
                ["grep", "-m 1", "SOLUTIONRCD", rot_log, os.linesep]
            ]

        cmd += [
            [RM_CMD, CCP4_SCRATCH, os.linesep],
            [EXPORT, "CCP4_SCR=" + self.ccp4_scr],
        ]
        amore_script = Script(directory=self.script_log_dir, prefix="amore_", stem=dat_model.pdb_code)
        for c in cmd:
            amore_script.append(" ".join(map(str, c)))
        amore_log = amore_script.path.rsplit(".", 1)[0] + ".log"
        amore_files = (amore_log, dat_model.dat_path)
        amore_script.write()
        return amore_script, amore_files

    def _generate_hklpck0(self):
        logger.info("Preparing files for AMORE rotation function")
        stdin = self.sortfun_stdin_template.format(f=self.mtz_obj.f, sigf=self.mtz_obj.sigf)
        hklpck0 = os.path.join(self.work_dir, "spmipch.hkl")
        cmd = [self.amore_exe, "hklin", self.mtz, "hklpck0", hklpck0]
        pyjob.cexec(cmd, stdin=stdin)
        return hklpck0

    @property
    def sortfun_stdin_template(self):
        return """TITLE   ** spmi  packing h k l F for crystal**
SORTFUN RESOL 100.  2.5
LABI FP={f}  SIGFP={sigf}"""

    @property
    def tabfun_stdin_template(self):
        return """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL {x} {y} {z} {a} {b} {c} ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0"""

    @property
    def rotfun_stdin_template(self):
        return """TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1
ROTFUN
GENE 1   RESO 100.0 {shres}  CELL_MODEL 80 75 65
CLMN CRYSTAL ORTH  1 RESO  20.0  {shres}  SPHERE   {intrad}
CLMN MODEL 1     RESO  20.0  {shres} SPHERE   {intrad}
ROTA  CROSS  MODEL 1  PKLIM {pklim}  NPIC {npic} STEP {step}"""

    @staticmethod
    def _rot_job_succeeded(amore_z_score):
        """Check values for job success"""
        return amore_z_score > 10

    def rot_succeeded_log(self, log):
        """Check a rotation search job for it's success

        Parameters
        ----------
        log : str
           The path to a log file

        Returns
        -------
        bool
           Success status of the rot run

        """
        if self.skip_mr or self.process_all:
            return False

        rot_prog, pdb = os.path.basename(log).replace(".log", "").split("_", 1)
        rotsearch_parser = simbad.parsers.rotsearch_parser.AmoreRotsearchParser(log)
        dat_model = [s for s in self.simbad_dat_files if pdb in s][0]
        score = simbad.core.amore_score.AmoreRotationScore(
            pdb,
            dat_model,
            rotsearch_parser.alpha,
            rotsearch_parser.beta,
            rotsearch_parser.gamma,
            rotsearch_parser.cc_f,
            rotsearch_parser.rf_f,
            rotsearch_parser.cc_i,
            rotsearch_parser.cc_p,
            rotsearch_parser.icp,
            rotsearch_parser.cc_f_z_score,
            rotsearch_parser.cc_p_z_score,
            rotsearch_parser.num_of_rot,
        )
        results = [score]
        try:
            job_succeeded = self._rot_job_succeeded(rotsearch_parser.cc_f_z_score)
        except TypeError:
            return False

        if job_succeeded and pdb not in self.tested:
            self.tested.append(pdb)
            output_dir = os.path.join(self.work_dir, "mr_search")
            mr = simbad.mr.MrSubmit(
                mtz=self.mtz,
                mr_program=self.mr_program,
                refine_program="refmac5",
                refine_type=None,
                refine_cycles=0,
                output_dir=output_dir,
                sgalternative="none",
                tmp_dir=self.tmp_dir,
                timeout=30,
            )
            mr.mute = True
            mr.submit_jobs(results, nproc=1, process_all=True, submit_qtype=self.submit_qtype, submit_queue=self.submit_queue)
            mr_log = os.path.join(output_dir, pdb, "mr", self.mr_program, pdb + "_mr.log")
            refmac_log = os.path.join(output_dir, pdb, "mr", self.mr_program, "refine", pdb + "_ref.log")
            if os.path.isfile(refmac_log):
                refmac_parser = simbad.parsers.refmac_parser.RefmacParser(refmac_log)
                if simbad.mr._refinement_succeeded(refmac_parser.final_r_fact, refmac_parser.final_r_free):
                    self.solution = True
                    return True
            if os.path.isfile(mr_log):
                if self.mr_program == "phaser":
                    phaser_parser = simbad.parsers.phaser_parser.PhaserParser(mr_log)
                    if simbad.mr._phaser_succeeded(phaser_parser.llg, phaser_parser.tfz):
                        self.solution = True
                        return True
        return False

    def progress_monitor(self):
        total_log_files = 0
        log_files = glob.glob(os.path.join(self.script_log_dir, '*.log'))
        for log in log_files:
            with open(log, 'r') as f:
                total_log_files += sum([1 for line in f.readlines() if " SOLUTIONRCD" in line])
        total_sh_files = len(glob.glob(os.path.join(self.script_log_dir, '*.sh')))
        percentage_complete = (total_log_files / total_sh_files) * 100
        if percentage_complete - self.progress >= 5:
            logger.info("Percentage complete: {:.1f}%".format(percentage_complete))
            self.progress = percentage_complete
