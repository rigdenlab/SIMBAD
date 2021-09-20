"""Module to run the phaser rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "15 April 2018"
__version__ = "0.4"

import glob
import logging
import os
import shutil
import uuid

logger = logging.getLogger(__name__)

from pyjob.script import ScriptCollector, Script

import simbad.db
import simbad.mr
import simbad.rotsearch
import simbad.core.dat_score
import simbad.core.phaser_score
import simbad.parsers.phaser_parser
import simbad.parsers.refmac_parser
import simbad.parsers.rotsearch_parser
import simbad.util
import simbad.util.pdb_util
import simbad.util.matthews_prob

from simbad.util import EXPORT, CMD_PREFIX, CCP4_SOURCE, CCP4_SCRATCH, MKDIR_CMD, RM_CMD


class PhaserRotationSearch(simbad.rotsearch._RotationSearch):
    """A class to perform the phaser rotation search
    Attributes
    ----------
    mtz : str
        The path to the input MTZ
    i : str
        Column label for I
    sigi : str
        Column label for SIGI
    work_dir : str
        The path to the working directory
    max_to_keep : int
        The maximum number of results to keep [default: 20]
    eid : int, optional
            The estimated sequence identity from which to calculate ermsd
    Examples
    --------
    >>> from simbad.rotsearch.phaser_search import PhaserRotationSearch
    >>> rotation_search = PhaserRotationSearch('<mtz>', '<mr_program>', '<tmp_dir>', '<work_dir>', '<max_to_keep>',
    ...                                        '<skip_mr>', '<eid>', '<process_all>')
    >>> rotation_search.run(
    ...     '<models_dir>', '<nproc>', '<min_solvent_content>', '<submit_nproc>', '<submit_qtype>',
    ...     '<submit_queue>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results
    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from phaser.
    """

    def __init__(self, mtz, mr_program, tmp_dir, work_dir, max_to_keep=20, skip_mr=False, eid=70, process_all=False, **kwargs):
        super(PhaserRotationSearch, self).__init__(mtz, mr_program, tmp_dir, work_dir,
                                                   max_to_keep=max_to_keep, skip_mr=skip_mr, process_all=process_all)
        self.eid = eid
        self.ccp4_scr = None
        self.script_log_dir = None

        self.columns = ['llg', 'rfz']
        self.progress = -5
        self.score_column = 'rfz'
        self.template_model = None
        self.template_tmp_dir = None

    def run(
        self,
        models_dir,
        nproc=2,
        min_solvent_content=20,
        submit_qtype=None,
        submit_queue=None,
        chunk_size=0,
        **kwargs
    ):
        """Run phaser rotation function on a directory of models
        Parameters
        ----------
        models_dir : str
            The directory containing the models to run the rotation search on
        nproc : int, optional
            The number of processors to run the job on
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
        self.submit_qtype = submit_qtype
        self.submit_queue = submit_queue

        self.simbad_dat_files = simbad.db.find_simbad_dat_files(models_dir)

        i = InputMR_DAT()
        i.setHKLI(self.mtz)
        i.setMUTE(True)
        run_mr_data = runMR_DAT(i)

        mat_prob = simbad.util.matthews_prob.MatthewsProbability(self.mtz_obj.cell.volume_per_image())

        dir_name = "simbad-tmp-" + str(uuid.uuid1())
        self.script_log_dir = os.path.join(self.work_dir, dir_name)
        os.mkdir(self.script_log_dir)

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
            solvent_fraction, n_copies = mat_prob.calculate_from_struct(pdb_struct)
            solvent_content = solvent_fraction * 100
            if solvent_content < min_solvent_content:
                msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                logger.debug(msg, name, min_solvent_content)
                continue
            mw_diff = abs(predicted_molecular_weight - pdb_struct.molecular_weight)

            info = simbad.core.dat_score.DatModelScore(name, dat_model, mw_diff, None, None, None, None,
                                                       solvent_fraction, n_copies)
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

            self.template_model = os.path.join(CCP4_SCRATCH, "{0}.pdb")

            collector = ScriptCollector(None)
            phaser_files = []
            for dat_model in sorted_dat_models[i: i + chunk_size]:
                script, run_file = self.generate_script(dat_model)
                collector.add(script)
                phaser_files.append(run_file)

            if len(phaser_files) > 0:
                logger.info("Running PHASER rotation functions")
                phaser_logs, dat_models = zip(*phaser_files)
                simbad.util.submit_chunk(
                    collector, self.script_log_dir, nproc, "simbad_phaser", submit_qtype, submit_queue, True,
                    self.progress_monitor,
                    self.rot_succeeded_log
                )

                for dat_model, phaser_log in zip(dat_models, phaser_logs):
                    base = os.path.basename(phaser_log)
                    pdb_code = base.replace("phaser_", "").replace(".log", "")
                    try:
                        phaser_rotation_parser = simbad.parsers.rotsearch_parser.PhaserRotsearchParser(phaser_log)
                        if phaser_rotation_parser.rfact:
                            phaser_rotation_parser.llg = 100
                            phaser_rotation_parser.rfz = 10
                        score = simbad.core.phaser_score.PhaserRotationScore(
                            pdb_code, dat_model, phaser_rotation_parser.llg, phaser_rotation_parser.rfz
                        )

                        if phaser_rotation_parser.rfz:
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
        logger.debug("Generating script to perform PHASER rotation " + "function on %s", dat_model.pdb_code)

        pdb_model = self.template_model.format(dat_model.pdb_code)
        template_rot_log = os.path.join(CCP4_SCRATCH, "{0}_rot.log")

        conv_py = "\"from simbad.db import convert_dat_to_pdb; convert_dat_to_pdb(r'{}', r'{}')\""
        conv_py = conv_py.format(dat_model.dat_path, pdb_model)

        rot_log = template_rot_log.format(dat_model.pdb_code)
        tmp_dir = self.template_tmp_dir.format(dat_model.pdb_code)

        phaser_cmd = [
            "simbad.rotsearch.phaser_rotation_search",
            "-eid",
            self.eid,
            "-hklin",
            self.mtz,
            "-f",
            self.mtz_obj.f,
            "-sigf",
            self.mtz_obj.sigf,
            "-i",
            self.mtz_obj.i,
            "-sigi",
            self.mtz_obj.sigi,
            "-pdbin",
            pdb_model,
            "-logfile",
            rot_log,
            "-solvent",
            dat_model.solvent,
            "-nmol",
            dat_model.nmol,
            "-work_dir",
            tmp_dir,
        ]
        phaser_cmd = " ".join(str(e) for e in phaser_cmd)

        source = simbad.util.source_ccp4()

        cmd = [
            [source],
            [EXPORT, "CCP4_SCR=" + tmp_dir],
            [MKDIR_CMD, CCP4_SCRATCH, os.linesep],
            [CMD_PREFIX, CCP4_SOURCE + "/bin/ccp4-python", "-c", conv_py, os.linesep],
            [CMD_PREFIX, CCP4_SOURCE + "/bin/ccp4-python", "-m", phaser_cmd, os.linesep],
            [RM_CMD, CCP4_SCRATCH, os.linesep],
            [EXPORT, "CCP4_SCR=" + self.ccp4_scr],
        ]

        phaser_script = Script(directory=self.script_log_dir, prefix="phaser_", stem=dat_model.pdb_code)
        for c in cmd:
            phaser_script.append(" ".join(map(str, c)))
        phaser_log = phaser_script.path.rsplit(".", 1)[0] + ".log"
        phaser_files = (phaser_log, dat_model.dat_path)
        phaser_script.write()
        return phaser_script, phaser_files

    @staticmethod
    def _rot_job_succeeded(phaser_rfz_score):
        """Check values for job success"""
        return phaser_rfz_score > 7

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
        rotsearch_parser = simbad.parsers.rotsearch_parser.PhaserRotsearchParser(log)
        dat_model = [s for s in self.simbad_dat_files if pdb in s][0]
        score = simbad.core.phaser_score.PhaserRotationScore(pdb, dat_model, rotsearch_parser.llg, rotsearch_parser.rfz)
        results = [score]
        if self._rot_job_succeeded(rotsearch_parser.rfz) or rotsearch_parser.rfact:
            if pdb not in self.tested:
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
                mr.submit_jobs(results, nproc=1, process_all=True, submit_qtype=self.submit_qtype,
                               submit_queue=self.submit_queue)
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
                total_log_files += sum([1 for line in f.readlines() if "EXIT STATUS: SUCCESS" in line])
        total_sh_files = len(glob.glob(os.path.join(self.script_log_dir, '*.sh')))
        percentage_complete = (total_log_files / total_sh_files) * 100
        if percentage_complete - self.progress >= 5:
            logger.info("Percentage complete: {:.1f}%".format(percentage_complete))
            self.progress = percentage_complete
