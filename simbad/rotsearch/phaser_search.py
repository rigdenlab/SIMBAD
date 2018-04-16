"""Module to run the phaser rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "15 April 2018"
__version__ = "0.4"

import logging
import os
import shutil
import uuid
logger = logging.getLogger(__name__)

import pyjob.misc

import simbad.db
import simbad.mr
import simbad.rotsearch
import simbad.score.dat_info
import simbad.score.phaser_score
import simbad.parsers.refmac_parser
import simbad.parsers.rotsearch_parser
import simbad.util.pdb_util
import simbad.util.mtz_util
import simbad.util.matthews_prob

from phaser import InputMR_DAT, runMR_DAT, InputCCA, runCCA
from simbad.util import EXPORT, CMD_PREFIX


class PhaserRotationSearch(object):
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
    Examples
    --------
    >>> from simbad.rotsearch.phaser_search import PhaserRotationSearch
    >>> rotation_search = PhaserRotationSearch('<mtz>', '<mr_program>', '<tmp_dir>', '<work_dir>', '<max_to_keep>',
    ...                                        '<skip_mr>')
    >>> rotation_search.run(
    ...     '<models_dir>', '<nproc>', '<min_solvent_content>', '<submit_qtype>',
    ...     '<submit_queue>', '<monitor>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results
    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from phaser.
    """

    def __init__(self, mtz, mr_program, tmp_dir, work_dir, max_to_keep=20, skip_mr=False, **kwargs):
        self.max_to_keep = max_to_keep
        self.mr_program = mr_program
        self.mtz = mtz
        self.tmp_dir = tmp_dir
        self.work_dir = work_dir
        self.skip_mr = skip_mr

        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.simbad_dat_files = None
        self.submit_qtype = None
        self.submit_queue = None
        self._search_results = None
        self.tested = []

    def run(self, models_dir, nproc=2, min_solvent_content=20, submit_qtype=None,
            submit_queue=None, monitor=None, chunk_size=0, **kwargs):
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
        monitor
        chunk_size : int, optional
            The number of jobs to submit at the same time
        Returns
        -------
        file
            log file for each model in the models_dir
        """
        self.submit_qtype = submit_qtype
        self.submit_queue = submit_queue
        self.f, self.sigf, self.i, self.sigi, _, _, _ = simbad.util.mtz_util.get_labels(self.mtz)

        self.simbad_dat_files = simbad.db.find_simbad_dat_files(models_dir)
        n_files = len(self.simbad_dat_files)

        i = InputMR_DAT()
        i.setHKLI(self.mtz)
        i.setMUTE(True)
        run_mr_data = runMR_DAT(i)

        sg = run_mr_data.getSpaceGroupName().replace(" ", "")
        cell = " ".join(map(str, run_mr_data.getUnitCell()))

        chunk_size = simbad.rotsearch.get_chunk_size(n_files, chunk_size)
        total_chunk_cycles = simbad.rotsearch.get_total_chunk_cycles(n_files, chunk_size)

        mat_coef = simbad.util.matthews_prob.MatthewsProbability(cell, sg)

        dir_name = "simbad-tmp-" + str(uuid.uuid1())
        script_log_dir = os.path.join(self.work_dir, dir_name)
        os.mkdir(script_log_dir)

        ccp4_scr = os.environ["CCP4_SCR"]
        default_tmp_dir = os.path.join(self.work_dir, 'tmp')
        if self.tmp_dir:
            template_tmp_dir = os.path.join(self.tmp_dir, dir_name + "-{0}")
        else:
            template_tmp_dir = os.path.join(default_tmp_dir, dir_name + "-{0}")

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
            pdb_struct = simbad.util.pdb_util.PdbStructure()
            pdb_struct.from_file(dat_model)
            solvent_fraction, n_copies = mat_coef.calculate_content_ncopies_from_struct(pdb_struct)
            solvent_content = solvent_fraction * 100
            if solvent_content < min_solvent_content:
                msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                logger.debug(msg, name, min_solvent_content)
                continue
            model_molecular_weight = pdb_struct.molecular_weight
            mw_diff = abs(predicted_molecular_weight - model_molecular_weight)

            info = simbad.score.dat_info.DatModelInfo(
                name, dat_model, mw_diff, None, None, None, None, solvent_fraction, n_copies
            )
            dat_models.append(info)

        sorted_dat_models = sorted(dat_models, key=lambda x: float(x.mw_diff), reverse=False)

        iteration_range = range(0, n_files, chunk_size)
        for cycle, i in enumerate(iteration_range):
            logger.info("Working on chunk %d out of %d", cycle + 1,
                        total_chunk_cycles)

            template_model = os.path.join("$CCP4_SCR", "{0}.pdb")

            phaser_files = []
            for dat_model in sorted_dat_models[i:i + chunk_size]:
                logger.debug("Generating script to perform PHASER rotation "
                             + "function on %s", dat_model.pdb_code)

                pdb_model = template_model.format(dat_model.pdb_code)
                template_rot_log = os.path.join("$CCP4_SCR", "{0}_rot.log")

                conv_py = "\"from simbad.db import convert_dat_to_pdb; convert_dat_to_pdb('{}', '{}')\""
                conv_py = conv_py.format(dat_model.dat_path, pdb_model)

                rot_log = template_rot_log.format(dat_model.pdb_code)
                tmp_dir = template_tmp_dir.format(dat_model.pdb_code)

                phaser_cmd = ["simbad.rotsearch.phaser_rotation_search",
                              "-hklin", self.mtz,
                              "-f", self.f,
                              "-sigf", self.sigf,
                              "-i", self.i,
                              "-sigi", self.sigi,
                              "-pdbin", pdb_model,
                              "-logfile", rot_log,
                              "-solvent", dat_model.solvent,
                              "-nmol", dat_model.nmol,
                              "-work_dir", tmp_dir
                              ]
                phaser_cmd = " ".join(str(e) for e in phaser_cmd)

                cmd = [
                    [EXPORT, "CCP4_SCR=" + tmp_dir],
                    ["mkdir", "-p", "$CCP4_SCR\n"],
                    [CMD_PREFIX, "$CCP4/bin/ccp4-python", "-c", conv_py, os.linesep],
                    [CMD_PREFIX, "$CCP4/bin/ccp4-python", "-m", phaser_cmd, os.linesep],
                    ["rm", "-rf", "$CCP4_SCR\n"],
                    [EXPORT, "CCP4_SCR=" + ccp4_scr],
                ]
                phaser_script = pyjob.misc.make_script(
                    cmd, directory=script_log_dir, prefix="phaser_", stem=dat_model.pdb_code
                )
                phaser_log = phaser_script.rsplit(".", 1)[0] + '.log'
                phaser_files += [(phaser_script, phaser_log, dat_model.dat_path)]

            results = []
            if len(phaser_files) > 0:
                logger.info("Running PHASER rotation functions")
                phaser_scripts, phaser_logs, dat_models = zip(*phaser_files)
                simbad.rotsearch.submit_chunk(phaser_scripts, script_log_dir, nproc, 'simbad_phaser',
                                              submit_qtype, submit_queue, monitor, self.rot_succeeded_log)

                for dat_model, phaser_log in zip(dat_models, phaser_logs):
                    base = os.path.basename(phaser_log)
                    pdb_code = base.replace("phaser_", "").replace(".log", "")
                    phaser_rotation_parser = simbad.parsers.rotsearch_parser.PhaserRotsearchParser(
                        phaser_log
                    )
                    score = simbad.score.phaser_score.PhaserRotationScore(pdb_code, dat_model,
                                                                          phaser_rotation_parser.llg,
                                                                          phaser_rotation_parser.rfz)

                    if phaser_rotation_parser.rfz:
                        results += [score]

            else:
                logger.critical("No structures to be trialled")

            self._search_results = results
            shutil.rmtree(script_log_dir)

            if os.path.isdir(default_tmp_dir):
                shutil.rmtree(default_tmp_dir)

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
        columns = [
            "llg", "rfz"
        ]
        summarize_result(self.search_results,
                         csv_file=csv_file, columns=columns)

    @property
    def search_results(self):
        return sorted(self._search_results, key=lambda x: float(x.llg), reverse=True)[:self.max_to_keep]

    @staticmethod
    def _rot_job_succeeded(phaser_llg_score):
        """Check values for job success"""
        return phaser_llg_score > 100

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
        if self.skip_mr:
            return False

        rot_prog, pdb = os.path.basename(log).replace('.log', '').split('_', 1)
        rotsearch_parser = simbad.parsers.rotsearch_parser.PhaserRotsearchParser(log)
        dat_model = [s for s in self.simbad_dat_files if pdb in s][0]
        score = simbad.score.phaser_score.PhaserRotationScore(
            pdb, dat_model, rotsearch_parser.llg, rotsearch_parser.rfz
        )
        results = [score]
        if self._rot_job_succeeded(rotsearch_parser.llg) and pdb not in self.tested:
            self.tested.append(pdb)
            output_dir = os.path.join(self.work_dir, "mr_search")
            mr = simbad.mr.MrSubmit(mtz=self.mtz,
                                    mr_program=self.mr_program,
                                    refine_program='refmac5',
                                    refine_type=None,
                                    refine_cycles=0,
                                    output_dir=output_dir,
                                    tmp_dir=self.tmp_dir,
                                    timeout=30)
            mr.mute = True
            mr.submit_jobs(results,
                           nproc=1,
                           process_all=True,
                           submit_qtype=self.submit_qtype,
                           submit_queue=self.submit_queue)
            refmac_log = os.path.join(output_dir, pdb, "mr", self.mr_program, "refine", pdb + "_ref.log")
            if os.path.isfile(refmac_log):
                refmac_parser = simbad.parsers.refmac_parser.RefmacParser(refmac_log)
                return simbad.rotsearch.mr_job_succeeded(refmac_parser.final_r_fact, refmac_parser.final_r_free)
        return False