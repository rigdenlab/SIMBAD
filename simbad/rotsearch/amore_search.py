"""Module to run the amore rotation search"""

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
import simbad.core.amore_score
import simbad.core.dat_score
import simbad.parsers.refmac_parser
import simbad.parsers.rotsearch_parser
import simbad.util.pdb_util
import simbad.util.mtz_util
import simbad.util.matthews_prob

from phaser import InputMR_DAT, runMR_DAT, InputCCA, runCCA
from simbad.util import EXPORT, CMD_PREFIX


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
    search_results : list
        The search results

    Examples
    --------
    >>> from simbad.rotsearch.amore_search import AmoreRotationSearch
    >>> rotation_search = AmoreRotationSearch('<amore_exe>', '<mtz>', '<mr_program>', '<tmp_dir>', '<work_dir>',
    ...                                       '<max_to_keep>')
    >>> rotation_search.run(
    ...     '<models_dir>', '<output_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>',
    ...     '<min_solvent_content>', '<submit_qtype>', '<submit_queue>', '<monitor>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results


    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from amore.

    """

    def __init__(self, mtz, mr_program, tmp_dir, work_dir, amore_exe=None, max_to_keep=20, skip_mr=False, **kwargs):

        self.amore_exe = amore_exe
        self.max_to_keep = max_to_keep
        self.mr_program = mr_program
        self.mtz = mtz
        self.skip_mr = skip_mr
        self.tmp_dir = tmp_dir
        self.work_dir = work_dir

        self.f = None
        self.sigf = None
        self.simbad_dat_files = None
        self.submit_qtype = None
        self.submit_queue = None
        self._search_results = None
        self.tested = []

    def run(self, models_dir, nproc=2, shres=3.0, pklim=0.5, npic=50, rotastep=1.0, min_solvent_content=20,
            submit_qtype=None, submit_queue=None, monitor=None, chunk_size=0, **kwargs):
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

        sol_calc = simbad.util.matthews_prob.SolventContent(cell, sg)

        dir_name = "simbad-tmp-" + str(uuid.uuid1())
        script_log_dir = os.path.join(self.work_dir, dir_name)
        os.mkdir(script_log_dir)

        hklpck0 = self._generate_hklpck0()

        ccp4_scr = os.environ["CCP4_SCR"]
        default_tmp_dir = os.path.join(self.work_dir, 'tmp')
        if self.tmp_dir:
            template_tmp_dir = os.path.join(self.tmp_dir, dir_name + "-{0}")
        else:
            template_tmp_dir = os.path.join(default_tmp_dir, dir_name + "-{0}")

        template_hklpck1 = os.path.join("$CCP4_SCR", "{0}.hkl")
        template_clmn0 = os.path.join("$CCP4_SCR", "{0}_spmipch.clmn")
        template_clmn1 = os.path.join("$CCP4_SCR", "{0}.clmn")
        template_mapout = os.path.join("$CCP4_SCR", "{0}_amore_cross.map")
        template_table1 = os.path.join("$CCP4_SCR", "{0}_sfs.tab")
        template_model = os.path.join("$CCP4_SCR", "{0}.pdb")
        template_rot_log = os.path.join("$CCP4_SCR", "{0}_rot.log")

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
            try:
                solvent_content = sol_calc.calculate_from_struct(pdb_struct)
                if solvent_content < min_solvent_content:
                    msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                    logger.debug(msg, name, min_solvent_content)
                    continue
            except ValueError:
                msg = "Skipping %s: Error calculating solvent content"
                logger.debug(msg, name)

            x, y, z, intrad = pdb_struct.integration_box
            model_molecular_weight = pdb_struct.molecular_weight
            mw_diff = abs(predicted_molecular_weight - model_molecular_weight)

            info = simbad.core.dat_score.DatModelScore(
                name, dat_model, mw_diff, x, y, z, intrad, solvent_content, None
            )
            dat_models.append(info)

        sorted_dat_models = sorted(dat_models, key=lambda x: float(x.mw_diff), reverse=False)

        iteration_range = range(0, n_files, chunk_size)
        for cycle, i in enumerate(iteration_range):
            logger.info("Working on chunk %d out of %d", cycle + 1,
                        total_chunk_cycles)

            amore_files = []
            for dat_model in sorted_dat_models[i:i + chunk_size]:
                logger.debug("Generating script to perform AMORE rotation "
                             + "function on %s", dat_model.pdb_code)

                pdb_model = template_model.format(dat_model.pdb_code)
                table1 = template_table1.format(dat_model.pdb_code)
                hklpck1 = template_hklpck1.format(dat_model.pdb_code)
                clmn0 = template_clmn0.format(dat_model.pdb_code)
                clmn1 = template_clmn1.format(dat_model.pdb_code)
                mapout = template_mapout.format(dat_model.pdb_code)

                conv_py = "\"from simbad.db import convert_dat_to_pdb; convert_dat_to_pdb('{}', '{}')\""
                conv_py = conv_py.format(dat_model.dat_path, pdb_model)

                tab_cmd = [self.amore_exe, "xyzin1", pdb_model, "xyzout1",
                           pdb_model, "table1", table1]
                tab_stdin = self.tabfun_stdin_template.format(
                    x=dat_model.x, y=dat_model.y, z=dat_model.z, a=90, b=90, c=120
                )

                rot_cmd = [self.amore_exe, 'table1', table1, 'HKLPCK1', hklpck1,
                           'hklpck0', hklpck0, 'clmn1', clmn1, 'clmn0', clmn0,
                           'MAPOUT', mapout]
                rot_stdin = self.rotfun_stdin_template.format(
                    shres=shres, intrad=dat_model.intrad, pklim=pklim, npic=npic,
                    step=rotastep
                )
                rot_log = template_rot_log.format(dat_model.pdb_code)

                tmp_dir = template_tmp_dir.format(dat_model.pdb_code)
                cmd = [
                    [EXPORT, "CCP4_SCR=" + tmp_dir],
                    ["mkdir", "-p", "$CCP4_SCR\n"],
                    [CMD_PREFIX, "$CCP4/bin/ccp4-python", "-c", conv_py, os.linesep],
                    tab_cmd + ["<< eof >", os.devnull],
                    [tab_stdin],
                    ["eof"], [os.linesep],
                    rot_cmd + ["<< eof >", rot_log],
                    [rot_stdin],
                    ["eof"], [os.linesep],
                    ["grep", "-m 1", "SOLUTIONRCD", rot_log, os.linesep],
                    ["rm", "-rf", "$CCP4_SCR\n"],
                    [EXPORT, "CCP4_SCR=" + ccp4_scr],
                ]
                amore_script = pyjob.misc.make_script(
                    cmd, directory=script_log_dir, prefix="amore_", stem=dat_model.pdb_code
                )
                amore_log = amore_script.rsplit(".", 1)[0] + '.log'
                amore_files += [(amore_script, tab_stdin, rot_stdin,
                                 amore_log, dat_model.dat_path)]

            results = []
            if len(amore_files) > 0:
                logger.info("Running AMORE tab/rot functions")
                amore_scripts, _, _, amore_logs, dat_models = zip(*amore_files)
                simbad.rotsearch.submit_chunk(amore_scripts, script_log_dir, nproc, 'simbad_amore',
                                              submit_qtype, submit_queue, monitor, self.rot_succeeded_log)

                for dat_model, amore_log in zip(dat_models, amore_logs):
                    base = os.path.basename(amore_log)
                    pdb_code = base.replace("amore_", "").replace(".log", "")
                    try:
                        rotsearch_parser = simbad.parsers.rotsearch_parser.AmoreRotsearchParser(
                            amore_log
                        )
                        score = simbad.core.amore_score.AmoreRotationScore(
                            pdb_code, dat_model, rotsearch_parser.alpha, rotsearch_parser.beta, rotsearch_parser.gamma,
                            rotsearch_parser.cc_f, rotsearch_parser.rf_f, rotsearch_parser.cc_i, rotsearch_parser.cc_p,
                            rotsearch_parser.icp, rotsearch_parser.cc_f_z_score, rotsearch_parser.cc_p_z_score,
                            rotsearch_parser.num_of_rot
                        )
                        if rotsearch_parser.cc_f_z_score:
                            results += [score]
                    except IOError:
                        pass

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
            "ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
            "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak"
        ]
        summarize_result(self.search_results,
                         csv_file=csv_file, columns=columns)

    def _generate_hklpck0(self):
        mtz_labels = simbad.util.mtz_util.GetLabels(self.mtz)
        logger.info("Preparing files for AMORE rotation function")
        stdin = self.sortfun_stdin_template.format(f=mtz_labels.f, sigf=mtz_labels.sigf)
        hklpck0 = os.path.join(self.work_dir, 'spmipch.hkl')
        cmd = [self.amore_exe, 'hklin', self.mtz, 'hklpck0', hklpck0]
        pyjob.cexec(cmd, stdin=stdin)
        return hklpck0

    def _write_stdin(self, directory, prefix, name, content):
        f = pyjob.misc.tmp_file(directory=directory, prefix=prefix,
                                stem=name, suffix=".stdin")
        with open(f, "w") as f_out:
            f_out.write(content)
        return f

    @property
    def search_results(self):
        return sorted(self._search_results, key=lambda x: float(x.CC_F_Z_score), reverse=True)[:self.max_to_keep]

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
        if self.skip_mr:
            return False

        rot_prog, pdb = os.path.basename(log).replace('.log', '').split('_', 1)
        rotsearch_parser = simbad.parsers.rotsearch_parser.AmoreRotsearchParser(log)
        dat_model = [s for s in self.simbad_dat_files if pdb in s][0]
        score = simbad.core.amore_score.AmoreRotationScore(
            pdb, dat_model, rotsearch_parser.alpha, rotsearch_parser.beta, rotsearch_parser.gamma,
            rotsearch_parser.cc_f, rotsearch_parser.rf_f, rotsearch_parser.cc_i, rotsearch_parser.cc_p,
            rotsearch_parser.icp, rotsearch_parser.cc_f_z_score, rotsearch_parser.cc_p_z_score,
            rotsearch_parser.num_of_rot
        )
        results = [score]
        if self._rot_job_succeeded(rotsearch_parser.cc_f_z_score) and pdb not in self.tested:
            self.tested.append(pdb)
            output_dir = os.path.join(self.work_dir, "mr_search")
            mr = simbad.mr.MrSubmit(mtz=self.mtz,
                                    mr_program=self.mr_program,
                                    refine_program='refmac5',
                                    refine_type=None,
                                    refine_cycles=0,
                                    output_dir=output_dir,
                                    sgalternative='none',
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