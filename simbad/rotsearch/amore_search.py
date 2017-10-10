"""Module to run the AMORE rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "10 Oct 2017"
__version__ = "0.3"

import logging
import math
import os
import shutil
import uuid
logger = logging.getLogger(__name__)

import pyjob
import pyjob.misc

import simbad.db
import simbad.rotsearch.amore_score
import simbad.parsers.rotsearch_parser
import simbad.util.pdb_util
import simbad.util.mtz_util
import simbad.util.matthews_coef

EXPORT = "SET" if os.name == "nt" else "export"


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
    >>> rotation_search = AmoreRotationSearch('<amore_exe>', '<mtz>', '<work_dir>', '<max_to_keep>')
    >>> rotation_search.run_pdb(
    ...     '<models_dir>', '<output_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>',
    ...     '<min_solvent_content>', '<submit_qtype>', '<submit_queue>', '<monitor>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results


    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from amore.

    """

    def __init__(self, amore_exe, mtz, work_dir, max_to_keep=20):
        self.amore_exe = amore_exe
        self.max_to_keep = max_to_keep
        self.mtz = mtz
        self.work_dir = work_dir

        self._search_results = None

    def run_pdb(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, models_dir, nproc=2, shres=3.0, pklim=0.5, npic=50,
            rotastep=1.0, min_solvent_content=20, submit_qtype=None,
            submit_queue=None, monitor=None, chunk_size=0):
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
        simbad_dat_files = simbad.db.find_simbad_dat_files(models_dir)

        space_group, _, cell = simbad.util.mtz_util.crystal_data(self.mtz)
        cell = " ".join(map(str, cell))

        chunk_size = AmoreRotationSearch.get_chunk_size(len(simbad_dat_files),
                                                        chunk_size)
        total_chunk_cycles = AmoreRotationSearch.get_total_chunk_cycles(len(simbad_dat_files),
                                                                        chunk_size)

        tmp_dir = os.path.join(os.environ["CCP4_SCR"],
                               "tmp-" + str(uuid.uuid4()))
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)

        hklpck0 = self._generate_hklpck0()
	template_hklpck1 = os.path.join(tmp_dir, "{}.hkl")
	template_clmn0 = os.path.join(tmp_dir, "{}_spmipch.clmn")
	template_clmn1 = os.path.join(tmp_dir, "{}.clmn")
	template_mapout = os.path.join(tmp_dir, "{}_amore_cross.map")
	template_table1 = os.path.join(tmp_dir, "{}_sfs.tab")
	template_model = os.path.join(tmp_dir, "{}.pdb")
        amore_temp_files = os.path.join(tmp_dir,
                                        os.path.basename(self.amore_exe) + "_*${PID1}")

        sol_calc = simbad.util.matthews_coef.SolventContent(cell,
                                                            space_group)

        iteration_range = range(0, len(simbad_dat_files), chunk_size)
        for cycle, i in enumerate(iteration_range):
            logger.info("Working on chunk %d out of %d", cycle + 1,
                        total_chunk_cycles)

            amore_files = []
            for dat_model in simbad_dat_files[i:i + chunk_size]:
                name = os.path.basename(dat_model.replace(".dat", ""))
                pdb_struct = simbad.util.pdb_util.PdbStructure(dat_model)
                solvent_content = sol_calc.calculate_from_struct(pdb_struct)
                if solvent_content < min_solvent_content:
                    msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                    logger.debug(msg, name, min_solvent_content)
                    continue
                x, y, z, intrad = pdb_struct.integration_box

                logger.debug("Generating script to perform AMORE rotation "
                             + "function on %s", name)

                pdb_model = template_model.format(name)
                table1 = template_table1.format(name)
                hklpck1 = template_hklpck1.format(name)
                clmn0 = template_clmn0.format(name)
                clmn1 = template_clmn1.format(name)
                mapout = template_mapout.format(name)

                conv_cmd = ["ccp4-python", "-c",
                            "\"import simbad.db; simbad.db.convert_dat_to_pdb('{}', '{}')\""]
                conv_cmd[-1] = conv_cmd[-1].format(dat_model, pdb_model)

                rot_stdin = self._write_stdin(
                    tmp_dir, "rotfun_", name,
                    self.rotfun_stdin_template.format(
                        shres=shres, intrad=intrad, pklim=pklim, npic=npic,
                        step=rotastep
                    )
                )
                rot_cmd = [self.amore_exe, 'table1', table1, 'HKLPCK1', hklpck1,
                           'hklpck0', hklpck0, 'clmn1', clmn1, 'clmn0', clmn0,
                           'MAPOUT', mapout]

                tab_stdin = self._write_stdin(
                    tmp_dir, "tabfun_", name,
                    self.tabfun_stdin_template.format(
                        x=x, y=y, z=z, a=90, b=90, c=120
                    )
                )
                tab_cmd = [self.amore_exe, "xyzin1", pdb_model, "xyzout1",
                           pdb_model, "table1", table1]

                removables = [clmn0, clmn1, hklpck1, table1, mapout, pdb_model,
                              amore_temp_files]
                amore_script = pyjob.misc.make_script(
                    [
                        [EXPORT, "CCP4_SCR=" + tmp_dir, os.linesep],
                        conv_cmd, os.linesep,
                        tab_cmd + ["<", tab_stdin, "&", os.linesep],
                        ["export PID1=$!", "&&", "wait", os.linesep],
                        rot_cmd + ["<", rot_stdin, os.linesep],
                        ["rm"] + removables
                    ],
                    directory=tmp_dir, prefix="amore_", stem=name
                )
                amore_log = amore_script.rsplit(".", 1)[0] + '.log'
                amore_files += [(amore_script, tab_stdin, rot_stdin,
                                 amore_log, dat_model)]

            results = []
            if len(amore_files) > 0:
                logger.info("Running AMORE tab/rot functions")
                amore_scripts, _, _, amore_logs, dat_models = zip(*amore_files)
                self.submit_chunk(amore_scripts, tmp_dir, nproc, 'simbad_amore',
                                  submit_qtype, submit_queue, monitor)

                for dat_model, amore_log in zip(dat_models, amore_logs):
                    base = os.path.basename(amore_log)
                    pdb_code = base.replace("amore_", "").replace(".log", "")
                    RP = simbad.parsers.rotsearch_parser.RotsearchParser(
                        amore_log
                    )
                    score = simbad.rotsearch.amore_score.AmoreRotationScore(
                        pdb_code, dat_model, RP.alpha, RP.beta, RP.gamma,
                        RP.cc_f, RP.rf_f, RP.cc_i, RP.cc_p, RP.icp,
                        RP.cc_f_z_score, RP.cc_p_z_score, RP.num_of_rot
                    )
                    if RP.cc_f_z_score:
                        results += [score]

            else:
                logger.critical("No structures to be trialled")

        self._search_results = results
        shutil.rmtree(tmp_dir)

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
        logger.info("Preparing files for AMORE rotation function")
        f, sigf, _, _, _, _, _ = simbad.util.mtz_util.get_labels(self.mtz)
        stdin = self.sortfun_stdin_template.format(f=f, sigf=sigf)
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
        return os.linesep.join([
            "TITLE   ** spmi  packing h k l F for crystal**",
            "SORTFUN RESOL 100.  2.5",
            "LABI FP={f}  SIGFP={sigf}",
        ])

    @property
    def tabfun_stdin_template(self):
        return os.linesep.join([
            "TITLE: Produce table for MODEL FRAGMENT",
            "TABFUN",
            "CRYSTAL {x} {y} {z} {a} {b} {c} ORTH 1",
            "MODEL 1 BTARGET 23.5",
            "SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0",
        ])

    @property
    def rotfun_stdin_template(self):
        return os.linesep.join([
            "TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1",
            "ROTFUN",
            "GENE 1   RESO 100.0 {shres}  CELL_MODEL 80 75 65",
            "CLMN CRYSTAL ORTH  1 RESO  20.0  {shres}  SPHERE   {intrad}",
            "CLMN MODEL 1     RESO  20.0  {shres} SPHERE   {intrad}",
            "ROTA  CROSS  MODEL 1  PKLIM {pklim}  NPIC {npic} STEP {step}"
        ])

    @staticmethod
    def get_chunk_size(total, size):
        return total if size == 0 else size

    @staticmethod
    def get_total_chunk_cycles(total, step):
        total_chunk_cycles, remainder = divmod(total, step)
        if remainder > 0:
            return total_chunk_cycles + 1
        else:
            return total_chunk_cycles

    @staticmethod
    def submit_chunk(chunk_scripts, run_dir, nproc, job_name, submit_qtype, submit_queue, monitor):
        """Submit jobs in small chunks to avoid using too much disk space
        
        Parameters
        ----------
        chunk_scripts : list
            List of scripts for each chunk
        nproc : int, optional
            The number of processors to run the job on
        job_name : str
            The name of the job to submit
        submit_qtype : str
            The cluster submission queue type - currently support SGE and LSF
        submit_queue : str
            The queue to submit to on the cluster

        """
        j = pyjob.Job(submit_qtype)
        j.submit(chunk_scripts, directory=run_dir, name=job_name,
                 nproc=nproc, queue=submit_queue, permit_nonzero=True)
        interval = int(math.log(len(chunk_scripts)) / 3)
        interval_in_seconds = interval if interval >= 5 else 5
        j.wait(interval=interval_in_seconds, monitor=monitor)
