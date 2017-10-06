"""Module to interact with pyrvapi"""

__author__ = "Adam Simpkin"
__date__ = "04 May 2017"
__version__ = "0.1"

from collections import OrderedDict

import json
import logging
import os
import pandas
import pyrvapi
import subprocess
import uuid
import urlparse

logger = logging.getLogger(__name__)


class RvapiMetadata(object):
    def __init__(self):
        self.first_tab_id = None
        self.xyz = OrderedDict()
        self.map = OrderedDict()
        self.dmap = OrderedDict()
        self.mtz = OrderedDict()
        self.ref_log = OrderedDict()
        self.mr_log = OrderedDict()

    @property
    def n_entries(self):
        assert len(self.xyz) == len(self.map) \
            and len(self.xyz) == len(self.dmap) \
            and len(self.xyz) == len(self.mtz) \
            and len(self.xyz) == len(self.ref_log) \
            and len(self.xyz) == len(self.mr_log)
        return len(self.xyz)

    def add_xyz(self, k, v):
        assert k not in self.xyz
        self.xyz[k] = v

    def add_mr_log(self, k, v):
        assert k not in self.mr_log
        self.mr_log[k] = v

    def add_ref_log(self, k, v):
        assert k not in self.ref_log
        self.ref_log[k] = v

    def add_map(self, k, v):
        assert k not in self.map
        self.map[k] = v

    def add_dmap(self, k, v):
        assert k not in self.dmap
        self.dmap[k] = v

    def add_mtz(self, k, v):
        assert k not in self.mtz
        self.mtz[k] = v

    def to_json(self):
        self.__dict__.update({"nEntries": self.n_entries})
        return json.dumps(self.__dict__)


class SimbadOutput(object):
    """Class to display the output of SIMBAD

    Attributes
    ----------
    webserver_uri : str
        The uri if run on a webserver
    display_gui : bool
        Option to prevent results being displayed
    logfile : str
        Path to the log file
    work_dir : str
        Path to the work directory [default: None]
    summary : bool
        Option to display summary tab [default: False]

    Examples
    --------
    >>> from simbad.util import pyrvapi_results
    >>> gui = pyrvapi_results.SimbadOutput()
    >>> gui.display_results(<'rvapi_file'>, <'webserver_uri'>, <'display__gui'>, 
    ...                     <'logfile'>, <'work_dir'>, <'summary'>)
    """
    _simbad_tooltips = {"PDB_code": "The 4 letter code representing the protein in the protein data bank",
                        "alt": "Alternate Niggli Cell",
                        "a": "Lattice parameter a",
                        "b": "Lattice parameter b",
                        "c": "Lattice parameter c",
                        "alpha": "Lattice parameter alpha",
                        "beta": "Lattice parameter beta",
                        "gamma": "Lattice parameter gamma",
                        "length_penalty": "The sum of the differences between lattice parameters a, b and c for the "
                        "model and the target",
                        "angle_penalty": "The sum of the differences between lattice parameters alpha, beta and gamma "
                        "for the model and the target",
                        "total_penalty": "The sum of the length penalty and the angle penalty",
                        "volume_difference": "The difference in volume between the query and reference unit cells",
                        "probability_score": "The probability that the structure corresponding to the total lattice "
                        "penalty will result in a solution",
                        "molrep_score": "MOLREP score for the Molecular Replacement solution",
                        "molrep_tfscore": "MOLREP translation function score for the Molecular Replacement solution",
                        "phaser_llg": "PHASER Log-likelihood gain for the Molecular Replacement solution",
                        "phaser_tfz": "PHASER Translation Function Z-score for the Molecular Replacement solution",
                        "phaser_rfz": "PHASER Rotational Function Z-score for the Molecular Replacement solution",
                        "final_r_fact": "R-fact score for REFMAC refinement of the Molecular Replacement solution",
                        "final_r_free": "R-free score for REFMAC refinement of the Molecular Replacement solution",
                        "peaks_over_6_rms": "Anomalous peaks over 6 RMS",
                        "peaks_over_6_rms_within_4a_of_model": "Anomalous peaks over 6 RMS within 4 Angstroms of the "
                        "Molecular Replacement solution",
                        "peaks_over_9_rms": "Anomalous peaks over 9 RMS",
                        "peaks_over_9_rms_within_4a_of_model": "Anomalous peaks over 9 RMS within 4 Angstroms of the "
                        "Molecular Replacement solution",
                        "ALPHA": "Lattice parameter alpha",
                        "BETA": "Lattice parameter beta",
                        "GAMMA": "Lattice parameter gamma",
                        "CC_F": "The correlation coefficient between the observed amplitudes for the crystal and the "
                        "calculated amplitudes for the model",
                        "RF_F": "The classic R factor between the observed amplitudes for the crystal and the "
                        "calculated amplitudes for the model",
                        "CC_I": "The correlation coefficient between the observed intensities for the crystal and the "
                        "sum of calculated intensities for all symmetry equivalents of the model",
                        "CC_P": "The Patterson correlation coefficient between the crystal and the model Pattersons "
                        "evaluated within the defined sphere centred on the Patterson origin",
                        "Icp": "",
                        "CC_F_Z_score": "Z-score of CC_F peaks",
                        "CC_P_Z_score": "Z-score of CC_P peaks",
                        "Number_of_rotation_searches_producing_peak": "Number of rotations searches which produce "
                        "each peak [out of 5]"}

    def __init__(self, work_dir):
        self.work_dir = work_dir
        self.jsrview_dir = os.path.join(work_dir, "jsrview")

        self.running = None
        self.webserver_uri = None
        self._webserver_start = None
        self.log_tab_id = None
        self.lattice_results_tab_id = None
        self.lattice_df = None
        self.contaminant_results_tab_id = None
        self.contaminant_df = None
        self.morda_db_results_tab_id = None
        self.morda_db_df = None
        self.summary_tab_id = None
        self.summary_tab_results_sec_id = None

        self.jscofe_mode = False
        self.rhs_tab_id = None
        self.rvapi_meta = RvapiMetadata()

    def _add_tab_to_pyrvapi(self, id, title, opened):
        if self.jscofe_mode:
            pyrvapi.rvapi_insert_tab(id, title, self.rhs_tab_id, opened)
        else:
            pyrvapi.rvapi_add_tab(id, title, opened)

    def create_log_tab(self, logfile):
        """Function to create log tab

        Parameters
        ----------
        logfile : str
            Path to the log file

        Returns
        -------
        str
            self.log_tab_id
        object
            Updating page containing log
        """
        if self.jscofe_mode or self.log_tab_id:
            return
        if not os.path.isfile(logfile):
            return False

        self.log_tab_id = "log_tab"
        logurl = self.fix_path(logfile)
        self._add_tab_to_pyrvapi(self.log_tab_id, "Log file", True)
        pyrvapi.rvapi_append_content(logurl, True, self.log_tab_id)
        return self.log_tab_id

    def _create_lattice_results_tab(self):
        """Function to create lattice results tab

        Returns
        -------
        str
            self.lattice_results_tab_id
        object
            Empty page to append lattice results to
        """
        if not self.lattice_results_tab_id:
            self.lattice_results_tab_id = "lattice_results_tab"
            self._add_tab_to_pyrvapi(self.lattice_results_tab_id,
                                     "Lattice Parameter Search Results", False)

    def _create_contaminant_results_tab(self):
        """Function to create contaminant results tab

        Returns
        -------
        str
            self.contaminant_results_tab_id
        object
            Empty page to append contaminant results to
        """
        if not self.contaminant_results_tab_id:
            self.contaminant_results_tab_id = "contaminants_results_tab"
            self._add_tab_to_pyrvapi(self.contaminant_results_tab_id,
                                     "Contaminant Search Results", False)

    def _create_morda_db_results_tab(self):
        """Function to create morda db results tab

        Returns
        -------
        str
            self.morda_db_results_tab_id
        object
            Empty page to append morda results to
        """
        if not self.morda_db_results_tab_id:
            self.morda_db_results_tab_id = "morda_db_results_tab"
            self._add_tab_to_pyrvapi(self.morda_db_results_tab_id,
                                     "MoRDa Database Search Results", False)

    def _create_summary_tab(self):
        """Function to create a summary tab

        Returns
        -------
        str
            self.summary_tab_id
        object
            Empty page to append summary to
        """
        if not self.summary_tab_id:
            self.summary_tab_id = "summary_tab"
            self._add_tab_to_pyrvapi(self.summary_tab_id, "Summary", True)

    def create_lattice_results_tab(self, lattice_results, lattice_mr_results):
        """Function to create the lattice results tab

        Parameters
        ----------
        lattice_results : str
            Path to the file containing the lattice results
        lattice_mr_results : str
            Path to the file containing the lattice MR results

        Returns
        -------
        object
            Page containing the results from the lattice parameter search
        """

        self._create_lattice_results_tab()

        if os.path.isfile(lattice_results):
            section_title = 'Lattice Parameter Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.lattice_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, True)

            table_title = "Lattice Parameter Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(lattice_results)
            self.create_table(df, table)

        if os.path.isfile(lattice_mr_results):
            section_title = 'Molecular Replacement Search Results'

            uid = str(uuid.uuid4())
            sec = section_title.replace(" ", "_") + uid
            tab = self.lattice_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, True)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(lattice_mr_results)
            self.create_table(df, table)

            section_title = 'Top 10 Lattice Parameter Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top 10 Lattice Parameter Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            self.lattice_df = df

            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(
                        self.work_dir, 'latt', 'mr_lattice', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(
                        mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(
                        mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

                    prefix = "{}_latt_".format(i)
                    self.output_result_files(download_sec, diff_map, ref_map,
                                             ref_mtz, ref_pdb, prefix=prefix)
                    self.output_log_files(logfile_sec, mr_log, ref_log,
                                          prefix=prefix)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def create_contaminant_results_tab(self, contaminant_results, contaminant_mr_results):
        """Function to create the contaminant results tab

        Parameters
        ----------
        contaminant_results : str
            Path to the file containing the contaminant results
        contaminant_mr_results : str
            Path to the file containing the contaminant MR results

        Returns
        -------
        object
            Page containing the results from the contaminant search
        """
        self._create_contaminant_results_tab()

        if os.path.isfile(contaminant_results):
            section_title = 'Contaminant AMORE Rotation Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.contaminant_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Contaminant AMORE Rotation Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(contaminant_results)
            self.create_table(df, table)

        if os.path.isfile(contaminant_mr_results):
            section_title = 'Molecular Replacement Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.contaminant_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(contaminant_mr_results)
            self.create_table(df, table)

            self.contaminant_df = df

            section_title = "Molecular Replacement Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(
                graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)

            section_title = 'Top 10 Contaminant Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top 10 Contaminant Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(
                        self.work_dir, 'cont', 'mr_contaminant', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(
                        mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(
                        mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

                    prefix = "{}_cont_".format(i)
                    self.output_result_files(download_sec, diff_map, ref_map,
                                             ref_mtz, ref_pdb, prefix=prefix)
                    self.output_log_files(logfile_sec, mr_log, ref_log,
                                          prefix=prefix)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def create_morda_db_results_tab(self, morda_db_results, morda_db_mr_results):
        """Function to create the MoRDa Database results tab

        Parameters
        ----------
        morda_db_results : str
            Path to the file containing the MoRDa db results
        morda_db_mr_results : str
            Path to the file containing the MoRDa db MR results

        Returns
        -------
        object
            Page containing the results from the MoRDa db search
        """
        self._create_morda_db_results_tab()

        if os.path.isfile(morda_db_results):
            section_title = 'MoRDa database AMORE Rotation Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.morda_db_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "MoRDa datbase AMORE Rotation Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(morda_db_results)
            self.create_table(df, table)

        if os.path.isfile(morda_db_mr_results):
            section_title = 'Molecular Replacement Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.morda_db_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(
                sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(morda_db_mr_results)
            self.create_table(df, table)

            self.morda_db_df = df

            section_title = "Molecular Replacement Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(
                graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)

            section_title = 'Top 10 MoRDa database Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top 10 MoRDa database Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(
                logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(
                        self.work_dir, 'morda', 'mr_morda', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(
                        mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(
                        mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(
                        mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(
                        mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

                    prefix = "{}_morda_".format(i)
                    self.output_result_files(download_sec, diff_map, ref_map,
                                             ref_mtz, ref_pdb, prefix=prefix)
                    self.output_log_files(logfile_sec, mr_log, ref_log,
                                          prefix=prefix)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def create_summary_tab(self):
        """Function to create the MoRDa Database results tab

        Returns
        -------
        object
            Page containing a summary of the best results from SIMBAD
        """
        self._create_summary_tab()

        if self.lattice_df is None:
            lattice_score = 1
        else:
            try:
                lattice_score = self.lattice_df['final_r_free'][0]
            except IndexError:
                lattice_score = 1

        if self.contaminant_df is None:
            contaminant_score = 1
        else:
            try:
                contaminant_score = self.contaminant_df['final_r_free'][0]
            except IndexError:
                contaminant_score = 1

        if self.morda_db_df is None:
            morda_db_score = 1
        else:
            try:
                morda_db_score = self.morda_db_df['final_r_free'][0]
            except IndexError:
                morda_db_score = 1

        if lattice_score == 1 and contaminant_score == 1 and morda_db_score == 1:
            section_title = 'SIMBAD Summary'
            uid = str(uuid.uuid4())
            sec = section_title.replace(" ", "_") + uid
            tab = self.summary_tab_id

            msg = "No solution was found by SIMBAD"

            pyrvapi.rvapi_add_section(
                sec, section_title, tab, 0, 0, 1, 1, True)
            pyrvapi.rvapi_add_text(msg, sec, 2, 0, 1, 1)
            return

        elif lattice_score <= contaminant_score and lattice_score <= morda_db_score:
            pdb_code = self.lattice_df.loc[0][0]
            r_fact = self.lattice_df['final_r_fact'][0]
            r_free = self.lattice_df['final_r_free'][0]
            mr_program = list(self.lattice_df)[1][0:6]
            mr_workdir = os.path.join(
                self.work_dir, 'latt', 'mr_lattice', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine',
                                   '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine',
                                   '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(
                mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

        elif contaminant_score <= lattice_score and contaminant_score <= morda_db_score:
            pdb_code = self.contaminant_df.loc[0][0]
            r_fact = self.contaminant_df['final_r_fact'][0]
            r_free = self.contaminant_df['final_r_free'][0]
            mr_program = list(self.contaminant_df)[1][0:6]
            mr_workdir = os.path.join(
                self.work_dir, 'cont', 'mr_contaminant', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine',
                                   '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine',
                                   '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(
                mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

        elif morda_db_score <= lattice_score and morda_db_score <= contaminant_score:
            pdb_code = self.morda_db_df.loc[0][0]
            r_fact = self.morda_db_df['final_r_fact'][0]
            r_free = self.morda_db_df['final_r_free'][0]
            mr_program = list(self.morda_db_df)[1][0:6]
            mr_workdir = os.path.join(
                self.work_dir, 'morda', 'mr_morda', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(
                mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine',
                                   '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine',
                                   '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(
                mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))

        else:
            logger.debug('Unexpected result')
            return

        section_title = 'SIMBAD Summary'
        uid = str(uuid.uuid4())
        sec = section_title.replace(" ", "_") + uid
        tab = self.summary_tab_id

        msg = 'The best search model found by SIMBAD was {0}. \
               This gave an R/Rfact of {1:.3f} and an R/Rfree of {2:.3f}. \
               An R/Rfree lower than 0.450 is indicative of a \
               solution. Values above this may also be indicative of a correct solution \
               but you should examine the maps through the graphical map viewer for \
               verification'.format(pdb_code, r_fact, r_free)

        pyrvapi.rvapi_add_section(
            sec, section_title, tab, 0, 0, 1, 1, True)
        pyrvapi.rvapi_add_text(msg, sec, 2, 0, 1, 1)

        section_title = 'Best SIMBAD result Downloads'
        uid = str(uuid.uuid4())
        download_sec = section_title.replace(" ", "_") + uid
        pyrvapi.rvapi_add_section(
            download_sec, section_title, tab, 0, 0, 1, 1, True)

        section_title = 'Best SIMBAD result Log Files'
        uid = str(uuid.uuid4())
        logfile_sec = section_title.replace(" ", "_") + uid
        pyrvapi.rvapi_add_section(
            logfile_sec, section_title, tab, 0, 0, 1, 1, False)

        prefix = "best_"
        self.output_result_files(download_sec, diff_map, ref_map,
                                 ref_mtz, ref_pdb, prefix=prefix)
        self.output_log_files(logfile_sec, mr_log, ref_log,
                              prefix=prefix)

    def output_result_files(self, sec, diff_map, ref_map, ref_mtz, ref_pdb, prefix=""):
        """Function to display the result files for the result

        Parameters
        ----------
        sec : str
            Section the output results files will be added to
        diff_map : str
            Path to the difference map
        ref_map : str
            Path to the refined map
        ref_mtz : str
            Path to the refined mtz
        ref_pdb : str
            Path to the refined pdb
        prefix : str
            Prefix for rvapi meta key

        Returns
        -------
        object
            Section containing the pdb and mtz for a result
        """
        title = "Electron density for {0}".format(
            os.path.basename(ref_pdb).split('_')[0])
        uid = str(uuid.uuid4())
        data = "dat" + uid

        pyrvapi.rvapi_add_data1(os.path.join(sec, data), title,
                                ref_pdb, "xyz", 2, 0, 1, 1, 1)
        self.rvapi_meta.add_xyz(prefix + os.path.basename(ref_pdb).replace("_refinement_output.pdb", ""),
                                self.rel_path(ref_pdb))

        pyrvapi.rvapi_append_to_data(data, ref_mtz, "hkl:map")
        self.rvapi_meta.add_mtz(prefix + os.path.basename(ref_mtz).replace("_refinement_output.mtz", ""),
                                self.rel_path(ref_mtz))

        pyrvapi.rvapi_append_to_data(data, ref_map, "hkl:ccp4_map")
        self.rvapi_meta.add_map(prefix + os.path.basename(ref_map).replace("_refmac_2fofcwt.map", ""),
                                self.rel_path(ref_map))

        pyrvapi.rvapi_append_to_data(data, diff_map, "hkl:ccp4_dmap")
        self.rvapi_meta.add_dmap(prefix + os.path.basename(diff_map).replace("_refmac_fofcwt.map", ""),
                                 self.rel_path(diff_map))

    def output_log_files(self, sec, mr_log, ref_log, prefix=""):
        """Function to display the log files for the result

        Parameters
        ----------
        sec : str
            Section the output logs will be added to
        mr_log : str
            Path to the output MR log
        ref_log : str
            Path to the output refinement log
        prefix : str
            Prefix for rvapi meta key

        Returns
        -------
        object
            Section containing mr and refinement logs
        """
        title = "Log files from {0}".format(
            os.path.basename(mr_log).split('_')[0])

        uid = str(uuid.uuid4())
        data = "dat" + uid

        pyrvapi.rvapi_add_data1(os.path.join(sec, data), title,
                                mr_log, "text", 2, 0, 1, 1, 0)
        self.rvapi_meta.add_mr_log(prefix + os.path.basename(mr_log).replace("_mr.log", ""),
                                   self.rel_path(mr_log))

        uid = str(uuid.uuid4())
        data = "dat" + uid

        pyrvapi.rvapi_add_data1(os.path.join(sec, data), "",
                                ref_log, "text", 2, 0, 1, 1, 0)
        self.rvapi_meta.add_ref_log(prefix + os.path.basename(ref_log).replace("_ref.log", ""),
                                    self.rel_path(ref_log))

    def create_table(self, df, table_id):
        """Function to create/display tables

        Parameters
        ----------
        df : :Pandas dataframe:
            Input Pandas dataframe contaiing data to be plotted
        table_id : str
            Table ID

        Returns
        -------
        object
            table containing the results from SIMBAD
        """
        num_labels = 0
        for i, l in enumerate(df):
            if i == 0:
                pyrvapi.rvapi_put_horz_theader(
                    table_id, "PDB_code", self._simbad_tooltips["PDB_code"], 0)
            else:
                pyrvapi.rvapi_put_horz_theader(
                    table_id, l, self._simbad_tooltips[l], i)
                num_labels = i

        ir = len(df)
        for i in range(0, ir):
            for j in range(num_labels + 1):
                if j == 0:
                    pyrvapi.rvapi_put_table_string(table_id,
                                                   '<a href="http://www.ebi.ac.uk/pdbe/entry/pdb/{0}" '
                                                   'target="_blank">{1}</a>'.format(df.loc[i][j][0:4],
                                                                                    df.loc[i][j]), i, j)
                else:
                    pyrvapi.rvapi_put_table_string(
                        table_id, str(df.loc[i][j]), i, j)

    @staticmethod
    def create_graphs(df, graph_sec, graph_widget):
        """Function to create/display graphs following MR

        df : :Pandas dataframe:
            Input Pandas dataframe containing data to be plotted
        graph_sec : str
            Section the output graph will be displayed in
        graph_widget : str
            Widget ID

        Returns
        -------
        object
            Section containing the graphic representation of results from SIMBAD
        """
        pyrvapi.rvapi_append_loggraph1(
            os.path.join(graph_sec, graph_widget))

        pyrvapi.rvapi_add_graph_data1(
            graph_widget + "/data1", " Scores Vs. Rank (by R-Free)")
        pyrvapi.rvapi_add_graph_dataset1(
            graph_widget + "/data1/x", "Rank", " (by R-Free)")

        pyrvapi.rvapi_add_graph_dataset1(
            graph_widget + "/data1/y1", "REFMAC R-Fact", "")
        pyrvapi.rvapi_add_graph_dataset1(
            graph_widget + "/data1/y2", "REFMAC R-Free", "")

        mr_program = list(df)[1][0:6]
        if mr_program == 'molrep':
            pyrvapi.rvapi_add_graph_dataset1(
                graph_widget + "/data1/y3", "MOLREP score", "")
            pyrvapi.rvapi_add_graph_dataset1(
                graph_widget + "/data1/y4", "MOLREP TF/sig", "")

        elif mr_program == 'phaser':
            pyrvapi.rvapi_add_graph_dataset1(
                graph_widget + "/data1/y3", "PHASER TFZ", "")
            pyrvapi.rvapi_add_graph_dataset1(
                graph_widget + "/data1/y4", "PHASER LLG", "")
            pyrvapi.rvapi_add_graph_dataset1(
                graph_widget + "/data1/y5", "PHASER RFZ", "")

        ir = len(df.index)
        for i in range(0, ir):
            if df['final_r_free'][i] < 0.7 and df['final_r_fact'][i] < 0.7:
                pyrvapi.rvapi_add_graph_int1(graph_widget + "/data1/x", i + 1)
                pyrvapi.rvapi_add_graph_real1(
                    graph_widget + "/data1/y1", df['final_r_fact'][i], "%g")
                pyrvapi.rvapi_add_graph_real1(
                    graph_widget + "/data1/y2", df['final_r_free'][i], "%g")
                if mr_program == 'molrep':
                    pyrvapi.rvapi_add_graph_real1(
                        graph_widget + "/data1/y3", df['molrep_score'][i], "%g")
                    pyrvapi.rvapi_add_graph_real1(
                        graph_widget + "/data1/y4", df['molrep_tfscore'][i], "%g")

                elif mr_program == 'phaser':
                    pyrvapi.rvapi_add_graph_real1(
                        graph_widget + "/data1/y3", df['phaser_tfz'][i], "%g")
                    pyrvapi.rvapi_add_graph_real1(
                        graph_widget + "/data1/y4", df['phaser_llg'][i], "%g")
                    pyrvapi.rvapi_add_graph_real1(
                        graph_widget + "/data1/y5", df['phaser_rfz'][i], "%g")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot1", "R-Fact/R-Free Vs. Rank", "Rank (by R-Free)",
                                      "R-Fact/R-Free")
        pyrvapi.rvapi_add_plot_line1(
            graph_widget + "/data1/plot1", "x", "y1")
        pyrvapi.rvapi_add_plot_line1(
            graph_widget + "/data1/plot1", "x", "y2")
        pyrvapi.rvapi_set_plot_xmin("plot1", "graphWidget1", -1.0)

        if mr_program == 'molrep':
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot2", "MOLREP score Vs. Rank", "Rank (by R-Free)",
                                          "MOLREP score")
            pyrvapi.rvapi_add_plot_line1(
                graph_widget + "/data1/plot2", "x", "y3")

            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot3", "MOLREP TF/sig Vs. Rank", "Rank (by R-Free)",
                                          "MOLREP TF/sig")
            pyrvapi.rvapi_add_plot_line1(
                graph_widget + "/data1/plot3", "x", "y4")

        elif mr_program == 'phaser':
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot2", "PHASER TFZ Vs. Rank", "Rank (by R-Free)",
                                          "PHASER TFZ")
            pyrvapi.rvapi_add_plot_line1(
                graph_widget + "/data1/plot2", "x", "y3")

            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot3", "PHASER LLG Vs. Rank", "Rank (by R-Free)",
                                          "PHASER LLG")
            pyrvapi.rvapi_add_plot_line1(
                graph_widget + "/data1/plot3", "x", "y4")

            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot4", "PHASER RFZ Vs. Rank", "Rank (by R-Free)",
                                          "PHASER RFZ")
            pyrvapi.rvapi_add_plot_line1(
                graph_widget + "/data1/plot4", "x", "y5")

    def display_results(self, rvapi_document, webserver_uri, display_gui, logfile, summary=False):
        """Function to display the results

        Parameters
        ----------
        rvapi_document : str
            A pre-created rvapi file
        webserver_uri : str
            The uri if run on a webserver
        display_gui : bool
            Option to prevent results being displayed
        logfile : str
            Path to the log file
        summary : bool
            Option to display summary tab [default: False]

        Returns
        -------
        object
            GUI displaying the results of SIMBAD
        file
            index.html containing the html of the GUI
        """
        if not display_gui:
            return

        if not self.running:
            # Infrastructure to run
            ccp4 = os.environ["CCP4"]
            share_jsrview = os.path.join(ccp4, "share", "jsrview")

            os.mkdir(self.jsrview_dir)

            if rvapi_document:
                pyrvapi.rvapi_restore_document2(rvapi_document)
                self.rhs_tab_id = pyrvapi.rvapi_get_meta()
                self.jscofe_mode = True
            else:
                pyrvapi.rvapi_init_document("SIMBAD_results", self.jsrview_dir, "SIMBAD Results", 1, 7, share_jsrview, None, None,
                                            None, None)
            if webserver_uri:
                self._webserver_start = len(self.jsrview_dir) + 1
                self.webserver_uri = webserver_uri
            else:
                # We start our own browser
                jsrview = os.path.join(ccp4, "libexec", "jsrview")
                subprocess.Popen(
                    [jsrview, os.path.join(self.jsrview_dir, "index.html")])
            pyrvapi.rvapi_add_header("SIMBAD Results")
            self.running = True

            if os.path.isfile(logfile):
                self.create_log_tab(logfile)

        lattice_results = os.path.join(
            self.work_dir, 'latt', 'lattice_search.csv')
        lattice_mr_results = os.path.join(
            self.work_dir, 'latt', 'lattice_mr.csv')
        if os.path.isfile(lattice_results) or os.path.isfile(lattice_mr_results):
            self.create_lattice_results_tab(
                lattice_results, lattice_mr_results)

        contaminant_results = os.path.join(
            self.work_dir, 'cont', 'rot_search.csv')
        contaminant_mr_results = os.path.join(
            self.work_dir, 'cont', 'cont_mr.csv')
        if os.path.isfile(contaminant_results) or os.path.isfile(contaminant_mr_results):
            self.create_contaminant_results_tab(
                contaminant_results, contaminant_mr_results)

        morda_db_results = os.path.join(
            self.work_dir, 'morda', 'rot_search.csv')
        morda_db_mr_results = os.path.join(
            self.work_dir, 'morda', 'morda_mr.csv')
        if os.path.isfile(morda_db_results) or os.path.isfile(morda_db_mr_results):
            self.create_morda_db_results_tab(
                morda_db_results, morda_db_mr_results)

        if summary:
            self.create_summary_tab()

        pyrvapi.rvapi_flush()

    def save_document(self, rvapi_document):
        self.rvapi_meta.first_tab_id = self.log_tab_id
        pyrvapi.rvapi_put_meta(self.rvapi_meta.to_json())
        pyrvapi.rvapi_store_document2(rvapi_document)

    def fix_path(self, path):
        if self.webserver_uri:
            return urlparse.urljoin(self.webserver_uri, path[self._webserver_start:])
        else:
            return path

    def rel_path(self, path):
        return os.path.join(".", os.path.relpath(path, self.work_dir.split(os.sep, 1)[0]))
