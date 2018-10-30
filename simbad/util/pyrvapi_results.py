"""Module to interact with pyrvapi"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "06 Oct 2017"
__version__ = "0.2"

import json
import logging
import os
import pandas
import pyrvapi
import subprocess
import uuid
import urlparse

from simbad.util import SIMBAD_PYRVAPI_SHAREDIR

logger = logging.getLogger(__name__)


class RvapiMetadata(object):
    """Storage container for metadata required by JsCoFe"""

    def __init__(self):
        self.results = []

    @property
    def n_results(self):
        return len(self.results)

    def add(self, e):
        self.results.append(e)

    def to_json(self):
        self.__dict__.update({"nResults": self.n_results})
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
    >>> gui = pyrvapi_results.SimbadOutput(<'rvapi_file'>, <'webserver_uri'>, <'display__gui'>,
    ...                                    <'logfile'>, <'work_dir'>)
    >>> gui.display_results(<'show_summary'>)

    """
    _simbad_tooltips = {
        "PDB_code":
        "The 4 letter code representing the protein in the protein data bank",
        "alt":
        "Alternate Niggli Cell",
        "a":
        "Lattice parameter a",
        "b":
        "Lattice parameter b",
        "c":
        "Lattice parameter c",
        "alpha":
        "Lattice parameter alpha",
        "beta":
        "Lattice parameter beta",
        "gamma":
        "Lattice parameter gamma",
        "length_penalty":
        "The sum of the differences between lattice parameters a, b and c for the "
        "model and the target",
        "angle_penalty":
        "The sum of the differences between lattice parameters alpha, beta and gamma "
        "for the model and the target",
        "total_penalty":
        "The sum of the length penalty and the angle penalty",
        "volume_difference":
        "The difference in volume between the query and reference unit cells",
        "probability_score":
        "The probability that the structure corresponding to the total lattice "
        "penalty will result in a solution",
        "molrep_score":
        "MOLREP score for the Molecular Replacement solution",
        "molrep_tfscore":
        "MOLREP translation function score for the Molecular Replacement solution",
        "phaser_llg":
        "PHASER Log-likelihood gain for the Molecular Replacement solution",
        "phaser_tfz":
        "PHASER Translation Function Z-score for the Molecular Replacement solution",
        "phaser_rfz":
        "PHASER Rotational Function Z-score for the Molecular Replacement solution",
        "final_r_fact":
        "R-fact score for REFMAC refinement of the Molecular Replacement solution",
        "final_r_free":
        "R-free score for REFMAC refinement of the Molecular Replacement solution",
        "dano_peak_height":
        "The Largest Anomalous peak found by ANODE",
        "nearest_atom":
        "The atom closest to the anomalous peak",
        "z_score":
        "Z-score calculated from all the anomalous peaks",
        "ALPHA":
        "Lattice parameter alpha",
        "BETA":
        "Lattice parameter beta",
        "GAMMA":
        "Lattice parameter gamma",
        "CC_F":
        "The correlation coefficient between the observed amplitudes for the crystal and the "
        "calculated amplitudes for the model",
        "RF_F":
        "The classic R factor between the observed amplitudes for the crystal and the "
        "calculated amplitudes for the model",
        "CC_I":
        "The correlation coefficient between the observed intensities for the crystal and the "
        "sum of calculated intensities for all symmetry equivalents of the model",
        "CC_P":
        "The Patterson correlation coefficient between the crystal and the model Pattersons "
        "evaluated within the defined sphere centred on the Patterson origin",
        "Icp":
        "",
        "CC_F_Z_score":
        "Z-score of CC_F peaks",
        "CC_P_Z_score":
        "Z-score of CC_P peaks",
        "Number_of_rotation_searches_producing_peak":
        "Number of rotations searches which produce each peak [out of 5]"
    }

    def __init__(self, rvapi_document, webserver_uri, display_gui, logfile, work_dir, ccp4i2_xml=None, tab_prefix=""):
        self.rvapi_document = rvapi_document
        self.webserver_uri = webserver_uri
        self.display_gui = display_gui
        self.logfile = logfile
        self.work_dir = work_dir
        self.ccp4i2 = bool(ccp4i2_xml)
        self.tab_prefix = tab_prefix

        self.jsrview_dir = None
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

        self.lattice_search_results_displayed = False
        self.contaminant_results_displayed = False
        self.morda_results_displayed = False

        self.jscofe_mode = False
        self.rhs_tab_id = None
        self.rvapi_meta = RvapiMetadata()

        if self.display_gui or self.ccp4i2:
            ccp4 = os.environ["CCP4"]
            share_jsrview = os.path.join(ccp4, "share", "jsrview")

            if self.rvapi_document:
                pyrvapi.rvapi_restore_document2(rvapi_document)
                self.rhs_tab_id = pyrvapi.rvapi_get_meta()
                self.jscofe_mode = True
                self.jsrview_dir = os.path.dirname(rvapi_document)
            else:
                self.jsrview_dir = os.path.join(work_dir, SIMBAD_PYRVAPI_SHAREDIR)
                os.mkdir(self.jsrview_dir)
                wintitle = "SIMBAD Results"

                if ccp4i2_xml:
                    self.init_from_ccp4i2_xml(ccp4i2_xml, self.jsrview_dir, share_jsrview, wintitle)
                else:
                    pyrvapi.rvapi_init_document("SIMBAD_results", self.jsrview_dir, wintitle, 1, 7, share_jsrview, None,
                                                None, None, None)
                    self.rvapi_document = os.path.join(self.jsrview_dir, "index.html")

            if webserver_uri:
                self._webserver_start = len(self.jsrview_dir) + 1
                self.jscofe_mode = True
            elif not ccp4i2_xml:
                # We start our own browser
                jsrview = os.path.join(ccp4, "libexec", "jsrview")
                subprocess.Popen([jsrview, os.path.join(self.jsrview_dir, "index.html")])

            pyrvapi.rvapi_add_header("SIMBAD Results")

            if os.path.isfile(logfile) and not self.ccp4i2:
                self.create_log_tab(logfile)

        pyrvapi.rvapi_flush()

    def init_from_ccp4i2_xml(self, ccp4i2_xml, pyrvapi_dir, share_jsrview, wintitle):
        """This code is largely stolen from Andrew Lebedev"""

        #// Document modes
        #define RVAPI_MODE_Silent  0x00100000
        #define RVAPI_MODE_Html    0x00000001
        #define RVAPI_MODE_Xmli2   0x00000002

        mode = pyrvapi.RVAPI_MODE_Html | bool(ccp4i2_xml) * pyrvapi.RVAPI_MODE_Xmli2

        #// Document layouts
        #define RVAPI_LAYOUT_Header   0x00000001
        #define RVAPI_LAYOUT_Toolbar  0x00000002
        #define RVAPI_LAYOUT_Tabs     0x00000004
        #define RVAPI_LAYOUT_Full     0x00000007

        xml_relpath = os.path.relpath(ccp4i2_xml, pyrvapi_dir) if ccp4i2_xml else None
        docid = 'TestRun'
        layout = pyrvapi.RVAPI_LAYOUT_Full
        html = 'index.html'

        pyrvapi.rvapi_init_document(
            docid,  # const char * docId      // mandatory
            pyrvapi_dir,  # const char * outDir     // mandatory
            wintitle,  # const char * winTitle   // mandatory
            mode,  # const int    mode       // mandatory
            layout,  # const int    layout     // mandatory
            share_jsrview,  # const char * jsUri      // needed
            None,  # const char * helpFName  // may be NULL
            html,  # const char * htmlFName  // may be NULL
            None,  # const char * taskFName  // may be NULL
            xml_relpath  # const char * xmli2FName // may be NULL
        )
        return

    def _add_tab_to_pyrvapi(self, id, title, opened):
        if self.jscofe_mode:
            self._insert_tab_to_pyrvapi(id, title, self.rhs_tab_id, opened)
        else:
            pyrvapi.rvapi_add_tab(id, title, opened)

    def _insert_tab_to_pyrvapi(self, id, title, other_tab_id, opened):
        pyrvapi.rvapi_insert_tab(id, title, other_tab_id, opened)

    def create_log_tab(self, logfile):
        """Function to create log tab

        Parameters
        ----------
        logfile : str
            Path to the log file

        Returns
        -------
        str
            Updating page containing log

        """
        if self.jscofe_mode or self.log_tab_id:
            return
        if not os.path.isfile(logfile):
            return False

        self.log_tab_id = self.tab_prefix + "log_tab"
        logurl = self.fix_path(logfile)
        self._add_tab_to_pyrvapi(self.log_tab_id, "Log file", True)
        pyrvapi.rvapi_append_content(logurl, True, self.log_tab_id)
        return self.log_tab_id

    def _create_lattice_results_tab(self):
        """Function to create lattice results tab"""
        if not self.lattice_results_tab_id:
            self.lattice_results_tab_id = self.tab_prefix + "lattice_results_tab"
            self._add_tab_to_pyrvapi(self.lattice_results_tab_id, "Lattice Parameter Search Results", False)

    def _create_contaminant_results_tab(self):
        """Function to create contaminant results tab"""
        if not self.contaminant_results_tab_id:
            self.contaminant_results_tab_id = self.tab_prefix + "contaminant_results_tab"
            self._add_tab_to_pyrvapi(self.contaminant_results_tab_id, "Contaminant Search Results", False)

    def _create_morda_db_results_tab(self):
        """Function to create morda db results tab"""
        if not self.morda_db_results_tab_id:
            self.morda_db_results_tab_id = self.tab_prefix + "morda_db_results_tab"
            self._add_tab_to_pyrvapi(self.morda_db_results_tab_id, "MoRDa Database Search Results", False)

    def _create_summary_tab(self):
        """Function to create a summary tab

        Returns
        -------
        str
            self.summary_tab_id
        object
            Empty page to append summary to
        """
        if self.summary_tab_id:
            return

        self.summary_tab_id = self.tab_prefix + "summary_tab"
        title = "Summary"
        opened = True
        if self.lattice_results_tab_id:
            self._insert_tab_to_pyrvapi(self.summary_tab_id, title, self.lattice_results_tab_id, opened)
        elif self.contaminant_results_tab_id:
            self._insert_tab_to_pyrvapi(self.summary_tab_id, title, self.contaminant_results_tab_id, opened)
        elif self.morda_db_results_tab_id:
            self._insert_tab_to_pyrvapi(self.summary_tab_id, title, self.morda_db_results_tab_id, opened)
        else:
            self._add_tab_to_pyrvapi(self.summary_tab_id, title, opened)

    def create_lattice_results_tab(self, lattice_results, lattice_mr_results, results_to_display):
        """Function to create the lattice results tab

        Parameters
        ----------
        lattice_results : str
            Path to the file containing the lattice results
        lattice_mr_results : str
            Path to the file containing the lattice MR results
        results_to_display : int
            Number of results to display

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

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)

            table_title = "Lattice Parameter Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(lattice_results)
            self.create_table(df, table)

        if os.path.isfile(lattice_mr_results):
            section_title = 'Molecular Replacement Search Results'

            uid = str(uuid.uuid4())
            sec = section_title.replace(" ", "_") + uid
            tab = self.lattice_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(lattice_mr_results)
            self.create_table(df, table)

            section_title = 'Top {0} Lattice Parameter Search Downloads'.format(results_to_display)
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top {0} Lattice Parameter Search Log Files'.format(results_to_display)
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            self.lattice_df = df

            for i in range(0, results_to_display):
                try:
                    pdb_code = df.loc[i][0]
                    mr_workdir = os.path.join(self.work_dir, 'output_files')
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, '{0}_refmac_fofcwt.map'.format(pdb_code))

                    pdb, mtz, map_, dmap, mr_log, ref_log = list(self.adjust_paths_of_files(
                        [ref_pdb, ref_mtz, ref_map, diff_map, mr_log, ref_log]
                    ))

                    self.store_entry_in_rvapi_meta(
                        i + 1, "latt", pdb_code, pdb, mtz, map_, dmap, False)
                    self.output_result_files(
                        download_sec, dmap, map_, mtz, pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def create_contaminant_results_tab(self, contaminant_results, contaminant_mr_results, results_to_display):
        """Function to create the contaminant results tab

        Parameters
        ----------
        contaminant_results : str
            Path to the file containing the contaminant results
        contaminant_mr_results : str
            Path to the file containing the contaminant MR results
        results_to_display : int
            Number of results to display

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

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Contaminant AMORE Rotation Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(contaminant_results)
            self.create_table(df, table)

            section_title = "AMORE Rotation Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)

        if os.path.isfile(contaminant_mr_results):
            section_title = 'Molecular Replacement Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.contaminant_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(contaminant_mr_results)
            self.create_table(df, table)

            self.contaminant_df = df

            section_title = 'Top {0} Contaminant Search Downloads'.format(results_to_display)
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top {0} Contaminant Search Log Files'.format(results_to_display)
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            for i in range(0, results_to_display):
                try:
                    pdb_code = df.loc[i][0]
                    mr_workdir = os.path.join(self.work_dir, 'output_files')
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, '{0}_refmac_fofcwt.map'.format(pdb_code))

                    pdb, mtz, map_, dmap, mr_log, ref_log = list(self.adjust_paths_of_files(
                        [ref_pdb, ref_mtz, ref_map, diff_map, mr_log, ref_log]
                    ))

                    self.store_entry_in_rvapi_meta(
                        i + 1, "cont", pdb_code, pdb, mtz, map_, dmap, False)
                    self.output_result_files(
                        download_sec, dmap, map_, mtz, pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def create_morda_db_results_tab(self, morda_db_results, morda_db_mr_results, results_to_display):
        """Function to create the MoRDa Database results tab

        Parameters
        ----------
        morda_db_results : str
            Path to the file containing the MoRDa db results
        morda_db_mr_results : str
            Path to the file containing the MoRDa db MR results
        results_to_display : int
            Number of results to display

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

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "MoRDa datbase AMORE Rotation Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(morda_db_results)
            self.create_table(df, table)

            section_title = "AMORE Rotation Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)

        if os.path.isfile(morda_db_mr_results):
            section_title = 'Molecular Replacement Search Results'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.morda_db_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, False)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(morda_db_mr_results)
            self.create_table(df, table)

            self.morda_db_df = df

            section_title = 'Top {0} MoRDa database Search Downloads'.format(results_to_display)
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Top {0} MoRDa database Search Log Files'.format(results_to_display)
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            for i in range(0, results_to_display):
                try:
                    pdb_code = df.loc[i][0]
                    mr_workdir = os.path.join(self.work_dir, 'output_files')
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, '{0}_refmac_fofcwt.map'.format(pdb_code))

                    pdb, mtz, map_, dmap, mr_log, ref_log = list(self.adjust_paths_of_files(
                        [ref_pdb, ref_mtz, ref_map, diff_map, mr_log, ref_log]
                    ))

                    self.store_entry_in_rvapi_meta(
                        i + 1, "full", pdb_code, pdb, mtz, map_, dmap, False)
                    self.output_result_files(
                        download_sec, dmap, map_, mtz, pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)

                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))

    def display_summary_tab(self):
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

        section_title = 'SIMBAD Summary'
        uid = str(uuid.uuid4())
        sec = section_title.replace(" ", "_") + uid
        tab = self.summary_tab_id

        if lattice_score == 1 and contaminant_score == 1 and morda_db_score == 1:
            msg = "No solution was found by SIMBAD"
            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)
            pyrvapi.rvapi_add_text(msg, sec, 2, 0, 1, 1)

        else:

            if lattice_score <= contaminant_score and lattice_score <= morda_db_score:
                pdb_code = self.lattice_df.loc[0][0]
                r_fact = self.lattice_df['final_r_fact'][0]
                r_free = self.lattice_df['final_r_free'][0]
                source = "latt"
            elif contaminant_score <= lattice_score and contaminant_score <= morda_db_score:
                pdb_code = self.contaminant_df.loc[0][0]
                r_fact = self.contaminant_df['final_r_fact'][0]
                r_free = self.contaminant_df['final_r_free'][0]
                source = "cont"
            elif morda_db_score <= lattice_score and morda_db_score <= contaminant_score:
                pdb_code = self.morda_db_df.loc[0][0]
                r_fact = self.morda_db_df['final_r_fact'][0]
                r_free = self.morda_db_df['final_r_free'][0]
                source = "morda"
            else:
                logger.debug('Unexpected result')
                return

            mr_workdir = os.path.join(self.work_dir, 'output_files', pdb_code)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, '{0}_ref.log'.format(pdb_code))
            ref_pdb = os.path.join(mr_workdir, '{0}_refinement_output.pdb'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, '{0}_refmac_2fofcwt.map'.format(pdb_code))
            ref_mtz = os.path.join(mr_workdir, '{0}_refinement_output.mtz'.format(pdb_code))
            diff_map = os.path.join(mr_workdir, '{0}_refmac_fofcwt.map'.format(pdb_code))

            msg = 'The best search model found by SIMBAD was {0}. \
                   This gave an R/Rfact of {1:.3f} and an R/Rfree of {2:.3f}. \
                   An R/Rfree lower than 0.450 is indicative of a \
                   solution. Values above this may also be indicative of a correct solution \
                   but you should examine the maps through the graphical map viewer for \
                   verification'.format(pdb_code, r_fact, r_free)

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)
            pyrvapi.rvapi_add_text(msg, sec, 2, 0, 1, 1)

            section_title = 'Best SIMBAD result Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)

            section_title = 'Best SIMBAD result Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)

            pdb, mtz, map_, dmap, mr_log, ref_log = list(
                self.adjust_paths_of_files([ref_pdb, ref_mtz, ref_map, diff_map, mr_log, ref_log]))
            for e in self.rvapi_meta.results:
                if e["name"] == pdb_code and e["source"] == source:
                    e["best"] = True
            self.output_result_files(download_sec, dmap, map_, mtz, pdb)
            self.output_log_files(logfile_sec, mr_log, ref_log)

    def output_result_files(self, sec, diff_map, ref_map, ref_mtz, ref_pdb):
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

        Returns
        -------
        object
            Section containing the pdb and mtz for a result
        """
        title = "Electron density for {0}".format(os.path.basename(ref_pdb).split('_')[0])

        data = "dat" + str(uuid.uuid4())

        pyrvapi.rvapi_add_data1(os.path.join(sec, data), title, ref_pdb, "xyz", 2, 0, 1, 1, 1)
        pyrvapi.rvapi_append_to_data(data, ref_mtz, "hkl:map")
        pyrvapi.rvapi_append_to_data(data, ref_map, "hkl:ccp4_map")
        pyrvapi.rvapi_append_to_data(data, diff_map, "hkl:ccp4_dmap")

    def output_log_files(self, sec, mr_log, ref_log):
        """Function to display the log files for the result

        Parameters
        ----------
        sec : str
            Section the output logs will be added to
        mr_log : str
            Path to the output MR log
        ref_log : str
            Path to the output refinement log

        Returns
        -------
        object
            Section containing mr and refinement logs
        """
        title = "Log files from {0}".format(os.path.basename(mr_log).split('_')[0])

        id = os.path.join(sec, "dat" + str(uuid.uuid4()))
        pyrvapi.rvapi_add_data1(id, title, mr_log, "text", 2, 0, 1, 1, 0)

        id = os.path.join(sec, "dat" + str(uuid.uuid4()))
        pyrvapi.rvapi_add_data1(id, "", ref_log, "text", 2, 0, 1, 1, 0)

    def create_table(self, df, table_id):
        """Function to create/display tables

        Parameters
        ----------
        df : :obj:`~pandas.DataFrame`
            Input  :obj:`~pandas.DataFrame` containing data to be plotted
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
                pyrvapi.rvapi_put_horz_theader(table_id, "PDB_code", self._simbad_tooltips["PDB_code"], 0)
            else:
                pyrvapi.rvapi_put_horz_theader(table_id, l, self._simbad_tooltips[l], i)
                num_labels = i

        ir = len(df)
        for i in range(0, ir):
            for j in range(num_labels + 1):
                if j == 0:
                    pyrvapi.rvapi_put_table_string(table_id, '<a href="http://www.ebi.ac.uk/pdbe/entry/pdb/{0}" '
                                                   'target="_blank">{1}</a>'.format(df.loc[i][j][0:4], df.loc[i][j]), i,
                                                   j)
                else:
                    pyrvapi.rvapi_put_table_string(table_id, str(df.loc[i][j]), i, j)

    @staticmethod
    def create_graphs(df, graph_sec, graph_widget):
        """Function to create/display graphs following MR

        df : :obj:`~pandas.DataFrame`
            Input :obj:`~pandas.DataFrame` containing data to be plotted
        graph_sec : str
            Section the output graph will be displayed in
        graph_widget : str
            Widget ID

        Returns
        -------
        object
            Section containing the graphic representation of Z-score results from SIMBAD
        """
        pyrvapi.rvapi_append_loggraph1(os.path.join(graph_sec, graph_widget))

        pyrvapi.rvapi_add_graph_data1(graph_widget + "/data1", "Scores Vs. Rank (by CC_F Z-score)")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/x", "Rank", "(by CC_F Z-score)")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y1", "CC_F", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y2", "RF_F", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y3", "CC_P", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y4", "CC_I", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y5", "CC_F Z-score", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y6", "CC_P Z-score", "")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y7", "Number_of_rotation_searches_producing_peak", "")

        ir = len(df.index)
        for i in range(0, ir):
            pyrvapi.rvapi_add_graph_int1(graph_widget + "/data1/x", i + 1)
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y1", df['CC_F'][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y2", df['RF_F'][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y3", df["CC_I"][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y4", df["CC_P"][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y5", df["CC_F_Z_score"][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y6", df["CC_P_Z_score"][i], "%g")
            pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y7",
                                          df["Number_of_rotation_searches_producing_peak"][i], "%g")

        # Create a range of graphs
        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot1", "All Z-scores Vs. Rank", "Rank (by CC_F Z-score)",
                                      "Z-Score")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot1", "x", "y5")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot1", "x", "y6")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot2", "CC_F Vs. Rank", "Rank (by CC_F Z-score)", "CC_F")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot2", "x", "y1")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot3", "RF_F Vs. Rank", "Rank (by CC_F Z-score)", "RF_F")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot3", "x", "y2")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot4", "CC_P Vs. Rank", "Rank (by CC_F Z-score)", "CC_P")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot4", "x", "y3")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot5", "CC_I Vs. Rank", "Rank (by CC_F Z-score)", "CC_I")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot5", "x", "y4")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot6", "CC_F Z-score Vs. Rank", "Rank (by CC_F Z-score)",
                                      "CC_F Z-score")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot6", "x", "y5")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot7", "CC_P Z-score Vs. Rank", "Rank (by CC_F Z-score)",
                                      "CC_P Z-score")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot7", "x", "y6")

        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot8", "Freq. of peak /5  Vs. Rank", "Rank (by CC_F Z-score)",
                                      "Freq. of peak /5")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot8", "x", "y7")

    def display_results(self, summarize, results_to_display):

        if self.display_gui or self.ccp4i2:
            if not self.lattice_search_results_displayed:
                lattice_results = os.path.join(self.work_dir, 'latt', 'lattice_search.csv')
                lattice_mr_results = os.path.join(self.work_dir, 'latt', 'lattice_mr.csv')
                if os.path.isfile(lattice_results) or os.path.isfile(lattice_mr_results):
                    self.create_lattice_results_tab(lattice_results, lattice_mr_results, results_to_display)
                    self.lattice_search_results_displayed = True

            if not self.contaminant_results_displayed:
                contaminant_results = os.path.join(self.work_dir, 'cont', 'rot_search.csv')
                contaminant_mr_results = os.path.join(self.work_dir, 'cont', 'cont_mr.csv')
                if os.path.isfile(contaminant_results) or os.path.isfile(contaminant_mr_results):
                    self.create_contaminant_results_tab(contaminant_results, contaminant_mr_results, results_to_display)
                    self.contaminant_results_displayed = True

            if not self.morda_results_displayed:
                morda_db_results = os.path.join(self.work_dir, 'morda', 'rot_search.csv')
                morda_db_mr_results = os.path.join(self.work_dir, 'morda', 'morda_mr.csv')
                if os.path.isfile(morda_db_results) or os.path.isfile(morda_db_mr_results):
                    self.create_morda_db_results_tab(morda_db_results, morda_db_mr_results, results_to_display)
                    self.morda_results_displayed = True

            if summarize:
                self.display_summary_tab()

            pyrvapi.rvapi_flush()

    def save_document(self):
        pyrvapi.rvapi_put_meta(self.rvapi_meta.to_json())
        pyrvapi.rvapi_store_document2(self.rvapi_document)
        pyrvapi.rvapi_keep_polling(True)

    def fix_path(self, path):
        if self.webserver_uri:
            return urlparse.urljoin(self.webserver_uri, path[self._webserver_start:])
        else:
            return path

    def rel_path_for_jscofe(self, path):
        return os.path.join("..", os.path.relpath(path, self.jsrview_dir))

    def adjust_paths_of_files(self, files):
        for f in files:
            if self.jscofe_mode:
                f = self.rel_path_for_jscofe(f)
            yield f

    def store_entry_in_rvapi_meta(self, rank, source, name, pdb, mtz, map_, dmap, best):
        entry = {
            "rank": rank,
            "source": source,
            "best": best,
            "name": name,
            "pdb": pdb,
            "mtz": mtz,
            "map": map_,
            "dmap": dmap
        }
        self.rvapi_meta.add(entry)
