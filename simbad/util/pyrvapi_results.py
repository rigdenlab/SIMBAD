"""Module to interact with pyrvapi"""

import os
import logging
import pandas
import pyrvapi
import subprocess
import uuid
import urlparse

__author__ = "Adam Simpkin"
__date__ = "04 May 2017"
__version__ = "0.1"

logger = logging.getLogger()

class SimbadOutput(object):
    """Class to display the output of SIMBAD"""

    _simbad_tooltips = {"PDB_code" : "The 4 letter code representing the protein in the protein databank",
                        "a" : "Lattice parameter a",
                        "b" : "Lattice parameter b",
                        "c" : "Lattice parameter c", 
                        "alpha" : "Lattice parameter alpha",
                        "beta" : "Lattice parameter beta",
                        "gamma" : "Lattice parameter gamma",
                        "length_penalty" : "The sum of the differences between lattice parameters a, b and c for the model and the target",
                        "angle_penalty" : "The sum of the differences between lattice parameters alpha, beta and gamma for the model and the target",
                        "total_penalty" : "The sum of the length penalty and the angle penalty",
                        "molrep_score" : "MOLREP score for the Molecular Replacement solution",
                        "molrep_tfscore" : "MOLREP translation function score for the Molecular Replacement solution",
                        "phaser_llg" : "PHASER Log-likelihood gain for the Molecular Replacement solution",
                        "phaser_tfz" : "PHASER Translation Function Z-score for the Molecular Replacement solution",
                        "phaser_rfz" : "PHASER Rotational Function Z-score for the Molecular Replacement solution",
                        "final_r_fact" : "R-fact score for REFMAC refinement of the Molecular Replacement solution",
                        "final_r_free" : "R-free score for REFMAC refinement of the Molecular Replacement solution",
                        "peaks_over_6_rms" : "Anomalous peaks over 6 RMS",
                        "peaks_over_6_rms_within_2A_of_model" : "Anomalous peaks over 6 RMS within 2 Angstroms of the Molecular Replacement solution",
                        "peaks_over_12_rms" : "Anomalous peaks over 12 RMS",
                        "peaks_over_12_rms_within_2A_of_model" : "Anomalous peaks over 12 RMS within 2 Angstroms of the Molecular Replacement solution",
                        "ALPHA" : "Lattice parameter alpha",
                        "BETA" : "Lattice parameter beta",
                        "GAMMA" : "Lattice parameter gamma",
                        "CC_F" : "The correlation coefficient between the observed amplitudes for the crystal and the calculated amplitudes for the model",
                        "RF_F" : "The classic R factor between the observed amplitudes for the crystal and the calculated amplitudes for the model.",
                        "CC_I" : "The correlation coefficient between the observed intensities for the crystal and the sum of calculated intensities for all symmetry equivalents of the model",
                        "CC_P" : "The Patterson correlation coefficient between the crystal and the model pattersons evaluated within the defined sphere centred on the Patterson origin.",
                        "Icp" : "",
                        "CC_F_Z_score" : "Z-score of CC_F peaks",
                        "CC_P_Z_score" : "Z-score of CC_P peaks",
                        "Number_of_rotation_searches_producing_peak" : "Number of rotations searches which produce each peak [out of 5]"}

    def __init__(self):
        self.running = None
        self.webserver_uri = None
        self.webserver_start = None
        self.log_tab_id = None
        self.lattice_results_tab_id = None
        self.lattice_df = None
        self.contaminant_results_tab_id = None
        self.contaminant_df = None
        self.morda_db_results_tab_id = None
        self.morda_db_df = None
        self.summary_tab_id = None
        self.summary_tab_results_sec_id = None
        return

    def create_log_tab(self, logfile):
        if self.log_tab_id:
            return
        if not os.path.isfile(logfile):
            return False

        self.log_tab_id = "log_tab"
        logurl = self.fix_path(logfile)
        pyrvapi.rvapi_add_tab(self.log_tab_id, "Log file", True)

        # Add watched (updatable) content to the log tab
        pyrvapi.rvapi_append_content(logurl, True, self.log_tab_id)

        return self.log_tab_id

    def _create_lattice_results_tab(self):
        if not self.lattice_results_tab_id:
            self.lattice_results_tab_id = "lattice_results_tab"
            pyrvapi.rvapi_add_tab(self.lattice_results_tab_id, "Lattice Parameter Search Results", False)
        return
    
    def _create_contaminant_results_tab(self):
        if not self.contaminant_results_tab_id:
            self.contaminant_results_tab_id = "contaminants_results_tab"
            pyrvapi.rvapi_add_tab(self.contaminant_results_tab_id, "Contaminant Search Results", False)
        return
    
    def _create_morda_db_results_tab(self):
        if not self.morda_db_results_tab_id:
            self.morda_db_results_tab_id = "morda_db_results_tab"
            pyrvapi.rvapi_add_tab(self.morda_db_results_tab_id, "MoRDa Database Search Results", False)
        return
    
    def _create_summary_tab(self):
        if not self.summary_tab_id:
            self.summary_tab_id = "summary_tab"
            pyrvapi.rvapi_add_tab(self.summary_tab_id, "Summary", True)
        return

    def create_lattice_results_tab(self, work_dir, lattice_results, lattice_mr_results):
        """Function to create the lattice results tab"""
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
            
            section_title = 'Top 10 Lattice Parameter Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)
            
            section_title = 'Top 10 Lattice Parameter Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)
            
            self.lattice_df = df
            
            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    run_dir = os.path.join(work_dir, 'jsrview')
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(work_dir, 'latt', 'mr_lattice', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
                    
                    self.output_result_files(download_sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)
                    
                    
                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))
                    pass
        return
    
    def create_contaminant_results_tab(self, work_dir, contaminant_results, contaminant_mr_results):
        """Function to create the contaminant results tab"""
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
            
            section_title = "Molecular Replacement Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)
            
            section_title = 'Top 10 Contaminant Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)
            
            section_title = 'Top 10 Contaminant Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)
            
            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    run_dir = os.path.join(work_dir, 'jsrview')
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(work_dir, 'cont', 'mr_contaminant', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
                    
                    self.output_result_files(download_sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)
                    
                    
                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))
                    pass
            
        return
    
    def create_morda_db_results_tab(self, work_dir, morda_db_results, morda_db_mr_results):
        """Function to create the MoRDa Database results tab"""
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
             
            self.morda_db_df  =df
             
            section_title = "Molecular Replacement Search Graphs"
            uid = str(uuid.uuid4())
            graph_sec = section_title.replace(" ", "_") + uid
            graph_widget = "graphWidget" + uid
            pyrvapi.rvapi_add_section(graph_sec, section_title, tab, 0, 0, 1, 1, True)
            self.create_graphs(df, graph_sec, graph_widget)
             
            section_title = 'Top 10 MoRDa database Search Downloads'
            uid = str(uuid.uuid4())
            download_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(download_sec, section_title, tab, 0, 0, 1, 1, True)
             
            section_title = 'Top 10 MoRDa database Search Log Files'
            uid = str(uuid.uuid4())
            logfile_sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(logfile_sec, section_title, tab, 0, 0, 1, 1, False)
             
            for i in range(0, 10):
                try:
                    pdb_code = df.loc[i][0]
                    run_dir = os.path.join(work_dir, 'jsrview')
                    mr_program = list(df)[1][0:6]
                    mr_workdir = os.path.join(work_dir, 'full', 'mr_full', pdb_code, 'mr', mr_program)
                    mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                    ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                    ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                    ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                    ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                    diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
                     
                    self.output_result_files(download_sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb)
                    self.output_log_files(logfile_sec, mr_log, ref_log)
                     
                     
                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))
                    pass
             
        return
    
    def create_summary_tab(self, work_dir):
        """Function to create the MoRDa Database results tab"""
        self._create_summary_tab()
        
        if self.lattice_df is None:
            lattice_score = 1
        else:
            lattice_score = self.lattice_df['final_r_free'][0]
        if self.contaminant_df is None:
            contaminant_score = 1
        else:
            contaminant_score = self.contaminant_df['final_r_free'][0]
        if self.morda_db_df is None:
            morda_db_score = 1
        else:
            morda_db_score = self.morda_db_df['final_r_free'][0]
        
        if lattice_score <= contaminant_score and lattice_score <= morda_db_score:
            pdb_code = self.lattice_df.loc[0][0]
            r_fact = self.lattice_df['final_r_fact'][0]
            r_free = self.lattice_df['final_r_free'][0]
            run_dir = os.path.join(work_dir, 'jsrview')
            mr_program = list(self.lattice_df)[1][0:6]
            mr_workdir = os.path.join(work_dir, 'latt', 'mr_lattice', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
            
        elif contaminant_score <= lattice_score and contaminant_score <= morda_db_score:
            pdb_code = self.contaminant_df.loc[0][0]
            r_fact = self.contaminant_df['final_r_fact'][0]
            r_free = self.contaminant_df['final_r_free'][0]
            run_dir = os.path.join(work_dir, 'jsrview')
            mr_program = list(self.contaminant_df)[1][0:6]
            mr_workdir = os.path.join(work_dir, 'cont', 'mr_contaminant', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
            
        elif morda_db_score <= lattice_score and morda_db_score <= contaminant_score:
            pdb_code = self.morda_db_df.loc[0][0]
            r_fact = self.morda_db_df['final_r_fact'][0]
            r_free = self.morda_db_df['final_r_free'][0]
            run_dir = os.path.join(work_dir, 'jsrview')
            mr_program = list(self.morda_db_df)[1][0:6]
            mr_workdir = os.path.join(work_dir, 'full', 'mr_full', pdb_code, 'mr', mr_program)
            mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
            ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
            ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
            ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
            ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
            diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
            
        else:
            logger.debug('Unexpected result')
            return
            
        section_title = 'SIMBAD Summary'
        uid = str(uuid.uuid4())
        sec = section_title.replace(" ", "_") + uid
        tab = self.summary_tab_id
        
        msg = 'The best search model found by SIMBAD was {0}. \
               This gave an R/Rfact of {1} and an R/Rfree of {2}. \
               An R/Rfree lower than 0.45 is indicative of a \
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
        
        self.output_result_files(download_sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb)
        self.output_log_files(logfile_sec, mr_log, ref_log)
        return

    def create_table(self, df, table_id):
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
                    pyrvapi.rvapi_put_table_string(table_id, 
                                                   '<a href="http://www.ebi.ac.uk/pdbe/entry/pdb/{0}" '
                                                   'target="_blank">{1}</a>'.format(df.loc[i][j][0:4], 
                                                                                    df.loc[i][j]), i, j)
                else:
                    pyrvapi.rvapi_put_table_string(table_id, str(df.loc[i][j]), i, j)
                    
    def create_graphs(self, df, graph_sec, graph_widget):
        pyrvapi.rvapi_append_loggraph1(os.path.join(graph_sec, graph_widget))
        
        pyrvapi.rvapi_add_graph_data1(graph_widget + "/data1"," Scores Vs. Rank (by R-Free)")
        pyrvapi.rvapi_add_graph_dataset1 ( graph_widget + "/data1/x","Rank"," (by R-Free)" )
        
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y1","REFMAC R-Fact","")
        pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y2","REFMAC R-Free","")
        
        mr_program = list(df)[1][0:6]
        if mr_program == 'molrep':
            pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y3","MOLREP score","")
            pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y4","MOLREP TF/sig","")
            
        elif mr_program == 'phaser':
            pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y3","PHASER TFZ","")
            pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y4","PHASER LLG","")
            pyrvapi.rvapi_add_graph_dataset1(graph_widget + "/data1/y5","PHASER RFZ","")
            
        
        ir = len(df.index)
        for i in range(0, ir):
            if df['final_r_free'][i] < 0.7 and df['final_r_fact'][i] < 0.7:
                    pyrvapi.rvapi_add_graph_int1(graph_widget + "/data1/x", i + 1)
                    pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y1", df['final_r_fact'][i],"%g")
                    pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y2", df['final_r_free'][i],"%g")
                    if mr_program == 'molrep':
                        pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y3", df['molrep_score'][i],"%g")
                        pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y4", df['molrep_tfscore'][i],"%g")
                    
                    elif mr_program == 'phaser':
                        pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y3", df['phaser_tfz'][i],"%g")
                        pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y4", df['phaser_llg'][i],"%g")
                        pyrvapi.rvapi_add_graph_real1(graph_widget + "/data1/y5", df['phaser_rfz'][i],"%g")
            
        
        pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot1","R-Fact/R-Free Vs. Rank", "Rank (by R-Free)", "R-Fact/R-Free")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot1","x","y1")
        pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot1","x","y2")
        pyrvapi.rvapi_set_plot_xmin("plot1","graphWidget1",-1.0)
     
        if mr_program == 'molrep':
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot2","MOLREP score Vs. Rank", "Rank (by R-Free)", "MOLREP score")
            pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot2","x","y3")
         
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot3","MOLREP TF/sig Vs. Rank", "Rank (by R-Free)", "MOLREP TF/sig")
            pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot3","x","y4")
            
        elif mr_program == 'phaser':
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot2","PHASER TFZ Vs. Rank", "Rank (by R-Free)", "PHASER TFZ")
            pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot2","x","y3")
         
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot3","PHASER LLG Vs. Rank", "Rank (by R-Free)", "PHASER LLG")
            pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot3","x","y4")
            
            pyrvapi.rvapi_add_graph_plot1(graph_widget + "/plot4","PHASER RFZ Vs. Rank", "Rank (by R-Free)", "PHASER RFZ")
            pyrvapi.rvapi_add_plot_line1(graph_widget + "/data1/plot4","x","y5")
        return
                    
    def output_result_files(self, sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb):
        uglymol = self.uglymol_html(run_dir, ref_pdb, ref_map, diff_map)
        
        title = "Electron density for {0}".format(os.path.basename(ref_pdb).split('_')[0])
        uid = str(uuid.uuid4())
        data = "dat" + uid
        
        pyrvapi.rvapi_add_data1(os.path.join(sec, data), title, ref_pdb, "xyz",2,0,1,1,1)
        pyrvapi.rvapi_append_to_data(data, ref_mtz, "hkl:map" )
    
        pyrvapi.rvapi_add_text('<a href="' + uglymol + '">Click here for Uglymol visualisation</a><br></br>',sec,2,0,1,1 )
        return
    
    def output_log_files(self, sec, mr_log, ref_log):
        title = "Log files from {0}".format(os.path.basename(mr_log).split('_')[0])
        uid = str(uuid.uuid4())
        data = "dat" + uid
        pyrvapi.rvapi_add_data1 (os.path.join(sec, data), title, mr_log, "text",2,0,1,1,0)
        uid = str(uuid.uuid4())
        data = "dat" + uid
        pyrvapi.rvapi_add_data1 (os.path.join(sec, data),"", ref_log,"text",2,0,1,1,0 )
        return

    def display_results(self, webserver_uri, no_gui, logfile, work_dir=None, summary=False):
        if no_gui:
            return

        if not self.running:
            # Infrastructure to run
            ccp4 = os.environ["CCP4"]
            share_jsrview = os.path.join(ccp4, "share", "jsrview")

            
            if not work_dir:
                logger.debug("Please specify working directory")
                return
            
            run_dir = os.path.join(work_dir, 'jsrview')
            if not os.path.isdir(run_dir):
                os.mkdir(run_dir)

            pyrvapi.rvapi_init_document("SIMBAD_results", run_dir, "SIMBAD Results", 1, 7, share_jsrview, None, None,
                                        None, None)
            if webserver_uri:
                # don't start browser and setup variables for the path on the webserver
                self._webserver_start = len(run_dir) + 1
                self.webserver_uri = webserver_uri
            else:
                # We start our own browser
                jsrview = os.path.join(ccp4, "libexec", "jsrview")
                subprocess.Popen([jsrview, os.path.join(run_dir, "index.html")])
            pyrvapi.rvapi_add_header("SIMBAD Results")
            self.running = True

        if os.path.isfile(logfile):
            self.create_log_tab(logfile)
        
        lattice_results = os.path.join(work_dir, 'latt/lattice_search.csv')
        lattice_mr_results = os.path.join(work_dir, 'latt/lattice_mr.csv')
        if os.path.isfile(lattice_results) or os.path.isfile(lattice_mr_results):
            self.create_lattice_results_tab(work_dir, lattice_results, lattice_mr_results)
            
        contaminant_results = os.path.join(work_dir, 'cont/rot_search.csv')
        contaminant_mr_results = os.path.join(work_dir, 'cont/cont_mr.csv')
        if os.path.isfile(contaminant_results) or os.path.isfile(contaminant_mr_results):
            self.create_contaminant_results_tab(work_dir, contaminant_results, contaminant_mr_results)
            
        morda_db_results = os.path.join(work_dir, 'full/rot_search.csv')
        morda_db_mr_results = os.path.join(work_dir, 'full/full_mr.csv')
        if os.path.isfile(morda_db_results) or os.path.isfile(morda_db_mr_results):
            self.create_morda_db_results_tab(work_dir, morda_db_results, morda_db_mr_results)
            
        if summary:
            self.create_summary_tab(work_dir)

        pyrvapi.rvapi_flush()

        return True

    def fix_path(self, path):
        """Ammend path so it's suitable for the webserver"""
        if self.webserver_uri:
            return urlparse.urljoin(self.webserver_uri, path[self._webserver_start:])
        else:
            return path
        
    def uglymol_html(self, run_dir, ref_pdb, ref_map, diff_map):
        html_out = os.path.join(run_dir, os.path.basename(ref_pdb).split('_')[0] + ".html")
        with open(html_out, "w") as f:
            f.write("""<!doctype html>
    <html lang="en">
    <head>
      <title>1mru - UglyMol</title>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, user-scalable=no">
      <meta name="theme-color" content="#333333">
      <style>
       canvas { display: block; }
       #hud {
         font: 14px sans-serif;
         color: #ddd;
         background-color: rgba(0,0,0,0.6);
         text-align: center;
         position: absolute;
         top: 10px;
         left: 50%;
         transform: translateX(-50%);
         padding: 2px 8px;
         border-radius: 5px;
         z-index: 9;
         white-space: pre-line;
       }
       #help {
         display: none;
         font: 16px sans-serif;
         color: #eee;
         background-color: rgba(0,0,0,0.7);
         position: absolute;
         left: 20px;
         top: 50%;
         transform: translateY(-50%);
         cursor: default;
         padding: 5px;
         border-radius: 5px;
         z-index: 9;
         white-space: pre-line;
       }
       #inset {
         width: 200px;
         height: 200px;
         background-color: #888;
         position: absolute;
         right: 0;
         bottom: 0;
         z-index: 2;
         display: none;
       }
       a { color: #59C; }
      </style>
    
      <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r79/three.min.js"></script>
      <script src="../output/uglymol.js"></script>
    
    </head>
    <body style="background-color: black">
      <div id="viewer" style="position: absolute; left:0px; top:0px;">
          <header id="hud" onmousedown="event.stopPropagation();"
                           onmousemove="event.stopPropagation();"
                           ondblclick="event.stopPropagation();"
                 >This is uglymol not coot. <a href="#"
                             onclick="V.toggle_help(); return false;"
                             >H shows help.</a></header>
        <footer id="help"></footer>
      </div>
      <div id="inset"></div>
      <script>
        V = new Viewer("viewer");""")
            f.write("""V.load_pdb("{pdb}");
        V.load_ccp4_maps("{map}", "{diff_map}");
        //V.show_nav("inset");
        V.render();
      </script>
    </body>
    </html>""".format(pdb=ref_pdb, map=ref_map, diff_map=diff_map))
            
        return self.fix_path(html_out)

