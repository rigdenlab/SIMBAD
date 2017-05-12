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
                        "phaser_tgz" : "PHASER Translation Function Z-score for the Molecular Replacement solution",
                        "phaser_rfz" : "PHASER Rotational Function Z-score for the Molecular Replacement solution",
                        "final_r_fact" : "R-fact score for REFMAC refinement of the Molecular Replacement solution",
                        "final_r_free" : "R-free score for REFMAC refinement of the Molecular Replacement solution",
                        "peaks_over_6_rms" : "Anomalous peaks over 6 RMS",
                        "peaks_over_6_rms_within_2A_of_model" : "Anomalous peaks over 6 RMS within 2 Angstroms of the Molecular Replacement solution",
                        "peaks_over_12_rms" : "Anomalous peaks over 12 RMS",
                        "peaks_over_12_rms_within_2A_of_model" : "Anomalous peaks over 12 RMS within 2 Angstroms of the Molecular Replacement solution"}

    def __init__(self):
        self.running = None
        self.webserver_uri = None
        self.webserver_start = None
        self.log_tab_id = None
        self.lattice_results_tab_id = None
        self.contaminant_results_tab_id = None
        self.morda_db_results_tab_id = None
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
            pyrvapi.rvapi_insert_tab(self.lattice_results_tab_id, "Lattice Parameter Search Results", 
                                     self.summary_tab_id, False)
        return

    def create_lattice_results_tab(self, work_dir, lattice_results, lattice_mr_results):
        """Note: results = lattice.csv"""

        # if not self.summary_tab_id:
        #     return

        if not lattice_results:
            return

        self._create_lattice_results_tab()

        if lattice_results:
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

        if lattice_mr_results:
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
            
            
            section_title = 'Lattice Parameter Search Downloads'
            uid = str(uuid.uuid4())
            sec = section_title.replace(" ", "_") + uid
            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)
            
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
                    
                    # Need to get make these files
                    ref_map = None
                    diff_map = None
                    
                    self.output_result_files(sec, run_dir, diff_map, ref_map, ref_mtz, ref_pdb)
                    
                    section_title = 'Lattice Parameter Search Log Files'
                    uid = str(uuid.uuid4())
                    sec = section_title.replace(" ", "_") + uid
                    pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, False)
                    
                    self.output_log_files(sec, mr_log, ref_log)
                    
                    
                except KeyError:
                    logger.debug("No result found at position %s", (i + 1))
                    pass
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

    def display_results(self, webserver_uri, no_gui, logfile, lattice_results=None, lattice_mr_results=None, work_dir=None):
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

        self.create_log_tab(logfile)
        self.create_lattice_results_tab(work_dir, lattice_results, lattice_mr_results)

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

if __name__ == "__main__":
    SR = SimbadOutput()

    webserver_uri = False
    no_gui = False
    logfile = '/Users/adamsimpkin/dev/test/SIMBAD_8/debug.log'
    lattice_results = '/Users/adamsimpkin/dev/test/SIMBAD_8/latt/lattice_search.csv'
    lattice_mr_results = '/Users/adamsimpkin/dev/test/SIMBAD_8/latt/lattice_mr.csv'
    work_dir = os.path.abspath(os.path.join(os.curdir, "SIMBAD_8"))

    SR.display_results(webserver_uri, no_gui, logfile, lattice_results, lattice_mr_results, work_dir)
