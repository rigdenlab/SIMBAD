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


class SimbadOutput(object):
    """Class to display the output of SIMBAD"""

    # To fill once I know exactly what I'll be displaying
    # Perhaps need lattice, contaminant and full tool tips
    _simbad_tooltips = {}

    def __init__(self):
        self.running = None
        self.webserver_uri = None
        self.webserver_start = None
        self.log_tab_id = None
        self.lattice_results_tab_id = None
        self.lattice_results_tab_sections = []
        self.contaminant_results_tab_id = None
        self.contaminant_results_tab_sections = []
        self.morda_db_results_tab_id = None
        self.morda_db_results_tab_sections = []
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

    def create_lattice_results_tab(self, lattice_results, lattice_mr_results):
        """Note: results = lattice.csv"""

        # if not self.summary_tab_id:
        #     return

        if not lattice_results:
            return

        self._create_lattice_results_tab()

        if lattice_results:
            section_title = 'Lattice Parameter Search'
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
            section_title = 'Lattice Parameter Search'
            uid = str(uuid.uuid4())

            sec = section_title.replace(" ", "_") + uid
            tab = self.lattice_results_tab_id
            table = "table" + uid

            pyrvapi.rvapi_add_section(sec, section_title, tab, 0, 0, 1, 1, True)

            table_title = "Molecular Replacement Search Results"
            pyrvapi.rvapi_add_table1(sec + "/" + table, table_title, 2, 0, 1, 1, 100)
            df = pandas.read_csv(lattice_mr_results)
            self.create_table(df, table)
        return

    @staticmethod
    def create_table(df, table_id):
        for i, l in enumerate(df):
            if i == 0:
                pyrvapi.rvapi_put_horz_theader(table_id, "PDB code", "", 0)
            else:
                pyrvapi.rvapi_put_horz_theader(table_id, l, "", i)
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

    def display_results(self, webserver_uri, no_gui, logfile, lattice_results=None, lattice_mr_results=None, run_dir=None):
        if no_gui:
            return

        logger = logging.getLogger()

        if not self.running:
            # Infrastructure to run
            ccp4 = os.environ["CCP4"]
            share_jsrview = os.path.join(ccp4, "share", "jsrview")

            if not run_dir:
                logger.debug("Please specify running directory for jsrview")
                return

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
        self.create_lattice_results_tab(lattice_results, lattice_mr_results)

        pyrvapi.rvapi_flush()

        return True

    def fix_path(self, path):
        """Ammend path so it's suitable for the webserver"""
        if self.webserver_uri:
            return urlparse.urljoin(self.webserver_uri, path[self._webserver_start:])
        else:
            return path

if __name__ == "__main__":
    SR = SimbadOutput()

    webserver_uri = False
    no_gui = False
    logfile = '/home/experiences/proxima2a/simpkin/dev/simbad_python/test/SIMBAD_RUN/SIMBAD.log'
    lattice_results = '/home/experiences/proxima2a/simpkin/dev/simbad_python/test/SIMBAD_RUN/lattice.csv'
    lattice_mr_results = '/home/experiences/proxima2a/simpkin/dev/simbad_python/test/SIMBAD_RUN/mr_results.csv'
    run_dir = os.path.abspath(os.path.join(os.curdir, "pyrvapi_tmp"))

    SR.display_results(webserver_uri, no_gui, logfile, lattice_results, lattice_mr_results, run_dir)
