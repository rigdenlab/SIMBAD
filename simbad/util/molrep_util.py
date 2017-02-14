'''
Wrapper to run MOLREP for molecular replacement

@author: hlasimpk

14/02/2017
'''

import os
import shutil
import logging

from simbad.util import simbad_util
from simbad.util import exit_util
from simbad.constants import CONTAMINANT_MODELS


### THIS NEEDS EDITING TO SUIT OUR PURPOSES###



if not "CCP4" in os.environ.keys():
    raise RuntimeError('CCP4 not found')

LOGGER = logging.getLogger(__name__)

class Molrep(object):
    '''A class to submit Molrep jobs to the local machine or a batch system on a cluster'''

    def __init__(self, optd=None, type=None):

        if optd:
            self.init(optd, type)

        return

    def init(self, optd, type):
        self.mtz = optd['mtz']

        if type == "lattice":
            self.results_file = os.path.join(optd.d["work_dir"], "lattice.csv")
            self.structure_dir = os.path.join(optd.d["work_dir"], "lattice_input_models")
            self.ncyc = optd.d["ncyc_lattice"]

        elif type == "contaminant":
            self.results_file = os.path.join(optd.d["work_dir"], "contaminant.csv")

            # Check to see if user defined contaminant database present, if not use the default
            if optd.d["Contaminant_database"]:
                self.structure_dir = optd.d["Contaminant_database"]
            else:
                self.structure_dir = CONTAMINANT_MODELS
            self.ncyc = optd.d["ncyc_contam"]

        elif type == "full":
            self.results_file = os.path.join(optd.d["work_dir"], "full.csv")

            if optd.d["MoRDa_database"]:
                self.structure_dir = optd.d["MoRDa_database"]
                self.ncyc = optd.d["ncyc_full"]
            else:
                msg = "Path to MoRDa database not specified"
                raise RuntimeError(msg)
                exit_util.exit_error(msg)
        else:
            msg = "Type of MOLREP run not specified"
            raise RuntimeError(msg)
            exit_util.exit_error(msg)

        return

    def set_logfile(self, filename):
        self.logfile = filename
        return

    def set_molrep_keywords(self, init, target_info, model, RFfile, SEARCH=""):

        ''' Set up the MOLREP keywords'''

        # Set the key words
        model.molrep_keywords["DOC"]="Y"
        model.molrep_keywords["LABIN"]='F=' + target_info.mtz_coldata['F'] \
                                       + ' SIGF=' + target_info.mtz_coldata['SIGF']
        if model.COPIES_SET:
            model.molrep_keywords["NMON"] = `model.no_copies`
        else:
            model.molrep_keywords["NMON"] = `target_info.no_of_mols/model.number_mols`

        model.molrep_keywords["SIM"] = `model.seqID[0]`

        # Use only basename for RFfile to stop Molrep crashing with long paths
        if os.path.isfile(os.path.join(model.model_directory, "mr", "molrep", os.path.basename(RFfile))) \
                and RFfile != os.path.join(model.model_directory, "mr", "molrep", os.path.basename(RFfile)):
            shutil.copyfile(RFfile, os.path.join(model.model_directory, "mr", "molrep", os.path.basename(RFfile)))
        model.molrep_keywords["FILE_T"] = os.path.basename(RFfile)

        # Add the spacegroup to search under if we are looking at the enantiomorph
        if SEARCH == "ENANT":
            model.molrep_keywords["NOSG"] = target_info.enant_SG_code

        # If a fixed component of the solution was provided we need to add the MODEL_2 keyword
        if init.keywords.FIXED:
            model.molrep_keywords["MODEL_2"] = init.molrep_fixed_PDBfile

    def run_molrep(self, optd):
        '''Run MOLREP to get MR solution for the top 20 lattice results'''
        count = 0
        with open(self.results_file, 'r') as f:
            for line in f:
                if count < self.jobs:
                    PDB_code = line.split(',')[0]


        return

    def molrep_cmd(self, model):
        '''Return the command to run MOLREP'''

        molrep = os.path.join(os.environ["CCP4"], "bin", "molrep")

        cmd =[
            molrep,
            "-f", "{0}".format(self.mtz),
            "-m", "{0}".format(model),
        ]

        return " ".join(cmd)


