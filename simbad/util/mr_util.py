'''
Class to run MR on SIMBAD results using code from MrBump

@author: hlasimpk
'''

import os
import sys

# Set up MrBUMP imports
if os.path.isdir(os.path.join(os.environ["CCP4"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
mrbump_incl = os.path.join(mrbump, "include")

sys.path.append(os.path.join(mrbump_incl, 'dev'))
sys.path.append(os.path.join(mrbump_incl, 'file_info'))
sys.path.append(os.path.join(mrbump_incl, 'modelling'))
sys.path.append(os.path.join(mrbump_incl, 'mr'))
sys.path.append(os.path.join(mrbump_incl, 'output'))
sys.path.append(os.path.join(mrbump_incl, 'seq_align'))
sys.path.append(os.path.join(mrbump_incl, 'ccp4'))
sys.path.append(os.path.join(mrbump_incl, 'structures'))
sys.path.append(os.path.join(mrbump_incl, 'tools'))
sys.path.append(os.path.join(mrbump_incl, 'initialisation'))
sys.path.append(os.path.join(mrbump_incl, 'building'))
sys.path.append(os.path.join(mrbump_incl, 'cluster'))
sys.path.append(os.path.join(mrbump_incl, 'dispatchers'))
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
sys.path.append(os.path.join(mrbump_incl, 'phasing'))

import MRBUMP_initialise
import MRBUMP_master
import molrepEXE

import phaser_util

# Needs to be written, just setting up how I want this to run

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