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

import cluster_run

import MRBUMP_initialise
import MRBUMP_master

import MRBUMP_Molrep
import MRBUMP_Phaser

"""Currently taking reading through Ronan's MR script in MrBump and determining the input variables needed by MRBUMP_Molrep/Phaser"""

class MR_submit(object):
    '''Class to submit MR job'''

    def __init__(self, optd):
        """NEED TO DECIDE WHAT TO DO HERE"""

        return

    def MR_is_done(self, init, mstat):
        """ A function that returns True if there are no more jobs to submit"""

        # NEED TO GIVE A SORTED MODEL LIST
        for model in mstat.sorted_model_list:

            # NEED TO GIVE MR PROGRAM TO USE
            for MRPROGRAM in init.keywords.MR_PROGRAM_LIST:

                # Set the MR program : NEED TO FIND OUT WHAT MSTAT IS HERE
                mstat.setMRPROGRAM(MRPROGRAM.upper())

                # DICTIONARY? Containing model path?
                MODEL = mstat.model_list[model]

                if mstat.MRPROGRAM == "MOLREP" and not MODEL.isMolrepSubmitted() or \
                   mstat.MRPROGRAM == "PHASER" and not MODEL.isPhaserSubmitted():
                    return False
            return True




class MR_cluster_submit(object):
    '''Class to run MR on a cluster'''

    def __init__(self, optd):
        self.input_file = None
        self.mr_keyfile = None
        self.mr_program = None
        self.refine_keyfile = None
        self.refine_program = "refmac5"

        self.parse_options(optd)

        cljob = cluster_run.ClusterJob()

        if self.input_file:
            # Parse input file
            cljob.parse_input(self.input_file)

        if self.mr_keyfile and self.mr_program and self.refine_keyfile and self.refine_program:
            # Run cluster job
            cljob.run(self.mr_program, self.refine_program, self.mr_keyfile, refine_keyfile=self.refine_keyfile)

    def parse_options(self, optd):
        '''Function to set up input files for the MR job'''

        # NEED TO FIND OUT ALL THE OPTIONS THAT NEED TO BE INPUT HERE AND ADD LOGIC TO DECIDE WHICH FILES TO INPUT BASED ON TYPE OF RUN

        self.input_file = None

        # input file needs to contain a list of input settings including:
        # DIRE = directory
        # SGIN = Space group in
        # HKL1 = HKL 1 input file
        # HKLR = HKL refined input file (file name HKL1 gets copied to)

        # If MR program == phaser
        # HKLO = phaser HKL output file

        # PDBO = MR PDB output
        # MRLO = MR log

        # REFH = refinement HKL out
        # REFP = refinement PDB out

        # ENAN = try enantiomorph forms < True | False >

        # FPIN = FP flag
        # SIGF = SIGF flag
        # FREE = Free flag
        # SOLV = solvent content
        # RESO = resolution

        # PDBI = PDB input


        self.mr_program = optd.d['MR_program']
        self.refine_program = optd.d['refine_program']

        # Need to find out what the mr/refine key files are and add them in here
        self.mr_keyfile = ""
        # Just a file containing MOLREP key words - full path needs to be given, can be empty for now

        self.refine_keyfile = ""
        # Just a file containing PHASER key words - full path needs to be given, can be empty for now

        return




# POTENTIAL CODE TO ADD TO __INIT__

# if type == "lattice":
#     self.results_file = os.path.join(optd.d["work_dir"], "lattice.csv")
#     self.structure_dir = os.path.join(optd.d["work_dir"], "lattice_input_models")
#     self.ncyc = optd.d["ncyc_lattice"]
#
# elif type == "contaminant":
#     self.results_file = os.path.join(optd.d["work_dir"], "contaminant.csv")
#
#     # Check to see if user defined contaminant database present, if not use the default
#     if optd.d["Contaminant_database"]:
#         self.structure_dir = optd.d["Contaminant_database"]
#     else:
#         self.structure_dir = CONTAMINANT_MODELS
#     self.ncyc = optd.d["ncyc_contam"]
#
# elif type == "full":
#     self.results_file = os.path.join(optd.d["work_dir"], "full.csv")
#
#     if optd.d["MoRDa_database"]:
#         self.structure_dir = optd.d["MoRDa_database"]
#         self.ncyc = optd.d["ncyc_full"]
#     else:
#         msg = "Path to MoRDa database not specified"
#         raise RuntimeError(msg)
#         exit_util.exit_error(msg)
# else:
#     msg = "Type of MOLREP run not specified"
#     raise RuntimeError(msg)
#     exit_util.exit_error(msg)

