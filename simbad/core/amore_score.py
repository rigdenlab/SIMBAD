"""Class to store AMORE rotation scores"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "10 Oct 2017"
__version__ = "0.1"

from collections import namedtuple


AmoreRotationScore = namedtuple(
    "AmoreRotationScore",
    [
        "pdb_code",
        "dat_path",
        "ALPHA",
        "BETA",
        "GAMMA",
        "CC_F",
        "RF_F",
        "CC_I",
        "CC_P",
        "Icp",
        "CC_F_Z_score",
        "CC_P_Z_score",
        "Number_of_rotation_searches_producing_peak",
    ],
)
