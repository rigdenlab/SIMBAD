"""Class to store dat file info"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.1"

from collections import namedtuple

DatModelScore = namedtuple("DatModelScore", ["pdb_code", "dat_path", "mw_diff", "x", "y", "z", "intrad", "solvent", "nmol"])
