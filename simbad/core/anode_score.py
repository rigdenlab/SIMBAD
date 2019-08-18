"""Class to store ANODE scores"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.1"

from collections import namedtuple

AnomScore = namedtuple("AnomScore", ["dano_peak_height", "nearest_atom"])
