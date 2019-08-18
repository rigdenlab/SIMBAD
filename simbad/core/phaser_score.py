"""Class to store PHASER rotation scores"""

__author__ = "Adam Simpkin"
__date__ = "27 Dec 2017"
__version__ = "0.1"

from collections import namedtuple

PhaserRotationScore = namedtuple("PhaserRotationScore", ["pdb_code", "dat_path", "llg", "rfz"])
