from enum import Enum


class MrPrograms(Enum):
    """Container for molecular replacement programs"""
    molrep = 'simbad.mr.molrep_mr'
    phaser = 'simbad.mr.phaser_mr'


class RefPrograms(Enum):
    """Container for refinement programs"""
    refmac5 = 'simbad.mr.refmac_refine'


class SGAlternatives(Enum):
    """Container for space group alternatives"""
    all = 'ALL'
    enant = 'HAND'
    none = 'NONE'