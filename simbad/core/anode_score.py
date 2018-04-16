"""Class to store ANODE scores"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.1"

from simbad.core import ScoreBase


class AnomScore(ScoreBase):
    """An anomalous phased fourier scoring class"""

    __slots__ = ("dano_peak_height", "nearest_atom")

    def __init__(self, dano_peak_height, nearest_atom):
        self.dano_peak_height = dano_peak_height
        self.nearest_atom = nearest_atom

    def __repr__(self):
        return "{0}(dano_peak_height={1} nearest_atom={2})".format(self.__class__.__name__,
                                                                   self.dano_peak_height,
                                                                   self.nearest_atom)

    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.score.anode_score.AnomScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}