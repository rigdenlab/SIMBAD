"""Class to store AMORE rotation scores"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "10 Oct 2017"
__version__ = "0.1"

from simbad.score import ScoreBase


class AmoreRotationScore(ScoreBase):
    """An amore rotation scoring class"""

    __slots__ = ("pdb_code", "dat_path", "ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
                 "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak")

    def __init__(self, pdb_code, dat_path, ALPHA, BETA, GAMMA, CC_F, RF_F, CC_I, CC_P, Icp,
                 CC_F_Z_score, CC_P_Z_score, Number_of_rotation_searches_producing_peak):
        self.pdb_code = pdb_code
        self.dat_path = dat_path
        self.ALPHA = ALPHA
        self.BETA = BETA
        self.GAMMA = GAMMA
        self.CC_F = CC_F
        self.RF_F = RF_F
        self.CC_I = CC_I
        self.CC_P = CC_P
        self.Icp = Icp
        self.CC_F_Z_score = CC_F_Z_score
        self.CC_P_Z_score = CC_P_Z_score
        self.Number_of_rotation_searches_producing_peak = Number_of_rotation_searches_producing_peak

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code} dat_path={dat_path} " \
                 "ALPHA={ALPHA} BETA={BETA} GAMMA={GAMMA} " \
                 "CC_F=CC_F RF_F={RF_F} CC_I={CC_I} CC_P={CC_P} Icp={Icp} " \
                 "CC_F_Z_score={CC_F_Z_score} CC_P_Z_score={CC_P_Z_score} " \
                 "Number_of_rotation_searches_producing_peak={Number_of_rotation_searches_producing_peak})"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`AmoreRotationScore <simbad.rotsearch.amore_score.AmoreRotationScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}
