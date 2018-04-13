"""Class to store MR scores"""

__author__ = "Adam Simpkin"
__date__ = "17 Oct 2017"
__version__ = "0.1"


class MrScore(object):
    """A molecular replacement scoring class"""

    __slots__ = ("pdb_code", "final_r_fact", "final_r_free", "molrep_score", "molrep_tfscore",
                 "phaser_tfz", "phaser_llg", "phaser_rfz", "dano_peak_height", "nearest_atom")

    def __init__(self, pdb_code):
        self.pdb_code = pdb_code
        self.molrep_score = None
        self.molrep_tfscore = None
        self.phaser_tfz = None
        self.phaser_llg = None
        self.phaser_rfz = None
        self.final_r_fact = 1.0
        self.final_r_free = 1.0
        self.dano_peak_height = None
        self.nearest_atom = None

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code} final_r_fact={final_r_fact} final_r_free={final_r_free}"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.mr.mr_score.MrScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}
