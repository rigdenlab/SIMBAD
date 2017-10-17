"""Class to store MR scores"""

__author__ = "Adam Simpkin"
__date__ = "17 Oct 2017"
__version__ = "0.1"


class MrScore(object):
    """A molecular replacement scoring class"""

    __slots__ = ("pdb_code", "final_r_fact", "final_r_free", "molrep_score", "molrep_tfscore",
                 "phaser_tfz", "phaser_llg", "phaser_rfz", "peaks_over_6_rms", "peaks_over_6_rms_within_4a_of_model",
                 "peaks_over_9_rms", "peaks_over_9_rms_within_4a_of_model")

    def __init__(self, pdb_code):
        self.pdb_code = pdb_code
        self.molrep_score = None
        self.molrep_tfscore = None
        self.phaser_tfz = None
        self.phaser_llg = None
        self.phaser_rfz = None
        self.final_r_fact = 1.0
        self.final_r_free = 1.0
        self.peaks_over_6_rms = None
        self.peaks_over_6_rms_within_4a_of_model = None
        self.peaks_over_9_rms = None
        self.peaks_over_9_rms_within_4a_of_model = None

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code}  final_r_fact={final_r_fact} final_r_free={final_r_free}"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`_MrScore <simbad.mr._MrScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}