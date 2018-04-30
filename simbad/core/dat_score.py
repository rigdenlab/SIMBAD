"""Class to store dat file info"""

__author__ = "Adam Simpkin"
__date__ = "16 April 2018"
__version__ = "0.1"

from simbad.core import ScoreBase


class DatModelScore(ScoreBase):
    """A dat model score storing class"""

    __slots__ = ("pdb_code", "dat_path", "mw_diff", "x", "y", "z", "intrad", "solvent", "nmol", "ermsd")

    def __init__(self, pdb_code, dat_path,  mw_diff, x, y, z, intrad, solvent, nmol, ermsd):
        self.pdb_code = pdb_code
        self.dat_path = dat_path
        self.mw_diff = mw_diff
        self.x = x
        self.y = y
        self.z = z
        self.intrad = intrad
        self.solvent = solvent
        self.nmol = nmol
        self.ermsd = ermsd

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code} dat_path={dat_path} mw_diff={mw_diff} x={x} y={y} z={z} intrad={intrad}" \
                 "solvent={solvent} nmol={nmol} ermsd={ermsd})"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`DatModelScore <simbad.score.dat_score.DatModelScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}