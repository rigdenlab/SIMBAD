"""Class to store Lattice scores"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "05 Mar 2017"
__version__ = "0.1"

from collections import namedtuple, OrderedDict

_fields = [
    "pdb_code",
    "pdb_path",
    "alt",
    "unit_cell",
    "volume_difference",
    "total_penalty",
    "length_penalty",
    "angle_penalty",
    "probability_score",
]


class LatticeSearchResult(namedtuple("LatticeSearchResult", _fields)):
    def _asdict(self):
        dictionary = OrderedDict()
        for k in _fields:
            if k == "unit_cell":
                for k, v in zip(["a", "b", "c", "alpha", "beta", "gamma"], self.unit_cell):
                    dictionary[k] = v
            else:
                dictionary[k] = getattr(self, k)
        return dictionary
