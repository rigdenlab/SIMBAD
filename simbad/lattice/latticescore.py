"""Class to store Lattice scores"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "05 Mar 2017"
__version__ = "0.1"


class LatticeSearchResult(object):
    """A basic lattice parameter scoring class"""

    __slots__ = ('pdb_code', 'alt', 'unit_cell', 'total_penalty', 'length_penalty', 'angle_penalty')

    def __init__(self, pdb_code, alt, unit_cell, total_penalty, length_penalty, angle_penalty):
        self.pdb_code = pdb_code
        self.alt = alt
        self.unit_cell = unit_cell
        self.total_penalty = total_penalty
        self.length_penalty = length_penalty
        self.angle_penalty = angle_penalty

    def __repr__(self):
        template = "{name}(pdb_code={pdb_code} alt={alt} unit_cell={unit_cell} total_penalty={total_penalty} " \
                   + "length_penalty={length_penalty} angle_penalty={angle_penalty}"
        return template.format(self.__class__.__name__, **{k: getattr(self, k) for k in self.__class__.__slots__})

    def _as_dict(self):
        """Convert the :obj:`_LatticeParameterScore <simbad.lattice.search._LatticeParameterScore>`
        object to a dictionary"""
        dictionary = {}
        for k in self.__slots__:
            if k == 'unit_cell':
                for k, v in zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], self.unit_cell):
                    dictionary[k] = v
            else:
                dictionary[k] = getattr(self, k)
        return dictionary

