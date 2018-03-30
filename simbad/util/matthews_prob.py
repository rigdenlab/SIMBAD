"""Module for running matthews probabilities to calculate solvent content"""

__author__ = "Adam Simpkin"
__date__ = "10 Oct 2017"
__version__ = "0.2"

import contextlib
import io
import sys

from cctbx.crystal import symmetry
from mmtbx.scaling.matthews import density_calculator
from mmtbx.scaling.matthews import matthews_rupp
from simbad.util.pdb_util import PdbStructure


@contextlib.contextmanager
def no_stdout():
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout


class SolventContent(object):
    def __init__(self, cell, sg):
        self.crystal_symmetry = symmetry(unit_cell=cell, space_group_symbol=sg)
        self.dens_calc = density_calculator(self.crystal_symmetry)

    def calculate_from_file(self, pdb):
        return self.calculate_from_struct(PdbStructure(pdb))

    def calculate_from_struct(self, struct):
        return self._calculate(struct.molecular_weight)

    def _calculate(self, mw):
        return self.dens_calc.solvent_fraction(mw, 0.74) * 100


class MatthewsProbability(object):
    def __init__(self, cell, sg):
        self.crystal_symmetry = symmetry(unit_cell=cell, space_group_symbol=sg)

    def calculate_content_ncopies_from_file(self, pdb):
        return self.calculate_content_ncopies_from_struct(PdbStructure(pdb))

    def calculate_content_ncopies_from_struct(self, struct):
        return self._calculate(struct.nres)

    def _calculate(self, nres):
        with no_stdout():
            result = matthews_rupp(self.crystal_symmetry, n_residues=nres)
            return result.solvent_content, result.n_copies
