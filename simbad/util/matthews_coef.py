"""Module for running matthews coefficients to calculate solvent content"""

__author__ = "Adam Simpkin"
__date__ = "07 Oct 2017"
__version__ = "0.1"

from cctbx.crystal import symmetry
from mmtbx.scaling.matthews import density_calculator
from mmtbx.scaling.matthews import matthews_rupp
from simbad.util import molecular_weight


class SolventContent(object):
    def __init__(self, cell, sg):
        self.crystal_symmetry = symmetry(unit_cell=cell, space_group_symbol=sg)
        self.dens_calc = density_calculator(self.crystal_symmetry)

    def calculate(self, pdb):
        return self.dens_calc.solvent_fraction(molecular_weight(pdb), 0.74) * 100


class MatthewsCoefficient(object):
    def __init__(self, cell, sg):
        self.crystal_symmetry = symmetry(unit_cell=cell, space_group_symbol=sg)

    def calculate_content_ncopies(self, nres):
        result = matthews_rupp(self.crystal_symmetry, n_residues=nres)
        return result.solvent_content, result.n_copies


def solvent_content(pdbin, cell_parameters, space_group):
    """Get the solvent content for an input pdb

    Parameters
    ----------
    pdbin : str
        Path to input PDB file
    cell_parameters : str
        The parameters describing the unit cell of the crystal
    space_group : str
        The space group of the crystal

    Returns
    -------
    float
        The solvent content
    """
    return SolventContent(cell_parameters, space_group).calculate(pdbin)


def matthews_coef(cell_parameters, space_group, nres):
    """Function to run matthews coefficient to decide if the model can fit in the unit cell

    Parameters
    ----------
    cell_parameters
        The parameters of the unit cell
    space_group
        The space group of the crystal
    nres
        The number of residues in input model

    Returns
    -------
    float
        predicted solvent content
    int
        predicted number of copies of protein copies
    """
    mathcoeff = MatthewsCoefficient(cell_parameters, space_group)
    return mathcoeff.calculate_content_ncopies(nres)
