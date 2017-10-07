"""Module for running matthews coefficients to calculate solvent content"""

__author__ = "Adam Simpkin"
__date__ = "07 Oct 2017"
__version__ = "0.1"

from simbad.util import molecular_weight
import cctbx.crystal
import mmtbx.scaling.matthews


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
    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=cell_parameters, space_group_symbol=space_group)
    dens_calc = mmtbx.scaling.matthews.density_calculator(crystal_symmetry)
    return dens_calc.solvent_fraction(molecular_weight(pdbin), 0.74) * 100


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

    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=cell_parameters, space_group_symbol=space_group)
    result = mmtbx.scaling.matthews.matthews_rupp(
        crystal_symmetry, n_residues=nres)
    return result.solvent_content, result.n_copies
