"""Various miscellaneous functions"""

__author__ = "Adam Simpkin, Felix Simkovic & Jens Thomas"
__date__ = "05 May 2017"
__version__ = "1.0"

import logging
import os
import string

import iotbx.pdb
import iotbx.pdb.amino_acid_codes

import mbkit.chemistry
import mbkit.dispatch.cexectools

logger = logging.getLogger(__name__)


def filename_append(filename=None, astr=None, directory=None, separator="_"):
    """Append astr to filename, before the suffix, and return the new filename."""
    dirname, fname = os.path.split(filename)
    name, suffix = os.path.splitext(fname)
    name = name + separator + astr + suffix
    if directory is None:
        directory = dirname
    return os.path.join(directory, name)


def molecular_weight(pdbin):
    """Function to calculate the molecular weight of a model

    Parameters
    ----------
    pdbin : str
       Path to input pdb

    Returns
    -------
    float
       Molecular weight of input model
    """

    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    # Define storage variables
    mw = 0
    hydrogen_atoms = 0

    # Collect all the data using the hierarchy
    for m in hierarchy.models():
        for c in m.chains():
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter and resseq != rg.resseq:
                        resseq = rg.resseq
                        try:
                            hydrogen_atoms += mbkit.chemistry.atomic_composition[ag.resname].H
                        except AttributeError:
                            logger.debug('Ignoring non-standard amino acid: %s', ag.resname)
                    for atom in ag.atoms():
                        if ag.resname.strip() == 'HOH' or ag.resname.strip() == 'WAT':
                            # Ignore water atoms
                            pass
                        else:
                            # Be careful, models might not have the last element column
                            if atom.element.strip():
                                aname = atom.element.strip()
                            else:
                                aname = atom.name.strip()
                                aname = aname.translate(None, string.digits)[0]
                            try:
                                mw += mbkit.chemistry.periodic_table[aname].atomic_mass * atom.occ
                            except AttributeError:
                                try:
                                    aname = ''.join([i for i in aname if not i.isdigit()])
                                    mw += mbkit.chemistry.periodic_table[aname].atomic_mass * atom.occ
                                except AttributeError:
                                    logger.debug('Ignoring non-standard atom type: %s', aname)

    mw += hydrogen_atoms * mbkit.chemistry.periodic_table['H'].atomic_mass

    return mw

