"""Useful manipulations on PDB files"""

from __future__ import division, print_function

__author__ = "Adam Simpkin, Jens Thomas & Felix Simkovic"
__date__ = "21 Apr 2017"
__version__ = "2.0"

import os
import iotbx.pdb
import iotbx.pdb.amino_acid_codes

three2one = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter


def _cache(pdbin):
    """Cache the PDB input file"""
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    crystal_symmetry = pdb_input.crystal_symmetry()
    hierarchy = pdb_input.construct_hierarchy()
    return pdb_input, hierarchy, crystal_symmetry


def _first_chain_only(h):
    """Remove everything from hierarchy but the first chain"""
    for i, m in enumerate(h.models()):
        if i != 0:
            h.remove_model(m)
    m = h.models()[0]
    for i, c in enumerate(m.chains()):
        if i != 0:
            m.remove_chain(c)


def _number_of_chains(h):
    """Returns the number of chains in a hierarchy"""
    # Only check for the first model if multi-model input
    for i, m in enumerate(h.models()):
        if i != 0:
            h.remove_model(m)
    m = h.models()[0]
    return len(m.chains())


def _number_of_residues(h):
    """Return the number of residues in a hierarchy"""
    nres = 0
    for m in h.models():
        for c in m.chains():
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in three2one and resseq != rg.resseq:
                        nres += 1
                        resseq = rg.resseq
    return nres


def _save(pdbout, hierarchy, crystal_symmetry=None, remarks=[]):
    """Save the CCTBX hierarchy to a file"""
    with open(pdbout, 'w') as f_out:
        for remark in remarks:
            f_out.write("REMARK %s" % remark + os.linesep)
        f_out.write(hierarchy.as_pdb_string(
            anisou=False, write_scale_records=True, crystal_symmetry=crystal_symmetry
        ))


def _select_chain(h, chain_idx):
    """Select atoms from hierarchy by index"""
    for i, m in enumerate(h.models()):
        if i != 0:
            h.remove_model(m)
    m = h.models()[0]
    for i, c in enumerate(m.chains()):
        if i != chain_idx:
            m.remove_chain(c)
    return h


def to_single_chain(pdbin, pdbout):
    """Condense a single-model multi-chain pdb to a single-chain pdb

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB

    """
    _, hierarchy, symmetry = _cache(pdbin)
    _first_chain_only(hierarchy)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry)


def number_of_chains(pdbin):
    """Return the number of chains in a multi-chain pdb

    Parameters
    ----------
    pdbin : str
        The path to the input PDB

    Returns
    -------
    int
        The number of chains
    """
    _, hierarchy, _ = _cache(pdbin)
    nchains = _number_of_chains(hierarchy)
    return nchains


def number_of_residues(pdbin):
    """Return the number of residues in a multi-chain pdb

    Parameters
    ----------
    pdbin : str
        The path to the input PDB
    chain_idx : int
        Specify a specific chain by index

    Returns
    -------
    int
        The number of residues in the PDB
    """

    _, hierarchy, _ = _cache(pdbin)
    nres = _number_of_residues(hierarchy)
    return nres