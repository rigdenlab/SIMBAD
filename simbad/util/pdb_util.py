import logging
import os
import string

import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import iotbx.pdb.fetch

import numpy as np

from simbad.chemistry import atomic_composition, periodic_table
from simbad.db import read_dat

logger = logging.getLogger(__name__)
three2one = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter


class PdbStructure(object):
    def __init__(self, pdbin):
        if pdbin.endswith(".dat"):
            pdb_str = read_dat(pdbin)
            self.pdb_input = iotbx.pdb.input(source_info=None, lines=pdb_str)
        else:
            self.pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
        self.hierarchy = self.pdb_input.construct_hierarchy()
        try:
            self.crystal_symmetry = self.pdb_input.crystal_symmetry()
        except AssertionError:
            logger.debug('Unable to generate crystal symmetry for %s', pdbin)
            self.crystal_symmetry = None

    @property
    def molecular_weight(self):
        mw = 0
        hydrogen_atoms = 0
        for m in self.hierarchy.models():
            for c in m.chains():
                for rg in c.residue_groups():
                    resseq = None
                    for ag in rg.atom_groups():
                        if ag.resname in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter \
                                and resseq != rg.resseq:
                            resseq = rg.resseq
                            try:
                                hydrogen_atoms += atomic_composition[ag.resname].H
                            except AttributeError:
                                logger.debug("Ignoring non-standard amino acid: %s", ag.resname)
                        for atom in ag.atoms():
                            if ag.resname.strip() == 'HOH' or ag.resname.strip() == 'WAT':
                                pass
                            else:
                                # Be careful, models might not have the last element column
                                if atom.element.strip():
                                    aname = atom.element.strip()
                                else:
                                    aname = atom.name.strip()
                                    aname = aname.translate(None, string.digits)[0]
                                try:
                                    mw += periodic_table[aname].atomic_mass * atom.occ
                                except AttributeError:
                                    try:
                                        aname = ''.join([i for i in aname if not i.isdigit()])
                                        mw += periodic_table[aname].atomic_mass * atom.occ
                                    except AttributeError:
                                        logger.debug("Ignoring non-standard atom type: %s", aname)

        mw += hydrogen_atoms * periodic_table['H'].atomic_mass
        return mw

    @property
    def integration_box(self):
        try:
            resolution = float(self.pdb_input.extract_remark_iii_records(2)[0].split()[3])
        except IndexError:
            resolution = 2.0
        chain = self.hierarchy.models()[0].chains()[0]
        xyz = np.zeros((chain.atoms_size(), 3))
        for i, atom in enumerate(chain.atoms()):
            xyz[i] = atom.xyz
        diffs = np.asarray([
            np.ptp(xyz[:, 0]),
            np.ptp(xyz[:, 1]),
            np.ptp(xyz[:, 2]),
        ])
        intrad = diffs.min() * 0.75
        x, y, z = diffs + intrad + resolution
        return x.item(), y.item(), z.item(), intrad.item()

    @property
    def nchains(self):
        return len(self.hierarchy.models()[0].chains())

    @property
    def nres(self):
        nres = 0
        for m in self.hierarchy.models():
            for c in m.chains():
                for rg in c.residue_groups():
                    resseq = None
                    for ag in rg.atom_groups():
                        if ag.resname in three2one and resseq != rg.resseq:
                            nres += 1
                            resseq = rg.resseq
        if nres == 0:
            # No standard amino acids in model, check for all e.g. DNA
            for m in self.hierarchy.models():
                for c in m.chains():
                    for _ in c.residue_groups():
                        nres += 1
        return nres

    def keep_first_chain_only(self):
        self.select_chain(0)

    def select_chain(self, chain_idx):
        for i, m in enumerate(self.hierarchy.models()):
            if i != 0:
                self.hierarchy.remove_model(m)
        m = self.hierarchy.models()[0]
        for i, c in enumerate(m.chains()):
            if i != chain_idx:
                m.remove_chain(c)

    def save(self, pdbout, remarks=[]):
        with open(pdbout, 'w') as f_out:
            for remark in remarks:
                f_out.write("REMARK %s" % remark + os.linesep)
            f_out.write(
                self.hierarchy.as_pdb_string(
                    anisou=False, write_scale_records=True, crystal_symmetry=self.crystal_symmetry))


def get_pdb_content(pdb_code):
    """Download iotbx data structure from pdb code

    Parameters
    ----------
    pdb_code : str
        4 letter pdb code

    Returns
    -------
    str
        Content of the downloaded file

    """
    import iotbx.pdb.fetch
    import urllib2

    # Work around until cctbx/cctbx_project#118 in release
    def cctbx_workaround(pdb_code):
        import ssl
        context = ssl._create_unverified_context()
        url_frame = "https://pdb-redo.eu/db/{0}/{0}_final.pdb"
        return urllib2.urlopen(url_frame.format(code), context=context)

    try:
        try:
            content = cctbx_workaround(pdb_code)
        except Exception:
            content = iotbx.pdb.fetch.fetch(pdb_code, data_type='pdb', format='pdb', mirror='pdbe')
        logger.debug("Downloaded PDB entry %s from %s", pdb_code, content.url)
        return content.read()
    except Exception as e:
        logger.critical(e)
        return None
