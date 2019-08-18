import gzip
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
    def __init__(self):
        self.pdb_input = None
        self.hierarchy = None
        self.crystal_symmetry = None

    @classmethod
    def from_file(cls, input_file):
        struct = cls()
        if input_file.endswith(".dat"):
            struct.pdb_input = iotbx.pdb.input(source_info=None, lines=read_dat(input_file))
        elif input_file.endswith(".pdb") or input_file.endswith(".ent"):
            struct.pdb_input = iotbx.pdb.pdb_input(file_name=input_file)
        elif input_file.endswith(".ent.gz"):
            with gzip.open(input_file, "rb") as f_in:
                pdb_str = f_in.read()
            struct.pdb_input = iotbx.pdb.input(source_info=None, lines=f_in.read())
        struct.hierarchy = struct.pdb_input.construct_hierarchy()
        struct.assert_hierarchy()
        struct.set_crystal_symmetry(input_file)
        return struct

    @classmethod
    def from_pdb_code(cls, pdb_code):
        struct = cls()
        content = struct.get_pdb_content(pdb_code)
        if content:
            struct.pdb_input = iotbx.pdb.input(source_info=None, lines=content)
            struct.hierarchy = struct.pdb_input.construct_hierarchy()
            struct.assert_hierarchy()
            struct.set_crystal_symmetry(pdb_code)
        else:
            raise RuntimeError
        return struct

    def set_crystal_symmetry(self, source):
        try:
            self.crystal_symmetry = self.pdb_input.crystal_symmetry()
        except AssertionError:
            logger.debug("Unable to generate crystal symmetry for %s", source)
            self.crystal_symmetry = None

    def assert_hierarchy(self):
        assert len(self.hierarchy.models()) > 0, "No models found in hierarchy"

    @staticmethod
    def get_pdb_content(pdb_code):
        import urllib2

        # Work around until cctbx/cctbx_project#118 in release
        def cctbx_workaround(pdb_code):
            import ssl

            context = ssl._create_unverified_context()
            url_frame = "https://pdb-redo.eu/db/{0}/{0}_final.pdb"
            return urllib2.urlopen(url_frame.format(pdb_code), context=context)

        try:
            try:
                content = cctbx_workaround(pdb_code)
            except Exception:
                content = iotbx.pdb.fetch.fetch(pdb_code, data_type="pdb", format="pdb", mirror="pdbe")
            logger.debug("Downloaded PDB entry %s from %s", pdb_code, content.url)
            return content.read()
        except Exception as e:
            logger.critical("Encountered problem downloading PDB %s: %s", pdb_code, e)
            return None

    @property
    def get_sequence_info(self):
        chain2data = {}
        unique_chains = []
        for c in set(self.hierarchy.models()[0].chains()):
            if not c.is_protein():
                continue
            got = False
            seq = ""
            for r in c.conformers()[0].residues():
                if any([not a.hetero for a in r.atoms()]):
                    if r.resname in three2one:
                        got = True
                        seq += three2one[r.resname]
            if got and seq not in unique_chains:
                chain2data[c.id] = seq
                unique_chains.append(seq)
        return chain2data

    @property
    def molecular_weight(self):
        mw = 0
        hydrogen_atoms = 0
        for c in set(self.hierarchy.models()[0].chains()):
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter and resseq != rg.resseq:
                        resseq = rg.resseq
                        try:
                            hydrogen_atoms += atomic_composition[ag.resname].H
                        except AttributeError:
                            logger.debug("Ignoring non-standard amino acid: %s", ag.resname)
                    for atom in ag.atoms():
                        if ag.resname.strip() == "HOH" or ag.resname.strip() == "WAT":
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
                                    aname = "".join([i for i in aname if not i.isdigit()])
                                    mw += periodic_table[aname].atomic_mass * atom.occ
                                except AttributeError:
                                    logger.debug("Ignoring non-standard atom type: %s", aname)

        mw += hydrogen_atoms * periodic_table["H"].atomic_mass
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
        diffs = np.asarray([np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]), np.ptp(xyz[:, 2])])
        intrad = diffs.min() * 0.75
        x, y, z = diffs + intrad + resolution
        return x.item(), y.item(), z.item(), intrad.item()

    @property
    def nchains(self):
        return len(self.hierarchy.models()[0].chains())

    @property
    def nres(self):
        nres = 0
        for c in set(self.hierarchy.models()[0].chains()):
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in three2one and resseq != rg.resseq:
                        nres += 1
                        resseq = rg.resseq
        if nres == 0:
            # No standard amino acids in model, check for all e.g. DNA
            for c in set(self.hierarchy.models()[0].chains()):
                for _ in c.residue_groups():
                    nres += 1
        return nres

    def keep_first_chain_only(self):
        self.select_chain_by_idx(0)

    def select_chain_by_idx(self, chain_idx):
        for i, m in enumerate(self.hierarchy.models()):
            if i != 0:
                self.hierarchy.remove_model(m)
        m = self.hierarchy.models()[0]
        for i, c in enumerate(m.chains()):
            if i != chain_idx:
                m.remove_chain(c)

    def select_chain_by_id(self, chain_id):
        for i, m in enumerate(self.hierarchy.models()):
            if i != 0:
                self.hierarchy.remove_model(m)
        m = self.hierarchy.models()[0]
        for c in m.chains():
            if c.id != chain_id or not c.is_protein():
                m.remove_chain(c)

    def standardize(self):
        for m in self.hierarchy.models():
            for c in m.chains():
                for rg in c.residue_groups():
                    for ag in rg.atom_groups():
                        for a in ag.atoms():
                            if a.element.strip().upper() == "H" or a.hetero:
                                ag.remove_atom(a)
                        if ag.atoms_size() == 0:
                            rg.remove_atom_group(ag)
                    if rg.atom_groups_size() == 0:
                        c.remove_residue_group(rg)
                if c.residue_groups_size() == 0:
                    m.remove_chain(c)
            if m.chains_size() == 0:
                self.hierarchy.remove_model(m)

    def save(self, pdbout, remarks=[]):
        with open(pdbout, "w") as f_out:
            for remark in remarks:
                f_out.write("REMARK %s" % remark + os.linesep)
            f_out.write(self.hierarchy.as_pdb_string(anisou=False, write_scale_records=True, crystal_symmetry=self.crystal_symmetry))
