import gemmi
import logging
import numpy as np
import os

from simbad.db import read_dat

logger = logging.getLogger(__name__)


class PdbStructure(object):
    def __init__(self):
        self.structure = None

    @classmethod
    def from_file(cls, input_file):
        struct = cls()
        if input_file.endswith(".dat"):
            struct.structure = gemmi.read_pdb_string(read_dat(input_file))
        elif input_file.endswith(".pdb") or input_file.endswith(".ent"):
            struct.structure = gemmi.read_structure(input_file)
        elif input_file.endswith(".ent.gz"):
            struct.structure = gemmi.read_structure(input_file)
        return struct

    @classmethod
    def from_pdb_code(cls, pdb_code):
        struct = cls()
        content = struct.get_pdb_content(pdb_code)
        if content:
            struct.structure = gemmi.read_pdb_string(content)
        else:
            raise RuntimeError
        return struct

    @staticmethod
    def get_pdb_content(pdb_code):
        import iotbx.pdb.fetch
        try:
            try:
                content = iotbx.pdb.fetch.fetch(pdb_code, data_type="pdb", format="pdb", mirror="pdb-redo")
            except urllib2.HTTPError:
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
        for chain in self.structure[0]: # Just the first model
            id = chain.name
            seq = ""
            for residue in chain.first_conformer(): # Only use the main conformer
                res_info = gemmi.find_tabulated_residue(residue.name)
                if res_info.is_amino_acid(): # Only amino acids
                    seq += res_info.one_letter_code
            if seq not in unique_chains:
                chain2data[id] = seq
                unique_chains.append(seq)
        return chain2data

    @property
    def molecular_weight(self):
        mw = 0
        hydrogen_atoms = 2 # For each end
        for chain in self.structure[0]:  # Just first model
            for residue in chain:
                res_info = gemmi.find_tabulated_residue(residue.name)
                # Many PDB files don't include hydrogens so account for them here
                hydrogen_atoms += res_info.hydrogen_count - 2 # Due to peptide bonds
                if res_info.is_standard(): # Only count standard amino acids
                    for atom in residue:
                        if atom.is_hydrogen: # Ignore as hydrogens already counted
                            pass
                        mw += atom.element.weight * atom.occ
        mw += gemmi.Element('H').weight * hydrogen_atoms
        return mw

    @property
    def integration_box(self):
        resolution = self.structure.resolution if self.structure.resolution > 0 else 2.0
        chain = self.structure[0][0]
        xyz = np.zeros((chain.count_atom_sites(), 3))
        count = 0
        for residue in chain:
            for atom in residue:
                xyz[count] = (atom.pos.x, atom.pos.y, atom.pos.z)
                count += 1
        diffs = np.asarray([np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]), np.ptp(xyz[:, 2])])
        intrad = diffs.min() * 0.75
        x, y, z = diffs + intrad + resolution
        return x.item(), y.item(), z.item(), intrad.item()

    @property
    def nchains(self):
        return len(self.structure[0])

    @property
    def nres(self):
        nres = 0
        for chain in self.structure[0]:
            for residue in chain:
                res_info = gemmi.find_tabulated_residue(residue.name)
                if res_info.is_standard():
                    nres += 1
        if nres == 0:
            # No standard amino acids in model, check for all e.g. DNA
            for chain in self.structure[0]:
                for _ in chain:
                    nres += 1
        return nres

    def keep_first_chain_only(self):
        self.select_chain_by_idx(0)


    def keep_first_model_only(self):
        models = [m.name for m in self.structure]
        for model in models[1:]:
            del self.structure[model]

    def remove_empty_models(self):
        for model in self.structure:
            if len(model.subchains()) == 0:
                del self.structure[model.name]

    def select_chain_by_idx(self, chain_idx):
        self.keep_first_model_only()
        model = self.structure[0]
        for i, chain in enumerate(model):
            if i != chain_idx:
                model.remove_chain(chain.name)

    def select_chain_by_id(self, chain_id):
        self.keep_first_model_only()
        model = self.structure[0]
        for chain in model:
            if chain.name != chain_id:
                model.remove_chain(chain.name)

    def standardize(self):
        self.structure.remove_hydrogens()
        self.structure.remove_ligands_and_waters()
        self.structure.remove_empty_chains()
        self.remove_empty_models()

    def save(self, pdbout, remarks=[]):
        self.structure.remove_ligands_and_waters()
        pdb_string = [line for line in self.structure.make_minimal_pdb().split('\n') if not line.startswith('ANISOU')]
        with open(pdbout, "w") as f_out:
            for remark in remarks:
                f_out.write("REMARK %s" % remark + os.linesep)
            for line in pdb_string:
                f_out.write(line + os.linesep)










