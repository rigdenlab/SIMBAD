#!/usr/bin/env ccp4-python
# encoding: utf-8

import logging
import os
import pandas as pd

LATTICE_ID = 'latt'
CONTAMINANT_ID = 'cont'
MORDA_ID = 'morda'

ID2STR = {
    LATTICE_ID: 'lattice',
    CONTAMINANT_ID: 'contaminant',
    MORDA_ID: 'MORDA',
}

logger = logging.getLogger(__name__)


class FileCollection(object):
    """Class for holding results files and annotations"""

    def __init__(self):
        self.mr_log = None
        self.ref_pdb = None
        self.ref_pdb_annotation = None
        self.ref_mtz_annotation = None
        self.ref_log = None
        self.ref_map = None
        self.diff_map = None


class SimbadResults(object):
    """Class to handle the results from a SIMBAD search

    Attributes
    ----------
    work_dir : str
        The SIMBAD work directory where SIMBAD was run
    lattice_results : :obj:`~pandas.DataFrame`
        :obj:`~pandas.DataFrame` containing the results of the lattice search
    lattice_mr_results : :obj:`~pandas.DataFrame`
        :obj:`~pandas.DataFrame` containing the results of the MR

    """

    def __init__(self, work_dir):
        assert os.path.isdir(work_dir), "Cannot find work_dir: %s" % work_dir
        self.work_dir = os.path.abspath(work_dir)

        # Create dataframes from results
        lattice_results = os.path.join(self.work_dir, LATTICE_ID, 'lattice_search.csv')
        self.lattice_results = pd.read_csv(lattice_results) if os.path.isfile(lattice_results) else None
        lattice_mr_results = os.path.join(self.work_dir, LATTICE_ID, 'lattice_mr.csv')
        self.lattice_mr_results = pd.read_csv(lattice_mr_results) if os.path.isfile(lattice_mr_results) else None

        contaminant_results = os.path.join(self.work_dir, CONTAMINANT_ID, 'rot_search.csv')
        self.contaminant_results = pd.read_csv(contaminant_results) if os.path.isfile(contaminant_results) else None
        contaminant_mr_results = os.path.join(self.work_dir, CONTAMINANT_ID, 'cont_mr.csv')
        self.contaminant_mr_results = pd.read_csv(contaminant_mr_results) if os.path.isfile(
            contaminant_mr_results) else None

        morda_db_results = os.path.join(self.work_dir, MORDA_ID, 'rot_search.csv')
        self.morda_db_results = pd.read_csv(morda_db_results) if os.path.isfile(morda_db_results) else None
        morda_db_mr_results = os.path.join(self.work_dir, MORDA_ID, 'morda_mr.csv')
        self.morda_db_mr_results = pd.read_csv(morda_db_mr_results) if os.path.isfile(morda_db_mr_results) else None

    def top_files(self, num_results=3):
        """Return a list of num_results FileCollection objects with the results of the best MR steps"""
        if self.morda_db_mr_results is not None:
            return self.get_files(MORDA_ID, num_results=num_results)
        elif self.contaminant_mr_results is not None:
            return self.get_files(CONTAMINANT_ID, num_results=num_results)
        elif self.lattice_mr_results is not None:
            return self.get_files(LATTICE_ID, num_results=num_results)
        else:
            None

    def get_files(self, search_type, num_results=10):
        """Return the best num_results results files from search_type

        Parameters
        ----------
        search_type : str
            The type of search undertaken

        Returns
        -------
        """

        if search_type == LATTICE_ID:
            df = self.lattice_mr_results
        elif search_type == CONTAMINANT_ID:
            df = self.contaminant_mr_results
        elif search_type == MORDA_ID:
            df = self.morda_db_mr_results
        else:
            raise RuntimeError("Unrecognised search type: %s" % search_type)

        results = []
        for i in range(0, num_results):
            try:
                fc = FileCollection()
                pdb_code = df.loc[i][0]
                mr_program = list(df)[1][0:6]
                mr_workdir = os.path.join(self.work_dir, search_type, 'mr_search', pdb_code, 'mr', mr_program)
                fc.mr_log = os.path.join(mr_workdir, '{0}_mr.log'.format(pdb_code))
                fc.ref_pdb = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.pdb'.format(pdb_code))
                fc.ref_pdb_annotation = 'PDB #{0} from REFMAC-refined result of the {1} search'.format(
                    i + 1, ID2STR[search_type])
                fc.ref_mtz = os.path.join(mr_workdir, 'refine', '{0}_refinement_output.mtz'.format(pdb_code))
                fc.ref_mtz_annotation = 'MTZ #{0} from REFMAC-refined result of the {1} search'.format(
                    i + 1, ID2STR[search_type])
                fc.ref_log = os.path.join(mr_workdir, 'refine', '{0}_ref.log'.format(pdb_code))
                fc.ref_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_2fofcwt.map'.format(pdb_code))
                fc.diff_map = os.path.join(mr_workdir, 'refine', '{0}_refmac_fofcwt.map'.format(pdb_code))
                results.append(fc)
            except KeyError:
                logger.debug("No result found at position %s", (i + 1))

        return results


if __name__ == "__main__":
    # Setup argument parser
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(dest="work_dir")

    # Process arguments
    args = parser.parse_args()

    SR = SimbadResults(args.work_dir)

    for fc in SR.top_files():
        logger.debug("FC %s %s %s", fc, fc.ref_log, fc.ref_pdb_annotation)
