"""Module to skim the PDB for similar unit cells"""

from __future__ import division

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "30 Jun 2017"
__version__ = "1.0"

import ast
import datetime
import glob
import logging
import numpy as np
import os

from simbad.core.lattice_score import LatticeSearchResult
import simbad.util.pdb_util

logger = logging.getLogger(__name__)


class LatticeSearch(object):
    """A class to do a search for PDB entries with similar unit cell dimensions."""

    class _ResultCache(object):
        def __init__(self):
            self._codes = set()
            self._data = []

        def __iter__(self):
            for e in self._data:
                yield e

        def append(self, score):
            code = score.pdb_code.upper()
            if code in self._codes:
                return
            self._codes.add(code)
            self._data.append(score)

    def __init__(self, lattice_db_fname, model_dir):
        """Initialize a new Lattice Search class

        Parameters
        ----------
        lattice_db_fname : str
           The path to the lattice database [pickle format]
        max_to_keep : int, optional
           The maximum number of results to keep [default: 20]

        """
        self._lattice_db_fname = None

        self.lattice_db_fname = lattice_db_fname
        self.model_dir = model_dir

        self.results = []

    @property
    def lattice_db_fname(self):
        """The path to the lattice database"""
        return self._lattice_db_fname

    @lattice_db_fname.setter
    def lattice_db_fname(self, lattice_db_fname):
        """Define the lattice database filename"""
        # Check how old the database is, if older than 90 days print
        # message to suggest updating
        timestamp = os.path.getmtime(lattice_db_fname)
        if (datetime.date.today() - datetime.date.fromtimestamp(timestamp)).days > 90:
            logger.info("Lattice database is older than 90 days, consider updating!\n" 'Use the command "simbad-database lattice" in your Terminal')
        self._lattice_db_fname = lattice_db_fname

    def search(self, space_group, unit_cell, tolerance=0.05, max_to_keep=50, max_penalty=12):
        """Search for similar Niggli cells

        Parameters
        ----------
        unit_cell : list, tuple
           The parameters of the unit cell
        space_group : str
           The space group
        tolerance : int, float, optional
           The tolerance applied for Niggli cell comparison [default: 0.05]
        max_to_keep : int, optional
           The top-N number of results to return [default: 50]
        max_penalty : int, optional
           The total penalty score over which results are ignored [default: 12]

        """
        space_group = self.check_sg(space_group)
        niggli_cell = self.calculate_niggli_cell(unit_cell, space_group)
        niggli_cell = np.array(niggli_cell)

        tol_niggli_cell = niggli_cell * tolerance

        results = self._ResultCache()
        with np.load(self.lattice_db_fname) as compressed:
            for entry in compressed["arr_0"]:
                pdb_code = "".join(chr(c) for c in entry[:4].astype("uint8"))
                pdb_path = os.path.join(self.model_dir, "{}.pdb".format(pdb_code))
                alt_cell = chr(int(entry[4])) if entry[4] != 0.0 else " "
                db_cell = entry[5:]

                if self.cell_within_tolerance(niggli_cell, db_cell, tol_niggli_cell):
                    total_pen, length_pen, angle_pen = self.calculate_penalties(niggli_cell, db_cell)
                    vol_diff = self.calculate_volume_difference(niggli_cell, db_cell)
                    if total_pen < max_penalty:
                        prob = self.calculate_probability(total_pen)
                        score = LatticeSearchResult(pdb_code, pdb_path, alt_cell, db_cell, vol_diff, total_pen, length_pen, angle_pen, prob)
                        results.append(score)

        results_sorted = sorted(results, key=lambda x: float(x.total_penalty), reverse=False)
        self.results = results_sorted[:max_to_keep]

    @classmethod
    def calculate_penalties(cls, query, reference):
        """Calculate the linear cell variation between unit cells

        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters

        Returns
        -------
        float
           Total penalty
        float
           Length penalty
        float
           Angle penalty
        """
        delta = np.absolute(query - reference)
        length_penalty = delta[:3].sum().item()
        angle_penalty = delta[3:].sum().item()
        total_penalty = length_penalty + angle_penalty
        return total_penalty, length_penalty, angle_penalty

    @classmethod
    def calculate_probability(cls, penalty_score):
        """Calculate the probability that a penalty score will give a solution

        Parameters
        ----------
        penalty_score : float
            The total penalty score calculate for a search result

        Returns
        -------
        float
            Probability score
        """
        x = -1.01 * penalty_score + 2.11
        return np.around(1 / (1 + np.exp(-x)), decimals=3)

    @classmethod
    def cell_within_tolerance(cls, query, reference, tolerance):
        """Compare two cells and determine if ``query`` is within ``reference`` cell parameter tolerance

        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters
        tolerance : list, tuple
           The tolerance cell parameter values

        Returns
        -------
        bool

        """
        return np.all(np.absolute(query - reference) <= tolerance)

    @classmethod
    def calculate_volume_difference(cls, query, reference):
        """Calculate the difference in volume between the query unit cell and the reference unit cell

        Parameters
        ----------
        query : list, tuple
           The query cell parameters
        reference : list, tuple
           The reference cell parameters

        Returns
        -------
        float
            The absolute difference in cell volumes
        """
        import cctbx.uctbx
        cell_volume_1 = cctbx.uctbx.unit_cell(list(query)).volume()
        cell_volume_2 = cctbx.uctbx.unit_cell(list(reference)).volume()
        return np.absolute(cell_volume_1 - cell_volume_2).round(decimals=3).item()

    @staticmethod
    def calculate_niggli_cell(unit_cell, space_group):
        """Calculate the parameters of the Niggli cell

        Parameters
        ----------
        unit_cell : list, tuple
           The parameters of the unit cell
        space_group : str
           The space group

        Returns
        -------
        list
           The Niggli cell parameters

        """
        import cctbx.crystal
        import cctbx.uctbx

        unit_cell = cctbx.uctbx.unit_cell(list(unit_cell))
        xs = cctbx.crystal.symmetry(unit_cell=unit_cell, space_group=space_group, correct_rhombohedral_setting_if_necessary=True)
        niggli_cell = xs.change_basis(xs.change_of_basis_op_to_niggli_cell()).unit_cell()
        niggli_cell = list(ast.literal_eval(str(niggli_cell)))
        logger.info("Niggli cell calculated as: [%s]", ", ".join(map(str, niggli_cell)))
        return niggli_cell

    @staticmethod
    def check_sg(sg):
        """Check the space group for known anomalies"""
        sg_conversion = {
            "A1": "P1",
            "B2": "B112",
            "C1211": "C2",
            "F422": "I422",
            "I21": "I2",
            "I1211": "I2",
            "P21212A": "P212121",
            "R3": "R3:R",
            "C4212": "P422",
        }
        return sg_conversion.get(sg, sg)

    def copy_results(self, source, destination):
        """Copy the results from a local copy of the PDB

        Parameters
        ----------
        source : str
           The path to copy results from
        destination : str
           The path to save results to 

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        IOError
            Search result not found in installed PDB

        """
        if not self.results:
            raise ValueError("No search results found/available")

        if not os.path.isdir(destination):
            msg = "Output directory does not exist: {0}".format(destination)
            raise ValueError(msg)

        prefix, ext = os.path.basename(glob.glob(os.path.join(source, "*", "*"))[0]).split(".", 1)
        prefix = prefix[:-4]

        to_del = []
        for count, result in enumerate(self.results):
            f_name = os.path.join(source, "{0}", "{1}{2}.{3}").format(result.pdb_code[1:3].lower(), prefix, result.pdb_code.lower(), ext)
            f_name_out = os.path.join(destination, "{0}.pdb".format(result.pdb_code))
            try:
                struct = simbad.util.pdb_util.PdbStructure.from_file(f_name)
                struct.standardize()
                struct.save(f_name_out)
            except IOError:
                logger.warning("Encountered problem copying PDB %s from %s - removing entry from list", result.pdb_code, source)
                to_del.append(count)

        for i in reversed(to_del):
            self.results.pop(i)

    def download_results(self, destination):
        """Download the results directly from the PDB

        Parameters
        ----------
        destination : str
           The path to save results to 

        Raises
        ------
        ValueError
           No search results found/available
        ValueError
           Output directory does not exist
        RuntimeError
            Unable to download PDB

        """
        if not self.results:
            raise ValueError("No search results found/available")

        if not os.path.isdir(destination):
            msg = "Output directory does not exist: {0}".format(destination)
            raise ValueError(msg)

        to_del = []
        for count, result in enumerate(self.results):
            try:
                f_name_out = os.path.join(destination, result.pdb_code + ".pdb")
                struct = simbad.util.pdb_util.PdbStructure.from_pdb_code(result.pdb_code)
                struct.standardize()
                struct.save(f_name_out)
            except RuntimeError:
                logger.warning("Encountered problem downloading PDB %s - removing entry from list", result.pdb_code)
                to_del.append(count)

        for i in reversed(to_del):
            self.results.pop(i)

    def summarize(self, csvfile):
        """Summarize the search results

        Parameters
        ----------
        csvfile : str
           The path to an output CSV file

        """
        from simbad.util import summarize_result

        columns = [
            "alt",
            "a",
            "b",
            "c",
            "alpha",
            "beta",
            "gamma",
            "length_penalty",
            "angle_penalty",
            "total_penalty",
            "volume_difference",
            "probability_score",
        ]
        summarize_result(self.results, csv_file=csvfile, columns=columns)
