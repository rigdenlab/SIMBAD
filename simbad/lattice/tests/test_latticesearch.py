"""Test functions for simbad.lattice.latticesearch.LatticeSearch"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import numpy as np
import unittest
from simbad.lattice.lattice_search import LatticeSearch

try:
    SHARE_DIR = os.environ['SIMBAD_ROOT']
except KeyError:
    from simbad.command_line import CCP4RootDirectory
    SHARE_DIR = os.path.join(str(CCP4RootDirectory()), "share", "simbad")


class Test(unittest.TestCase):
    """Unit test"""

    @classmethod
    def setUpClass(cls):
        lattice_db = os.path.join(SHARE_DIR, "static", "niggli_database.npz")
        cls.LS = LatticeSearch(lattice_db, os.getcwd())

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_search_1(self):
        """Test case for LatticeSearch.search"""

        # Process the data from the toxd test case
        space_group = "P212121"
        unit_cell = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]

        self.LS.search(space_group, unit_cell)
        results = self.LS.results

        # Take the name of the top result (should be toxd)
        for i, r in enumerate(results):
            if i == 0:
                data = r.pdb_code

        reference_data = "1DTX"

        self.assertEqual(data, reference_data)

    def test_calculate_penalties_1(self):
        """Test case for LatticeSearch.calculate_penalties"""

        # Same cells
        unit_cell_1 = np.array([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.array([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])

        data = self.LS.calculate_penalties(unit_cell_1, unit_cell_2)
        reference_data = (0.0, 0.0, 0.0)

        self.assertEqual(data, reference_data)

    def test_calculate_penalties_2(self):
        """Test case for LatticeSearch.calculate_penalties"""

        # Different cells
        unit_cell_1 = np.array([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.array([41.34, 123.01, 93.23, 120.00, 90.00, 89.00])

        data = self.LS.calculate_penalties(unit_cell_1, unit_cell_2)
        reference_data = (217.56, 186.56, 31.0)

        self.assertEqual(data, reference_data)

    def test_calculate_probability_1(self):
        """Test case for LatticeSearch.calculate_probability"""

        score = 0.0
        data = self.LS.calculate_probability(score)
        reference_data = 0.892

        self.assertEqual(data, reference_data)

    def test_calculate_probability_2(self):
        """Test case for LatticeSearch.calculate_probability"""

        score = 0.25
        data = self.LS.calculate_probability(score)
        reference_data = 0.865

        self.assertEqual(data, reference_data)

    def test_cell_within_tolerance_1(self):
        """Test case for LatticeSearch.cell_within_tolerance"""

        # Same cells
        unit_cell_1 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        tolerance = unit_cell_1 * 0.05

        data = self.LS.cell_within_tolerance(unit_cell_1, unit_cell_2, tolerance)
        reference_data = True

        self.assertEqual(data, reference_data)

    def test_cell_within_tolerance_2(self):
        """Test case for LatticeSearch.cell_within_tolerance"""

        # One parameter beyond 0.05 tolerance
        unit_cell_1 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.asarray([69.16, 38.73, 23.19, 90.00, 90.00, 90.00])
        tolerance = unit_cell_1 * 0.05

        data = self.LS.cell_within_tolerance(unit_cell_1, unit_cell_2, tolerance)
        reference_data = False

        self.assertEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_calculate_volume_difference_1(self):
        """Test case for LatticeSearch.calculate_volume_difference"""

        # Same cells
        unit_cell_1 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])

        data = self.LS.calculate_volume_difference(unit_cell_1, unit_cell_2)
        reference_data = 0.00

        self.assertEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_calculate_volume_difference_2(self):
        """Test case for LatticeSearch.calculate_volume_difference"""

        # Same cells
        unit_cell_1 = np.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = np.asarray([63.28, 38.73, 29.01, 90.00, 90.00, 90.00])

        data = self.LS.calculate_volume_difference(unit_cell_1, unit_cell_2)
        reference_data = 5012.925

        self.assertEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_calculate_niggli_cell_1(self):
        """Test case for LatticeSearch.calculate_niggli_cell"""

        space_group = "P212121"
        unit_cell = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]

        data = self.LS.calculate_niggli_cell(unit_cell, space_group)
        reference_data = [23.19, 38.73, 73.58, 90.0, 90.0, 90.0]

        self.assertEqual(data, reference_data)

    def test_check_sg_1(self):
        """Test case for LatticeSearch.check_sg"""

        sg = "A1"
        data = self.LS.check_sg(sg)
        reference_data = "P1"

        self.assertEqual(data, reference_data)


if __name__ == "__main__":
    unittest.main()
