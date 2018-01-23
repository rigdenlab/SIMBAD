"""Test functions for simbad.lattice.latticesearch.LatticeSearch"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import numpy
import unittest
import simbad
from simbad.lattice.latticesearch import LatticeSearch


class Test(unittest.TestCase):
    """Unit test"""
    
    @classmethod
    def setUpClass(cls):
        lattice_db = simbad.LATTICE_DB
        cls.LS = LatticeSearch(lattice_db, os.getcwd())
    
    def test_search_1(self):
        """Test case for LatticeSearch.search"""
        
        # Process the data from the toxd test case
        space_group = 'P212121'
        unit_cell = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]

        self.LS.search(space_group, unit_cell)
        results = self.LS.results
        
        # Take the name of the top result (should be toxd)
        for i, r in enumerate(results):
            if i == 0:
                data = r.pdb_code
        
        reference_data = "1DTX"
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_penalty_1(self):
        """Test case for LatticeSearch.calculate_penalty"""
        
        # Same cells
        unit_cell_1 = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]
        unit_cell_2 = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]
        
        data = self.LS.calculate_penalty(unit_cell_1, unit_cell_2)
        reference_data = (0.0, 0.0, 0.0) 
        
        self.assertEqual(data, reference_data)
    
    def test_calculate_penalty_2(self):
        """Test case for LatticeSearch.calculate_penalty"""
        
        # Different cells
        unit_cell_1 = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]
        unit_cell_2 = [41.34, 123.01, 93.23, 120.00, 90.00, 89.00]
        
        data = self.LS.calculate_penalty(unit_cell_1, unit_cell_2)
        reference_data = (217.56, 186.56, 31.0)
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_probability_1(self):
        """Test case for LatticeSearch.calculate_probability"""
        
        score = 0.0
        data = self.LS.calculate_probability(score)
        reference_data = 1
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_probability_2(self):
        """Test case for LatticeSearch.calculate_probability"""
        
        score = 0.25
        data = self.LS.calculate_probability(score)
        reference_data = 0.902
        
        self.assertEqual(data, reference_data)
        
    def test_cell_within_tolerance_1(self):
        """Test case for LatticeSearch.cell_within_tolerance"""
        
        # Same cells
        unit_cell_1 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        tolerance = unit_cell_1 * 0.05
        
        data = self.LS.cell_within_tolerance(unit_cell_1, unit_cell_2, tolerance)
        reference_data = True
        
        self.assertEqual(data, reference_data)
        
    def test_cell_within_tolerance_2(self):
        """Test case for LatticeSearch.cell_within_tolerance"""
        
        # One parameter beyond 0.05 tolerance
        unit_cell_1 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = numpy.asarray([69.16, 38.73, 23.19, 90.00, 90.00, 90.00])
        tolerance = unit_cell_1 * 0.05
        
        data = self.LS.cell_within_tolerance(unit_cell_1, unit_cell_2, tolerance)
        reference_data = False
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_volume_difference_1(self):
        """Test case for LatticeSearch.calculate_volume_difference"""
        
        # Same cells
        unit_cell_1 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        
        data = self.LS.calculate_volume_difference(unit_cell_1, unit_cell_2)
        reference_data = 0.00 
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_volume_difference_2(self):
        """Test case for LatticeSearch.calculate_volume_difference"""
        
        # Same cells
        unit_cell_1 = numpy.asarray([73.58, 38.73, 23.19, 90.00, 90.00, 90.00])
        unit_cell_2 = numpy.asarray([63.28, 38.73, 29.01, 90.00, 90.00, 90.00])
        
        data = self.LS.calculate_volume_difference(unit_cell_1, unit_cell_2)
        reference_data = 5012.925
        
        self.assertEqual(data, reference_data)
        
    def test_calculate_niggli_cell_1(self):
        """Test case for LatticeSearch.calculate_niggli_cell"""
        
        space_group = 'P212121'
        unit_cell = [73.58, 38.73, 23.19, 90.00, 90.00, 90.00]
        
        data = self.LS.calculate_niggli_cell(unit_cell, space_group)
        reference_data = [23.19, 38.73, 73.58, 90.0, 90.0, 90.0]
        
        self.assertEqual(data, reference_data)
        
    def test_check_sg_1(self):
        """Test case for LatticeSearch.check_sg"""
        
        sg = 'A1'
        data = self.LS.check_sg(sg)
        reference_data = 'P1'
        
        self.assertEqual(data, reference_data)

    def test_pdb_in_results_1(self):
        """Test case for LatticeSearch.pdb_in_results"""

        from simbad.lattice.latticescore import LatticeSearchResult
        results = []
        score = LatticeSearchResult("4Z88", "", "", "62.420,62.492,163.64,90.234,90.000,119.962", 4.164, 3.892, 0.272, 23300, 0.1800)
        results.append(score)

        self.assertTrue(self.LS.pdb_in_results("4Z88", results)) 
        
if __name__ == "__main__":
    unittest.main()

