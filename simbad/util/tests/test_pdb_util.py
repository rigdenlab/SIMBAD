"""Test functions for simbad.util.pdb_util"""

__author__ = "Adam Simpkin"
__date__ = "19 Jan 2018"

import os
import unittest
import simbad.util.pdb_util

from simbad.command_line import CCP4RootDirectory

CCP4ROOT = str(CCP4RootDirectory())

class Test(unittest.TestCase):
    """Unit test"""

    def test_calculate_integration_box(self):
        """Test case for PdbStructure.integration_box"""
 
        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure(input_model)
        data = pdb_struct.integration_box
        reference_data = (48.306749999999994, 56.73474999999999, 48.589749999999995, 19.84575)

        self.assertEqual(data, reference_data)
    
    def test_molecular_weight_1(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 7147.307000000012
        
        self.assertEqual(data, reference_data)
        
    def test_molecular_weight_2(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "rnase", "rnase.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 21333.267999999967
        
        self.assertEqual(data, reference_data)
        
    def test_molecular_weight_3(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "data", "3a22.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 132496.16599998937
        
        self.assertEqual(data, reference_data)

if __name__ == "__main__":
    unittest.main()
