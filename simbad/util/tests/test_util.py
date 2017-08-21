"""Test functions for simbad.util"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import os
import unittest
import simbad.util

from simbad.command_line import ccp4_root


class Test(unittest.TestCase):
    """Unit test"""
    
    def test_molecular_weight_1(self):
        """Test case for simbad.util.molecular_weight"""
        
        input_model = os.path.join(ccp4_root(), "examples", "toxd", "toxd.pdb")
        data = simbad.util.molecular_weight(input_model)
        reference_data = 7147.307000000012
        
        self.assertEqual(data, reference_data)
        
    def test_molecular_weight_2(self):
        """Test case for simbad.util.molecular_weight"""
        
        input_model = os.path.join(ccp4_root(), "examples", "rnase", "rnase.pdb")
        data = simbad.util.molecular_weight(input_model)
        reference_data = 21333.267999999967
        
        self.assertEqual(data, reference_data)
        
    def test_molecular_weight_3(self):
        """Test case for simbad.util.molecular_weight"""
        
        input_model = os.path.join(ccp4_root(), "examples", "data", "3a22.pdb")
        data = simbad.util.molecular_weight(input_model)
        reference_data = 132496.16599998937
        
        self.assertEqual(data, reference_data)

if __name__ == "__main__":
    unittest.main()