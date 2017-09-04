"""Test functions for util.mtz_util"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.command_line import ccp4_root
from simbad.util import mtz_util


class Test(unittest.TestCase):
    """Unit test"""
    
    def test_crystal_data_1(self):
        """Test case for mtz_util.crystal_data"""
        
        input_mtz = os.path.join(ccp4_root(), "examples", "toxd", "toxd.mtz")
        data = mtz_util.crystal_data(input_mtz)
        
        reference_data = ('P212121', '2.30', '73.58 38.73 23.19 90.00 90.00 90.00')
        
        self.assertEqual(data, reference_data)
        
    def test_crystal_data_2(self):
        """Test case for mtz_util.crystal_data"""
        
        input_mtz = os.path.join(ccp4_root(), "examples", "rnase", "rnase25.mtz")
        data = mtz_util.crystal_data(input_mtz)
        
        reference_data = ('P212121', '2.50', '64.90 78.32 38.79 90.00 90.00 90.00')
        
        self.assertEqual(data, reference_data)
        
    def test_get_labels_1(self):
        """Test case for mtz_util.get_labels"""
        
        input_mtz = os.path.join(ccp4_root(), "examples", "toxd", "toxd.mtz")
        data = mtz_util.get_labels(input_mtz)
        
        reference_data = ('FTOXD3', 'SIGFTOXD3', 'ANAU20', 'SIGANAU20', 'FreeR_flag')
        
        self.assertEqual(data, reference_data)
        
    def test_get_labels_2(self):
        """Test case for mtz_util.get_labels"""
        
        input_mtz = os.path.join(ccp4_root(), "examples", "rnase", "rnase25.mtz")
        data = mtz_util.get_labels(input_mtz)
        
        reference_data = ('FNAT', 'SIGFNAT', None, None, 'FreeR_flag')
        
        self.assertEqual(data, reference_data)
        
if __name__ == "__main__":
    unittest.main()