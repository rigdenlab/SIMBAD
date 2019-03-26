"""Test functions for util.mtz_util"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.command_line import CCP4RootDirectory 
from simbad.util import mtz_util

CCP4ROOT = str(CCP4RootDirectory())


class Test(unittest.TestCase):
    """Unit test"""
    
    def test_crystal_data_1(self):
        """Test case for mtz_util.crystal_data"""
        
        input_mtz = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.mtz")
        data = mtz_util.crystal_data(input_mtz)
        
        reference_data = ('P212121', 
                2.300205240684743, 
                (73.58200073242188, 38.733001708984375, 23.18899917602539, 90.0, 90.0, 90.0)
        )
        
        self.assertEqual(data, reference_data)
        
    def test_crystal_data_2(self):
        """Test case for mtz_util.crystal_data"""
        
        input_mtz = os.path.join(CCP4ROOT, "examples", "rnase", "rnase25.mtz")
        data = mtz_util.crystal_data(input_mtz)
        
        reference_data = ('P212121', 
                2.4999665357495098,
                (64.89700317382812, 78.322998046875, 38.79199981689453, 90.0, 90.0, 90.0)
        )
        
        self.assertEqual(data, reference_data)
        
    def test_get_labels_1(self):
        """Test case for mtz_util.get_labels"""
        
        input_mtz = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        temp_log = os.path.join(os.getcwd(), "input.log")
        mtz_util.ctruncate(input_mtz, temp_mtz)
        mtz_labels = mtz_util.GetLabels(temp_mtz)
        os.remove(temp_mtz)
        os.remove(temp_log)

        data = (mtz_labels.f, mtz_labels.sigf, mtz_labels.free)

        reference_data = ('FTOXD3', 'SIGFTOXD3', 'FreeR_flag')
        
        self.assertEqual(data, reference_data)
        
    def test_get_labels_2(self):
        """Test case for mtz_util.get_labels"""
        
        input_mtz = os.path.join(CCP4ROOT, "examples", "rnase", "rnase25F+F-.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        temp_log = os.path.join(os.getcwd(), "input.log")
        mtz_util.ctruncate(input_mtz, temp_mtz)
        mtz_labels = mtz_util.GetLabels(temp_mtz)
        os.remove(temp_mtz)
        os.remove(temp_log)

        data = (mtz_labels.f, mtz_labels.sigf,
                mtz_labels.fplus, mtz_labels.sigfplus,
                mtz_labels.fminus, mtz_labels.sigfminus,
                mtz_labels.free)
        reference_data = ('FNAT', 'SIGFNAT',
                          'FPTNCD25(+)', 'SIGFPTNCD25(+)',
                          'FPTNCD25(-)', 'SIGFPTNCD25(-)',
                          'FreeR_flag')
        
        self.assertEqual(data, reference_data)

    def test_change_space_group_1(self):
        """Test case for mtz_util.ExperimentalData.change_space_group"""
        input_mtz = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        mtz_util.reindex(input_mtz, temp_mtz, '18')

        data = mtz_util.crystal_data(temp_mtz)
        reference_data = ('P21212',
                          2.300205240684743,
                          (73.58200073242188, 38.733001708984375, 23.18899917602539, 90.0, 90.0, 90.0)
                          )
        self.assertEqual(data, reference_data)
        
if __name__ == "__main__":
    unittest.main()
