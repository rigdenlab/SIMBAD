"""Test functions for simbad.util.pdb_util"""

__author__ = "Adam Simpkin"
__date__ = "19 Jan 2018"

import os
import unittest
import simbad.util.pdb_util
import urllib2

from simbad.command_line import CCP4RootDirectory

CCP4ROOT = str(CCP4RootDirectory())


def internet_on():
    try:
        urllib2.urlopen('http://www.google.com', timeout=1)
        return True
    except urllib2.URLError:
        return False


class Test(unittest.TestCase):
    """Unit test"""

    def test_calculate_integration_box(self):
        """Test case for PdbStructure.integration_box"""
 
        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure()
        pdb_struct.from_file(input_model)
        data = pdb_struct.integration_box
        reference_data = (48.306749999999994, 56.73474999999999, 48.589749999999995, 19.84575)

        self.assertEqual(data, reference_data)
    
    def test_molecular_weight_1(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure()
        pdb_struct.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 7147.307000000012
        
        self.assertAlmostEqual(data, reference_data)
        
    def test_molecular_weight_2(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "rnase", "rnase.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure()
        pdb_struct.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 21333.267999999967
        
        self.assertAlmostEqual(data, reference_data)
        
    def test_molecular_weight_3(self):
        """Test case for PdbStructure.molecular_weight"""
        
        input_model = os.path.join(CCP4ROOT, "examples", "data", "3a22.pdb")
        pdb_struct = simbad.util.pdb_util.PdbStructure()
        pdb_struct.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 132496.16599998937
        
        self.assertAlmostEqual(data, reference_data)

    def test_standardise_1(self):
        """Test case for PdbStructure.standardise"""

        if internet_on():
            pdb_code = "1RGH"
            pdb_struct = simbad.util.pdb_util.PdbStructure()
            pdb_struct.from_pdb_code(pdb_code)
            pdb_struct.standardize()
            pdb_struct.save("tmp.pdb")

            data = []
            count = 0
            with open("tmp.pdb", 'r') as f:
                for line in f:
                    if count < 10:
                        data.append(line)
                        count += 1
            os.unlink("tmp.pdb")

            reference_data = ['CRYST1   64.820   78.560   39.050  90.00  90.00  90.00 P 21 21 21\n',
                              'SCALE1      0.015427  0.000000  0.000000        0.00000\n',
                              'SCALE2      0.000000  0.012729  0.000000        0.00000\n',
                              'SCALE3      0.000000  0.000000  0.025608        0.00000\n',
                              'ATOM      1  N   ASP A   1      44.945  12.883   8.824  1.00 19.69           N\n',
                              'ATOM      2  CA  ASP A   1      45.082  12.376   7.484  1.00 22.49           C\n',
                              'ATOM      3  C   ASP A   1      46.551  12.322   7.040  1.00 20.41           C\n',
                              'ATOM      4  O   ASP A   1      47.413  12.994   7.593  1.00 22.53           O\n',
                              'ATOM      5  CB  ASP A   1      44.383  13.428   6.514  1.00 34.50           C\n',
                              'ATOM      6  CG  ASP A   1      42.999  13.726   7.124  1.00 45.87           C\n']

            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_1")

    def test_standardise_2(self):
        """Test case for PdbStructure.standardise"""

        if internet_on():
            pdb_code = "1LNI"
            pdb_struct = simbad.util.pdb_util.PdbStructure()
            pdb_struct.from_pdb_code(pdb_code)
            pdb_struct.standardize()
            pdb_struct.save("tmp.pdb")

            data = []
            count = 0
            with open("tmp.pdb", 'r') as f:
                for line in f:
                    if count < 10:
                        data.append(line)
                        count += 1
            os.unlink("tmp.pdb")

            reference_data = ['CRYST1   64.200   77.800   38.280  90.00  90.00  90.00 P 21 21 21\n',
                              'SCALE1      0.015576  0.000000  0.000000        0.00000\n',
                              'SCALE2      0.000000  0.012853  0.000000        0.00000\n',
                              'SCALE3      0.000000  0.000000  0.026123        0.00000\n',
                              'ATOM      1  N   ASP A   1      44.437  13.199   8.551  1.00 14.05           N\n',
                              'ATOM      2  CA  ASP A   1      44.452  12.952   7.124  1.00 14.52           C\n',
                              'ATOM      3  C   ASP A   1      45.912  12.766   6.723  1.00 12.79           C\n',
                              'ATOM      4  O   ASP A   1      46.786  13.402   7.273  1.00 17.99           O\n',
                              'ATOM      5  CB  ASP A   1      43.905  14.169   6.325  1.00 22.32           C\n',
                              'ATOM      6  CG  ASP A   1      42.555  14.639   6.688  1.00 28.30           C\n']

            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_2")

if __name__ == "__main__":
    unittest.main()
