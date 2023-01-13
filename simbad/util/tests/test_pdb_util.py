"""Test functions for simbad.util.pdb_util"""

__author__ = "Adam Simpkin"
__date__ = "19 Jan 2018"

import os
import numpy as np
import sys
import unittest

from simbad.util.pdb_util import PdbStructure

if sys.version_info.major < 3:
    from urllib2 import urlopen
    from urllib2 import URLError
else:
    from urllib.request import urlopen
    from urllib.error import URLError

try:
    ROOT_DIR = os.environ['SIMBAD_ROOT']
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "test_data")
except KeyError:
    from simbad.command_line import CCP4RootDirectory
    ROOT_DIR = str(CCP4RootDirectory())
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "examples")


def internet_on():
    try:
        urlopen("http://www.google.com", timeout=1)
    except URLError:
        return False
    else:
        return True


class Test(unittest.TestCase):
    """Unit test"""

    def test_molecular_weight_1(self):
        """Test case for PdbStructure.molecular_weight"""

        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 6855.978639999951

        self.assertAlmostEqual(np.round(data, 0), np.round(reference_data, 0))

    def test_molecular_weight_2(self):
        """Test case for PdbStructure.molecular_weight"""

        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 21163.001080955823

        self.assertAlmostEqual(np.round(data, 0), np.round(reference_data, 0))

    def test_molecular_weight_3(self):
        """Test case for PdbStructure.molecular_weight"""

        input_model = os.path.join(EXAMPLE_DIR, "data", "3a22.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.molecular_weight
        reference_data = 128535.91044924183

        self.assertAlmostEqual(np.round(data, 0), np.round(reference_data, 0))

    def test_calculate_integration_box(self):
        """Test case for PdbStructure.integration_box"""

        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.integration_box
        reference_data = (48.306749999999994, 56.73474999999999, 48.589749999999995, 19.84575)

        self.assertEqual(data, reference_data)

    def test_nchains_1(self):
        """Test case for PdbStructure.nchains"""
        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.nchains
        reference_data = 1

        self.assertEqual(data, reference_data)

    def test_nchains_2(self):
        """Test case for PdbStructure.nchains"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.nchains
        reference_data = 2

        self.assertEqual(data, reference_data)

    def test_nres_1(self):
        """Test case for PdbStructure.nres"""
        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.nres
        reference_data = 59

        self.assertEqual(data, reference_data)

    def test_nres_2(self):
        """Test case for PdbStructure.nres"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        data = pdb_struct.nres
        reference_data = 192

        self.assertEqual(data, reference_data)

    def test_keep_first_chain(self):
        """Test case for PdbStructure.keep_first_chain_only"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        pdb_struct.keep_first_chain_only()
        data = pdb_struct.nchains
        reference_data = 1
        self.assertEqual(data, reference_data)

        data = pdb_struct.structure[0][0].name
        reference_data = 'A'
        self.assertEqual(data, reference_data)

    def test_select_chain_by_idx_1(self):
        """Test case for PdbStructure.select_chain_by_idx"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        pdb_struct.select_chain_by_idx(0)
        data = pdb_struct.structure[0][0].name
        reference_data = 'A'
        self.assertEqual(data, reference_data)

        data = pdb_struct.nchains
        reference_data = 1
        self.assertEqual(data, reference_data)

    def test_select_chain_by_idx_2(self):
        """Test case for PdbStructure.select_chain_by_idx"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        pdb_struct.select_chain_by_idx(1)
        data = pdb_struct.structure[0][0].name
        reference_data = 'B'
        self.assertEqual(data, reference_data)

        data = pdb_struct.nchains
        reference_data = 1
        self.assertEqual(data, reference_data)

    def test_select_chain_by_id_1(self):
        """Test case for PdbStructure.select_chain_by_id"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        pdb_struct.select_chain_by_id('A')
        data = pdb_struct.structure[0][0].name
        reference_data = 'A'
        self.assertEqual(data, reference_data)

        data = pdb_struct.nchains
        reference_data = 1
        self.assertEqual(data, reference_data)

    def test_select_chain_by_id_2(self):
        """Test case for PdbStructure.select_chain_by_id"""
        input_model = os.path.join(EXAMPLE_DIR, "rnase", "rnase.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        pdb_struct.select_chain_by_id('B')
        data = pdb_struct.structure[0][0].name
        reference_data = 'B'
        self.assertEqual(data, reference_data)

        data = pdb_struct.nchains
        reference_data = 1
        self.assertEqual(data, reference_data)

    def test_select_residues(self):
        """Test case for PdbStructure.select_residues"""
        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        pdb_struct = PdbStructure.from_file(input_model)
        seqid_range = range(0, 5)
        pdb_struct.select_residues(to_keep_idx=seqid_range)
        data = pdb_struct.nres
        reference_data = 5
        self.assertEqual(data, reference_data)

        data = [res.seqid.num for res in pdb_struct.structure[0][0]]
        reference_data = [1, 2, 3, 4, 5]
        self.assertListEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_standardise_1(self):
        """Test case for PdbStructure.standardise"""

        if internet_on():
            pdb_code = "1RGH"
            pdb_struct = PdbStructure.from_pdb_code(pdb_code)
            pdb_struct.standardize()
            pdb_struct.save("tmp.pdb")

            data = []
            count = 0
            with open("tmp.pdb", "r") as f:
                for line in f:
                    if count < 10:
                        data.append(line.replace('\r\n', '\n'))
                        count += 1
            os.unlink("tmp.pdb")
            
            reference_data = [
                    'CRYST1   64.820   78.560   39.050  90.00  90.00  90.00 P 21 21 21    0          \n', 
                    'ATOM      1  N   ASP A   1      45.025  12.898   8.789  1.00 20.71           N  \n',
                    'ATOM      2  CA  ASP A   1      45.101  12.554   7.385  1.00 21.13           C  \n',
                    'ATOM      3  C   ASP A   1      46.569  12.441   6.994  1.00 18.35           C  \n',
                    'ATOM      4  O   ASP A   1      47.443  13.016   7.624  1.00 23.34           O  \n',
                    'ATOM      5  CB  ASP A   1      44.385  13.622   6.515  1.00 27.06           C  \n',
                    'ATOM      6  CG  ASP A   1      42.911  13.800   6.841  1.00 40.10           C  \n',
                    'ATOM      7  OD1 ASP A   1      42.374  12.978   7.610  1.00 48.03           O  \n',
                    'ATOM      8  OD2 ASP A   1      42.292  14.768   6.325  1.00 52.96           O  \n',
                    'ATOM      9  N   VAL A   2      46.810  11.649   5.958  1.00 20.56           N  \n'
            ]
            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_1")

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_standardise_2(self):
        """Test case for PdbStructure.standardise"""

        if internet_on():
            pdb_code = "1LNI"
            pdb_struct = PdbStructure.from_pdb_code(pdb_code)
            pdb_struct.standardize()
            pdb_struct.save("tmp.pdb")

            data = []
            count = 0
            with open("tmp.pdb", "r") as f:
                for line in f:
                    if count < 10:
                        data.append(line.replace('\r\n', '\n'))
                        count += 1
            os.unlink("tmp.pdb")

            reference_data = [
                    'CRYST1   64.200   77.800   38.280  90.00  90.00  90.00 P 21 21 21    0          \n',
                    'ATOM      1  N   ASP A   1      44.433  13.208   8.560  1.00 13.64           N  \n',
                    'ATOM      2  CA  ASP A   1      44.447  13.009   7.129  1.00 14.33           C  \n',
                    'ATOM      3  C   ASP A   1      45.887  12.787   6.694  1.00 13.28           C  \n',
                    'ATOM      4  O   ASP A   1      46.792  13.384   7.254  1.00 15.83           O  \n',
                    'ATOM      5  CB  ASP A   1      43.901  14.240   6.407  1.00 18.69           C  \n',
                    'ATOM      6  CG  ASP A   1      42.483  14.625   6.750  1.00 22.80           C  \n',
                    'ATOM      7  OD1 ASP A   1      41.912  14.023   7.666  1.00 25.27           O  \n',
                    'ATOM      8  OD2 ASP A   1      41.945  15.532   6.100  1.00 31.97           O  \n',
                    'ATOM      9  N   VAL A   2      46.061  11.960   5.672  1.00 14.23           N  \n'
            ]
            

            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_2")


if __name__ == "__main__":
    unittest.main()
