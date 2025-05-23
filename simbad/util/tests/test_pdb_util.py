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
                'CRYST1   64.820   78.560   39.050  90.00  90.00  90.00 P 21 21 21    8          \n', 
                'MTRIX1   1  0.972240  0.229670  0.044780      -33.28425    1                    \n', 
                'MTRIX2   1 -0.042920  0.363160 -0.930740       21.22728    1                    \n', 
                'MTRIX3   1 -0.230030  0.902980  0.362930       18.74696    1                    \n', 
                'ATOM      1  N   ASP A   1      44.945  12.883   8.824  1.00 19.69           N  \n', 
                'ATOM      2  CA  ASP A   1      45.082  12.376   7.484  1.00 22.49           C  \n', 
                'ATOM      3  C   ASP A   1      46.551  12.322   7.040  1.00 20.41           C  \n', 
                'ATOM      4  O   ASP A   1      47.413  12.994   7.593  1.00 22.53           O  \n', 
                'ATOM      5  CB  ASP A   1      44.383  13.428   6.514  1.00 34.50           C  \n', 
                'ATOM      6  CG  ASP A   1      42.999  13.726   7.124  1.00 45.87           C  \n'
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
                'CRYST1   64.200   77.800   38.280  90.00  90.00  90.00 P 21 21 21    8          \n', 
                'ATOM      1  N   ASP A   1      44.437  13.199   8.551  1.00 14.05           N  \n', 
                'ATOM      2  CA  ASP A   1      44.452  12.952   7.124  1.00 14.52           C  \n', 
                'ATOM      3  C   ASP A   1      45.912  12.766   6.723  1.00 12.79           C  \n', 
                'ATOM      4  O   ASP A   1      46.786  13.402   7.273  1.00 17.99           O  \n', 
                'ATOM      5  CB  ASP A   1      43.905  14.169   6.325  1.00 22.32           C  \n', 
                'ATOM      6  CG  ASP A   1      42.555  14.639   6.688  1.00 28.30           C  \n', 
                'ATOM      7  OD1 ASP A   1      42.038  15.680   6.186  1.00 42.93           O  \n', 
                'ATOM      8  OD2 ASP A   1      41.915  14.081   7.579  1.00 26.74           O  \n', 
                'ATOM      9  N   VAL A   2      46.081  11.961   5.697  1.00 14.42           N  \n'
            ]

            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_2")


if __name__ == "__main__":
    unittest.main()
