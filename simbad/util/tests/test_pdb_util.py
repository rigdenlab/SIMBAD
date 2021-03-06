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
                "CRYST1   64.820   78.560   39.050  90.00  90.00  90.00 P 21 21 21               \n",
                "ATOM      1  N   ASP A   1      45.029  12.905   8.790  1.00 20.30           N  \n",
                "ATOM      2  CA  ASP A   1      45.094  12.556   7.371  1.00 20.54           C  \n",
                "ATOM      3  CB  ASP A   1      44.395  13.636   6.493  1.00 25.79           C  \n",
                "ATOM      4  CG  ASP A   1      43.026  14.022   7.007  1.00 37.35           C  \n",
                "ATOM      5  OD1 ASP A   1      42.364  13.133   7.618  1.00 46.11           O  \n",
                "ATOM      6  OD2 ASP A   1      42.616  15.195   6.812  1.00 52.03           O  \n",
                "ATOM      7  C   ASP A   1      46.579  12.457   6.999  1.00 18.38           C  \n",
                "ATOM      8  O   ASP A   1      47.443  13.026   7.624  1.00 23.12           O  \n",
                "ATOM      9  N   VAL A   2      46.822  11.665   5.980  1.00 20.00           N  \n",
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
                "CRYST1   64.200   77.800   38.280  90.00  90.00  90.00 P 21 21 21               \n",
                "ATOM      1  N   ASP A   1      44.446  13.220   8.556  1.00 13.50           N  \n",
                "ATOM      2  CA  ASP A   1      44.430  12.980   7.111  1.00 14.06           C  \n",
                "ATOM      3  CB  ASP A   1      43.907  14.215   6.393  1.00 19.18           C  \n",
                "ATOM      4  CG  ASP A   1      42.510  14.608   6.756  1.00 23.17           C  \n",
                "ATOM      5  OD1 ASP A   1      41.884  14.041   7.647  1.00 25.19           O  \n",
                "ATOM      6  OD2 ASP A   1      41.979  15.572   6.153  1.00 34.30           O  \n",
                "ATOM      7  C   ASP A   1      45.890  12.779   6.690  1.00 13.21           C  \n",
                "ATOM      8  O   ASP A   1      46.785  13.367   7.265  1.00 15.90           O  \n",
                "ATOM      9  N   VAL A   2      46.068  11.965   5.678  1.00 14.01           N  \n",
            ]

            self.assertEqual(data, reference_data)
        else:
            print("Skipping test_standardise_2")


if __name__ == "__main__":
    unittest.main()
