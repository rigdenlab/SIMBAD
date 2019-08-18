"""Test functions for simbad.db"""

__author__ = "Adam Simpkin"
__date__ = "19 Jan 2018"

import os
import unittest
import simbad
import simbad.db

from simbad.command_line import CCP4RootDirectory

CCP4ROOT = str(CCP4RootDirectory())


class Test(unittest.TestCase):
    """Unit test"""

    def test_from_dat(self):
        """Test case for simbad.db._from_dat"""
        simbad_db = simbad.CONTAMINANT_MODELS
        input_dat = os.path.join(simbad_db, "CHICK", "LYSC_CHICK", "P6122", "2fbb.dat")
        with open(input_dat, "r") as f_in:
            output_str = simbad.db._from_dat(f_in)

        data = output_str.split("\n")[0]
        reference_data = "HEADER    HYDROLASE                               09-DEC-05   2FBB              "

        self.assertEqual(data, reference_data)

    def test_to_dat(self):
        """Test case for simbad.db._to_dat"""
        input_pdb = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        output_dat = os.path.join(os.getcwd(), "test.dat")
        with open(input_pdb, "r") as f_in, open(output_dat, "wb") as f_out:
            f_out.write(simbad.db._to_dat(f_in))

        self.assertTrue(simbad.db.is_valid_dat(output_dat))

        with open(input_pdb, "r") as f_in:
            input_str = f_in.read().split("\n")[0]

        output_str = simbad.db.read_dat(output_dat).split("\n")[0]

        self.assertEqual(input_str, output_str)

        os.unlink(output_dat)

    def test_str_to_dat(self):
        """Test case for simbad.db._str_to_dat"""
        input_str = "TEST"
        output_dat = os.path.join(os.getcwd(), "test.dat")
        with open(output_dat, "wb") as f_out:
            f_out.write(simbad.db._str_to_dat(input_str))

        self.assertTrue(simbad.db.is_valid_dat(output_dat))

        output_str = simbad.db.read_dat(output_dat)

        self.assertEqual(input_str, output_str)

        os.unlink(output_dat)

    def test_find_simbad_dat_files(self):
        """Test case for simbad.db.find_simbad_dat_files"""
        simbad_db = simbad.CONTAMINANT_MODELS
        test_dat_db = os.path.join(simbad_db, "CHICK", "LYSC_CHICK", "P6122")
        data = os.path.basename(simbad.db.find_simbad_dat_files(test_dat_db)[0])
        reference_data = "2fbb.dat"

        self.assertEqual(data, reference_data)

    def test_convert_pdb_to_dat(self):
        """Test case for simbad.db.convert_pdb_to_dat"""
        input_pdb = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        output_dat = os.path.join(os.getcwd(), "test.dat")
        simbad.db.convert_pdb_to_dat(input_pdb, output_dat)

        self.assertTrue(os.path.isfile(output_dat))
        self.assertTrue(simbad.db.is_valid_dat(output_dat))

        os.unlink(output_dat)

    def test_convert_dat_to_pdb(self):
        """Test case for simbad.db.convert_dat_to_pdb"""
        simbad_db = simbad.CONTAMINANT_MODELS
        input_dat = os.path.join(simbad_db, "CHICK", "LYSC_CHICK", "P6122", "2fbb.dat")
        output_pdb = os.path.join(os.getcwd(), "test.pdb")
        simbad.db.convert_dat_to_pdb(input_dat, output_pdb)

        self.assertTrue(os.path.isfile(output_pdb))
        os.unlink(output_pdb)

    def test_is_valid_dat(self):
        """Test case for simbad.db.is_valid_dat"""
        simbad_db = simbad.CONTAMINANT_MODELS
        input_dat = os.path.join(simbad_db, "CHICK", "LYSC_CHICK", "P6122", "2fbb.dat")
        data = simbad.db.is_valid_dat(input_dat)

        self.assertTrue(data)

    def test_read_dat(self):
        """Test case for simbad.db.read_dat"""
        simbad_db = simbad.CONTAMINANT_MODELS
        input_dat = os.path.join(simbad_db, "CHICK", "LYSC_CHICK", "P6122", "2fbb.dat")
        output_str = simbad.db.read_dat(input_dat)

        data = output_str.split("\n")[0]
        reference_data = "HEADER    HYDROLASE                               09-DEC-05   2FBB              "

        self.assertEqual(data, reference_data)
