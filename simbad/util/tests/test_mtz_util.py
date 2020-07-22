"""Test functions for util.mtz_util"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.parsers import mtz_parser
from simbad.util import mtz_util

SIMBAD_ROOT = os.environ['SIMBAD_ROOT']


class Test(unittest.TestCase):
    """Unit test"""

    def test_crystal_data_1(self):
        """Test case for mtz_util.crystal_data"""

        input_mtz = os.path.join(SIMBAD_ROOT, "test_data", "toxd.mtz")
        mtz_obj = mtz_parser.MtzParser(input_mtz)
        data = (
            "".join(mtz_obj.spacegroup_symbol.encode("ascii").split()),
            round(mtz_obj.resolution, 2),
            round(mtz_obj.cell.a, 2),
            round(mtz_obj.cell.b, 2),
            round(mtz_obj.cell.c, 2),
            round(mtz_obj.cell.alpha, 2),
            round(mtz_obj.cell.beta, 2),
            round(mtz_obj.cell.gamma, 2)
        )

        reference_data = (
            "P212121",
            2.30,
            73.58,
            38.73,
            23.19,
            90.00,
            90.00,
            90.00
        )

        self.assertEqual(data[0], reference_data[0])
        self.assertAlmostEqual(data[1], reference_data[1])
        self.assertAlmostEqual(data[2], reference_data[2])
        self.assertAlmostEqual(data[3], reference_data[3])
        self.assertAlmostEqual(data[4], reference_data[4])
        self.assertAlmostEqual(data[5], reference_data[5])
        self.assertAlmostEqual(data[6], reference_data[6])
        self.assertAlmostEqual(data[7], reference_data[7])

    def test_crystal_data_2(self):
        """Test case for mtz_util.crystal_data"""

        input_mtz = os.path.join(SIMBAD_ROOT, "test_data", "rnase25.mtz")
        mtz_obj = mtz_parser.MtzParser(input_mtz)
        data = (
            "".join(mtz_obj.spacegroup_symbol.encode("ascii").split()),
            round(mtz_obj.resolution, 2),
            round(mtz_obj.cell.a, 2),
            round(mtz_obj.cell.b, 2),
            round(mtz_obj.cell.c, 2),
            round(mtz_obj.cell.alpha, 2),
            round(mtz_obj.cell.beta, 2),
            round(mtz_obj.cell.gamma, 2)
        )

        reference_data = (
            "P212121",
            2.50,
            64.90,
            78.32,
            38.79,
            90.00,
            90.00,
            90.00
        )

        self.assertEqual(data[0], reference_data[0])
        self.assertAlmostEqual(data[1], reference_data[1])
        self.assertAlmostEqual(data[2], reference_data[2])
        self.assertAlmostEqual(data[3], reference_data[3])
        self.assertAlmostEqual(data[4], reference_data[4])
        self.assertAlmostEqual(data[5], reference_data[5])
        self.assertAlmostEqual(data[6], reference_data[6])
        self.assertAlmostEqual(data[7], reference_data[7])

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_get_labels_1(self):
        """Test case for mtz_util.get_labels"""

        input_mtz = os.path.join(SIMBAD_ROOT, "test_data", "toxd.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        temp_log = os.path.join(os.getcwd(), "input.log")
        mtz_util.ctruncate(input_mtz, temp_mtz)
        mtz_obj = mtz_parser.MtzParser(temp_mtz)
        mtz_obj.parse()
        os.remove(temp_mtz)
        os.remove(temp_log)

        data = (mtz_obj.f, mtz_obj.sigf, mtz_obj.free)

        reference_data = ("FTOXD3", "SIGFTOXD3", "FreeR_flag")

        self.assertEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_get_labels_2(self):
        """Test case for mtz_util.get_labels"""

        input_mtz = os.path.join(SIMBAD_ROOT, "test_data", "rnase25F+F-.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        temp_log = os.path.join(os.getcwd(), "input.log")
        mtz_util.ctruncate(input_mtz, temp_mtz)
        mtz_obj = mtz_parser.MtzParser(temp_mtz)
        mtz_obj.parse()
        os.remove(temp_mtz)
        os.remove(temp_log)

        data = (mtz_obj.f, mtz_obj.sigf, mtz_obj.f_plus, mtz_obj.sigf_plus, mtz_obj.f_minus,
                mtz_obj.sigf_minus, mtz_obj.free)
        reference_data = ("FNAT", "SIGFNAT", "FPTNCD25(+)", "SIGFPTNCD25(+)",
                          "FPTNCD25(-)", "SIGFPTNCD25(-)", "FreeR_flag")

        self.assertEqual(data, reference_data)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_change_space_group_1(self):
        """Test case for mtz_util.ExperimentalData.change_space_group"""
        input_mtz = os.path.join(SIMBAD_ROOT, "test_data", "toxd.mtz")
        temp_mtz = os.path.join(os.getcwd(), "input.mtz")
        mtz_util.reindex(input_mtz, temp_mtz, "18")
        mtz_obj = mtz_parser.MtzParser(temp_mtz)
        data = (
            "".join(mtz_obj.spacegroup_symbol.encode("ascii").split()),
            round(mtz_obj.resolution, 2),
            round(mtz_obj.cell.a, 2),
            round(mtz_obj.cell.b, 2),
            round(mtz_obj.cell.c, 2),
            round(mtz_obj.cell.alpha, 2),
            round(mtz_obj.cell.beta, 2),
            round(mtz_obj.cell.gamma, 2)
        )
        reference_data = (
            "P21212",
            2.30,
            73.58,
            38.73,
            23.19,
            90.00,
            90.00,
            90.00
        )

        self.assertEqual(data[0], reference_data[0])
        self.assertAlmostEqual(data[1], reference_data[1])
        self.assertAlmostEqual(data[2], reference_data[2])
        self.assertAlmostEqual(data[3], reference_data[3])
        self.assertAlmostEqual(data[4], reference_data[4])
        self.assertAlmostEqual(data[5], reference_data[5])
        self.assertAlmostEqual(data[6], reference_data[6])
        self.assertAlmostEqual(data[7], reference_data[7])


if __name__ == "__main__":
    unittest.main()
