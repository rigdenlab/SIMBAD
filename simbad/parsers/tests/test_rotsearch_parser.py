"""Test functions for simbad.parsers.rotsearch_parser"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import tempfile
import unittest

from simbad.parsers import rotsearch_parser


class Test(unittest.TestCase):
    def test_rotsearch_parser_1(self):
        content = """
          CROSS:Generate HKLPCK1 from MODEL FRAGMENT 1
            ITAB  ALPHA    BETA   GAMMA    TX      TY      TZ    CC_F RF_F CC_I CC_P Icp  Z_F  Z_M No
 SOLUTIONRCD   1    0.46   88.40  116.18  0.0000  0.0000  0.0000  7.9 60.7 10.3  3.8   1  3.9  2.2  1
 SOLUTIONRCD   1   58.40   90.00  296.07  0.0000  0.0000  0.0000  7.9 60.7 10.4  5.3   2  3.9  3.1  3
 SOLUTIONRCD   1   23.21   90.00   36.25  0.0000  0.0000  0.0000  7.8 60.9 10.5  4.8   5  3.8  2.8  1
        """
        rotsearch_log = tempfile.NamedTemporaryFile("w", delete=False)
        rotsearch_log.write(content)
        rotsearch_log.close()

        rp = rotsearch_parser.AmoreRotsearchParser(rotsearch_log.name)
        self.assertEqual(rp.alpha, 0.46)
        self.assertEqual(rp.beta, 88.40)
        self.assertEqual(rp.gamma, 116.18)
        self.assertEqual(rp.cc_f, 7.9)
        self.assertEqual(rp.rf_f, 60.7)
        self.assertEqual(rp.cc_i, 10.3)
        self.assertEqual(rp.cc_p, 3.8)
        self.assertEqual(rp.icp, 1)
        self.assertEqual(rp.cc_f_z_score, 3.9)
        self.assertEqual(rp.cc_p_z_score, 2.2)
        self.assertEqual(rp.num_of_rot, 1)

    def test_rotsearch_parser_2(self):
        content = """   Rotation Function Results
   Top1: ENSEMBLE PDB EULER 54.157 37.995 130.126 RF=4.8 RFZ=3.27

   Rotation Function Table: PDB
   ----------------------------
   (Z-scores from Fast Rotation Function)
   #SET        Top    (Z)      Second    (Z)       Third    (Z)
   1          4.82   3.27        3.94   3.17        2.76   3.04
        """

        rotsearch_log = tempfile.NamedTemporaryFile("w", delete=False)
        rotsearch_log.write(content)
        rotsearch_log.close()

        rp = rotsearch_parser.PhaserRotsearchParser(rotsearch_log.name)
        self.assertEqual(rp.llg, 4.82)
        self.assertEqual(rp.rfz, 3.27)


if __name__ == "__main__":
    unittest.main()
