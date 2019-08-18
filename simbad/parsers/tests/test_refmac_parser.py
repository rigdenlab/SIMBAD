"""Test functions for simbad.parsers.refmac_parser"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import os
import tempfile
import unittest

from simbad.parsers import refmac_parser


class Test(unittest.TestCase):
    def test_refmac_parser_1(self):
        content = """
 $TEXT:Result: $$ Final results $$
                      Initial    Final
           R factor    0.2666   0.2666
             R free    0.2659   0.2659
     Rms BondLength    0.0290   0.0290
      Rms BondAngle    3.3998   3.3998
     Rms ChirVolume    0.2495   0.2495
        """
        refmac_log = tempfile.NamedTemporaryFile("w", delete=False)
        refmac_log.write(content)
        refmac_log.close()

        rp = refmac_parser.RefmacParser(refmac_log.name)
        self.assertEqual(rp.init_r_free, 0.2659)
        self.assertEqual(rp.init_r_fact, 0.2666)
        self.assertEqual(rp.final_r_free, 0.2659)
        self.assertEqual(rp.final_r_fact, 0.2666)


if __name__ == "__main__":
    unittest.main()
