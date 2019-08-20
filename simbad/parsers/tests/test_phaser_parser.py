"""Test functions for simbad.parsers.phaser_parser"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import os
import tempfile
import unittest

from simbad.parsers import phaser_parser


class Test(unittest.TestCase):
    def test_phaser_parser_1(self):
        content = """
   SOLU SET  RF*0
   SOLU SET
   SOLU SET
   SOLU SET  RF*0 TF*0 LLG=    1419 TFZ==28.9 PAK=0 LLG=1419 TFZ==28.9
        """

        phaser_log = tempfile.NamedTemporaryFile("w", delete=False)
        phaser_log.write(content)
        phaser_log.close()

        pp = phaser_parser.PhaserParser(phaser_log.name)

        self.assertEqual(pp.llg, 1419)
        self.assertEqual(pp.tfz, 28.9)
        self.assertEqual(pp.rfz, None)


if __name__ == "__main__":
    unittest.main()
