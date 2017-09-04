"""Test functions for simbad.parsers.molrep_parser"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import os
import tempfile
import unittest

from simbad.parsers import molrep_parser

class Test(unittest.TestCase):
    """Unit test"""
    
    def test_molrep_parser(self):
        """Test case for simbad.parsers.molrep_parser.MolrepParser.parse"""
        
        # Make tmp file containing key lines from MOLREP log
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write("""
  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
    1  11   1   60.73 -114.70  137.64  0.202  0.150 -0.191   3.77  0.613  0.187""")
        molrep_log.close()
        
        MP = molrep_parser.MolrepParser(molrep_log.name)
        MP.parse(molrep_log.name)
        
        self.assertEqual(MP.score, 0.187)
        self.assertEqual(MP.tfscore, 3.77)
        self.assertEqual(MP.wrfac, 0.613)

if __name__ == "__main__":
    unittest.main()