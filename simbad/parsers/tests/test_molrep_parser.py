import tempfile
import unittest

from simbad.parsers.molrep_parser import MolrepParser


class Test(unittest.TestCase):
    def test_molrep_parser_1(self):
        content = """
 corrF =   0.7025
 TF/sig       =     9.34
 Final CC     =   0.7025
 Packing_Coef =   1.0000
 Contrast     =    12.39

 After stick correction:
  Move closer to origin
  I_sym_operator    :    4
  new position(frac):  -0.111   0.168   0.032


  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
    1   1   1   89.00   90.61  179.24 -0.111  0.168  0.032   9.34  0.417  0.703
        """
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, 0.7025)
        self.assertEqual(mp.tfscore, 9.34)
        self.assertEqual(mp.contrast, 12.39)
        self.assertEqual(mp.wrfac, 0.417)

    def test_molrep_parser_2(self):
        content = ""
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, None)
        self.assertEqual(mp.tfscore, None)
        self.assertEqual(mp.contrast, None)
        self.assertEqual(mp.wrfac, None)

    def test_molrep_parser_3(self):
        content = """
 corrF =   0.8164
 TF/sig       =    23.04
 Final CC     =   0.8164
 Packing_Coef =   1.0000
 Contrast     =    26.24

 After stick correction:
  Move closer to origin
  I_sym_operator    :    6
  new position(frac):   0.298   0.331  -0.046


  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
    1   1   1    0.00    0.00   59.84  0.298  0.331 -0.046  23.04  0.374  0.816
        """
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, 0.8164)
        self.assertEqual(mp.tfscore, 23.04)
        self.assertEqual(mp.contrast, 26.24)
        self.assertEqual(mp.wrfac, 0.374)


if __name__ == "__main__":
    unittest.main()
