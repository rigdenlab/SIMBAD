import os
import tempfile
import unittest

from simbad.parsers.molrep_parser import MolrepParser


class Test(unittest.TestCase):
    def test_molrep_parser_1(self):
        content = """
  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
    1  11   1   60.73 -114.70  137.64  0.202  0.150 -0.191   3.77  0.613  0.187
        """
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, 0.187)
        self.assertEqual(mp.tfscore, 3.77)
        self.assertEqual(mp.wrfac, 0.613)

    def test_molrep_parser_2(self):
        content = ""
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, None)
        self.assertEqual(mp.tfscore, None)
        self.assertEqual(mp.wrfac, None)

    def test_molrep_parser_3(self):
        content = """
  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
    1   7   1   46.37 -161.29  134.95  0.376  0.149 -0.0971002.14  0.614  0.306
        """
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        mp = MolrepParser(molrep_log.name)
        self.assertEqual(mp.score, 0.306)
        self.assertEqual(mp.tfscore, 1002.14)
        self.assertEqual(mp.wrfac, 0.614)

    def test_molrep_parser_4(self):
        content = """
  Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
Lorem ipsum dolor sit amet ...
        """
        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()
        with self.assertRaises(ValueError):
            mp = MolrepParser(molrep_log.name)



if __name__ == "__main__":
    unittest.main()
