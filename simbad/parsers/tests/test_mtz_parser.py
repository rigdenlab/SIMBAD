import os
import numpy as np
import unittest

from simbad.parsers.mtz_parser import MtzParser

try:
    ROOT_DIR = os.environ['SIMBAD_ROOT']
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "test_data")
except KeyError:
    from simbad.command_line import CCP4RootDirectory
    ROOT_DIR = str(CCP4RootDirectory())
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "examples")


class Test(unittest.TestCase):
    """Unit test"""
    def test_mtz_parser_1(self):
        """Test case for MtzParser"""
        input_mtz = os.path.join(EXAMPLE_DIR, "toxd", "toxd.mtz")
        mp = MtzParser(input_mtz)
        mp.parse()
        self.assertEqual(np.round(mp.resolution, 1), 2.3)
        self.assertEqual(np.round(mp.cell.a, 4), 73.5820)
        self.assertEqual(np.round(mp.cell.b, 4), 38.7330)
        self.assertEqual(np.round(mp.cell.c, 4), 23.1890)
        self.assertEqual(np.round(mp.cell.alpha, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.beta, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.gamma, 4), 90.0000)
        self.assertEqual(mp.spacegroup_symbol, "P 21 21 21")
        self.assertEqual(mp.nreflections, 3235)
        self.assertEqual(mp.f, "FTOXD3")
        self.assertEqual(mp.sigf, "SIGFTOXD3")
        self.assertEqual(mp.free, "FreeR_flag")

    def test_mtz_parser_2(self):
        """Test case for MtzParser"""
        input_mtz = os.path.join(EXAMPLE_DIR, "rnase", "rnase25.mtz")
        mp = MtzParser(input_mtz)
        mp.parse()
        self.assertEqual(np.round(mp.resolution, 1), 2.5)
        self.assertEqual(np.round(mp.cell.a, 4), 64.8970)
        self.assertEqual(np.round(mp.cell.b, 4), 78.3230)
        self.assertEqual(np.round(mp.cell.c, 4), 38.7920)
        self.assertEqual(np.round(mp.cell.alpha, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.beta, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.gamma, 4), 90.0000)
        self.assertEqual(mp.spacegroup_symbol, "P 21 21 21")
        self.assertEqual(mp.nreflections, 7262)
        self.assertEqual(mp.f, "FNAT")
        self.assertEqual(mp.sigf, "SIGFNAT")
        self.assertEqual(mp.free, "FreeR_flag")

    def test_mtz_parser_3(self):
        """Test case for MtzParser"""
        input_mtz = os.path.join(EXAMPLE_DIR, "rnase", "rnase25F+F-.mtz")
        mp = MtzParser(input_mtz)
        mp.parse()
        self.assertEqual(np.round(mp.resolution, 1), 2.5)
        self.assertEqual(np.round(mp.cell.a, 4), 64.8970)
        self.assertEqual(np.round(mp.cell.b, 4), 78.3230)
        self.assertEqual(np.round(mp.cell.c, 4), 38.7920)
        self.assertEqual(np.round(mp.cell.alpha, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.beta, 4), 90.0000)
        self.assertEqual(np.round(mp.cell.gamma, 4), 90.0000)
        self.assertEqual(mp.spacegroup_symbol, "P 21 21 21")
        self.assertEqual(mp.nreflections, 7262)
        self.assertEqual(mp.f, "FNAT")
        self.assertEqual(mp.sigf, "SIGFNAT")
        self.assertEqual(mp.f_plus, "FPTNCD25(+)")
        self.assertEqual(mp.sigf_plus, "SIGFPTNCD25(+)")
        self.assertEqual(mp.f_minus, "FPTNCD25(-)")
        self.assertEqual(mp.sigf_minus, "SIGFPTNCD25(-)")
        self.assertEqual(mp.free, "FreeR_flag")



