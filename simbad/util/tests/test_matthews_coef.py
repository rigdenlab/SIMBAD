"""Test functions for rotsearch.amore_search.AmoreRotationSearch"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.command_line import ccp4_root
from simbad.util import matthews_coef


class Test(unittest.TestCase):
    """Unit test"""

    def test_solvent_content(self):
        """Test case for matthews_coef.solvent_content"""

        input_model = os.path.join(ccp4_root(), "examples", "toxd", "toxd.pdb")
        unit_cell = '73.58 38.73 23.19 90.00 90.00 90.00'
        space_group = 'P212121'
        data = matthews_coef.solvent_content(input_model, unit_cell, space_group)

        reference_data = 46.82229046138755

        self.assertEqual(data, reference_data)

    def test_matthews_coef(self):
        """Test case for matthews_coef.matthews_coef"""

        unit_cell = '73.58 38.73 23.19 90.00 90.00 90.00'
        space_group = 'P212121'
        nres = 25

        data = matthews_coef.matthews_coef(unit_cell, space_group, nres)

        reference_data = (0.5814862630712596, 2)

        self.assertEqual(data, reference_data)

if __name__ == "__main__":
    unittest.main()