"""Test functions for simbad.util.matthews_prob"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.command_line import CCP4RootDirectory
from simbad.util import matthews_prob

CCP4ROOT = str(CCP4RootDirectory())


class Test(unittest.TestCase):
    """Unit test"""

    def test_solvent_content(self):
        """Test case for matthews_prob.SolventContent.calculate_from_file"""

        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        unit_cell = "73.58 38.73 23.19 90.00 90.00 90.00"
        space_group = "P212121"
        SC = matthews_prob.SolventContent(unit_cell, space_group)
        data = SC.calculate_from_file(input_model)

        reference_data = 46.82229046138755

        self.assertEqual(data, reference_data)

    def test_matthews_prob(self):
        """Test case for matthews_prob.MatthewsProbability.calculate_content_ncopies_from_file"""

        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        unit_cell = "73.58 38.73 23.19 90.00 90.00 90.00"
        space_group = "P212121"
        MC = matthews_prob.MatthewsProbability(unit_cell, space_group)
        data = MC.calculate_content_ncopies_from_file(input_model)

        reference_data = (0.5061537904240863, 1)

        self.assertEqual(data, reference_data)


if __name__ == "__main__":
    unittest.main()
