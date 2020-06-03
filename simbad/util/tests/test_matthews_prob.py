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
        volume = 16522.4616729
        SC = matthews_prob.SolventContent(volume)
        data = SC.calculate_from_file(input_model)

        reference_data = 46.79124862303431

        self.assertAlmostEqual(data, reference_data)

    def test_matthews_prob(self):
        """Test case for matthews_prob.MatthewsProbability.calculate_from_file"""

        input_model = os.path.join(CCP4ROOT, "examples", "toxd", "toxd.pdb")
        volume = 16522.4616729
        MC = matthews_prob.MatthewsProbability(volume)
        data = MC.calculate_from_file(input_model)

        reference_data = (0.4679124862303431, 1)

        self.assertAlmostEqual(data[0], reference_data[0])
        self.assertAlmostEqual(data[1], reference_data[1])


if __name__ == "__main__":
    unittest.main()
