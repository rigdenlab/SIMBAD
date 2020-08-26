"""Test functions for simbad.util.matthews_prob"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import numpy as np
import unittest
from simbad.util import matthews_prob

try:
    ROOT_DIR = os.environ['SIMBAD_ROOT']
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "test_data")
except KeyError:
    from simbad.command_line import CCP4RootDirectory
    ROOT_DIR = str(CCP4RootDirectory())
    EXAMPLE_DIR = os.path.join(ROOT_DIR, "examples")


class Test(unittest.TestCase):
    """Unit test"""

    def test_solvent_content(self):
        """Test case for matthews_prob.SolventContent.calculate_from_file"""

        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        volume = 16522.4616729
        SC = matthews_prob.SolventContent(volume)
        data = SC.calculate_from_file(input_model)

        reference_data = 0.48960068050640637

        self.assertAlmostEqual(np.round(data, 3), np.round(reference_data, 3))

    def test_matthews_prob(self):
        """Test case for matthews_prob.MatthewsProbability.calculate_from_file"""

        input_model = os.path.join(EXAMPLE_DIR, "toxd", "toxd.pdb")
        volume = 16522.4616729
        MC = matthews_prob.MatthewsProbability(volume)
        data = MC.calculate_from_file(input_model)

        reference_data = (0.48960068050640637, 1)

        self.assertAlmostEqual(np.round(data[0], 3), np.round(reference_data[0], 3))
        self.assertAlmostEqual(data[1], reference_data[1])


if __name__ == "__main__":
    unittest.main()
