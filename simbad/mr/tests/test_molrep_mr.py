"""Test functions for simbad.mr.molrep_mr"""

__author__ = "Adam Simpkin"
__date__ = "20 Feb 2018"

import os
import tempfile
import unittest
import simbad
import simbad.mr.molrep_mr


class Test(unittest.TestCase):
    """Unit test"""

    def test_check_contrast_1(self):
        """Test case for simbad.mr.molrep_mr.check_contrast"""

        content = """
        corrF =   0.6346
        TF/sig       =    20.29
        Final CC     =   0.6346
        Packing_Coef =   1.0000
        Contrast     =    13.93
        CC_for_fixed_model:  0.6320

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         10   2   1  124.85  134.13  119.00  0.780 -0.786  0.043  20.29  0.546  0.635
        """

        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()

        data = simbad.mr.molrep_mr.check_contrast(molrep_log.name)
        os.unlink(molrep_log.name)
        reference_data = 13.93

        self.assertEqual(data, reference_data)

    def test_check_contrast_2(self):
        """Test case for simbad.mr.molrep_mr.check_contrast"""

        content = """
        corrF =   0.5306
        TF/sig       =    11.12
        Final CC     =   0.5306
        Packing_Coef =   1.0000
        Contrast is not available

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         2   3   1   54.32  134.28  121.04  0.004 -0.183  0.176  11.12  0.580  0.531

        corrF =   0.5641
        TF/sig       =    13.38
        Final CC     =   0.5641
        Packing_Coef =   1.0000
        Contrast is not available
        CC_for_fixed_model:  0.5306

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         4   1   1   90.00   -0.90  180.00  0.338 -0.424  0.221  13.38  0.561  0.564

        corrF =   0.6012
        TF/sig       =    14.44
        Final CC     =   0.6012
        Packing_Coef =   1.0000
        Contrast is not available
        CC_for_fixed_model:  0.5641

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         6   2   1  124.85  134.13  119.00  0.280 -0.282  0.041  14.44  0.543  0.601

        corrF =   0.6320
        TF/sig       =    22.52
        Final CC     =   0.6320
        Packing_Coef =   1.0000
        Contrast is not available
        CC_for_fixed_model:  0.6012

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         8   1   1  180.00    0.00  178.20  0.162 -0.425  0.274  22.52  0.538  0.632

        corrF =   0.6346
        TF/sig       =    20.29
        Final CC     =   0.6346
        Packing_Coef =   1.0000
        Contrast     =    13.93
        CC_for_fixed_model:  0.6320

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         10   2   1  124.85  134.13  119.00  0.780 -0.786  0.043  20.29  0.546  0.635

        corrF =   0.6201
        TF/sig       =     4.15
        Final CC     =   0.6201
        Packing_Coef =   1.0000
        Contrast     =     1.93
        CC_for_fixed_model:  0.6346

        WARNING: program can not improve current model
                result is "molrep.crd" with    10 monomers
                Try to use dimer
        """

        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()

        data = simbad.mr.molrep_mr.check_contrast(molrep_log.name)
        os.unlink(molrep_log.name)
        reference_data = 13.93

        self.assertEqual(data, reference_data)

    def test_check_contrast_3(self):
        """Test case for simbad.mr.molrep_mr.check_contrast"""

        content = """
        corrF =   0.5306
        TF/sig       =    11.12
        Final CC     =   0.5306
        Packing_Coef =   1.0000
        Contrast is not available

        Nmon RF  TF   theta    phi     chi   tx     ty     tz     TF/sg  wRfac  Score
         2   3   1   54.32  134.28  121.04  0.004 -0.183  0.176  11.12  0.580  0.531
        """

        molrep_log = tempfile.NamedTemporaryFile("w", delete=False)
        molrep_log.write(content)
        molrep_log.close()

        data = simbad.mr.molrep_mr.check_contrast(molrep_log.name)
        os.unlink(molrep_log.name)
        reference_data = 0.0

        self.assertEqual(data, reference_data)


if __name__ == "__main__":
    unittest.main()
