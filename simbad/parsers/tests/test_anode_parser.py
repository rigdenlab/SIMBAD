"""Test functions for simbad.parsers.anode_parser"""

__author__ = "Adam Simpkin"
__date__ = "17 Aug 2017"

import tempfile
import unittest

from simbad.parsers import anode_parser


class Test(unittest.TestCase):
    def test_anode_parser(self):
        content = """ANODE - ANOmalous DEnsity analysis - version 2013/1
 ===================================================

 Command line: anode -b10.0 -d1.0 -h80 -m20 -r5.0 -n99.0 -s4.0 4WKI

  16 threads running in parallel on 16 CPUs

  1909 Atoms read from PDB file 4WKI.pdb
 Cell:   40.037   68.193   74.013  90.00  90.00  90.00
 Space group: P 21 21 21

   16870 Reflections read from file 4WKI_fa.hkl
         Highest resolution =   1.600 A


 Averaged anomalous densities (sigma)

  39.51  ZN_ZN
   7.71  CA_CA
   3.23  SD_MET
   2.53  SG_CYS
   1.69  CL2_3PW
   1.60  N18_3PW
   1.49  C4_3PW
   1.29  C9_3PW
   1.09  CB_TRP
   1.08  OH_TYR
   0.91  O25_3PW
   0.89  CZ_TYR
   0.69  N_TRP
   0.68  CA_TRP
   0.68  O_ARG
   0.61  C11_3PW
   0.61  C12_3PW
   0.51  C3_3PW
   0.50  O_ILE
   0.44  CH2_TRP
   
    Strongest unique anomalous peaks

          X        Y        Z   Height(sig)  SOF     Nearest atom

 ZN1  -0.14568  0.10187  0.31976   39.59    1.000    0.049  ZN_A:ZN501
 ZN2  -0.26378 -0.06058  0.37985    9.30    1.000    0.149  CA_A:CA504
 ZN3   0.34895  0.18006  0.48898    7.78    1.000    0.138  CA_A:CA502
 ZN4   0.40560  0.16274  0.52466    6.57    1.000    0.279  CA_A:CA503
 ZN5  -0.00596 -0.12591  0.45283    5.44    1.000    0.295  SD_A:MET230
 ZN6   0.17604  0.01967  0.46853    4.33    1.000    2.242  CG2_A:VAL224
 ZN7   0.01378  0.08476  0.39220    4.32    1.000    1.685  C_A:GLU362
 ZN8  -0.07702  0.06676  0.45484    4.27    1.000    1.911  N_A:ILE347
 ZN9   0.32379  0.28692  0.35855    4.25    1.000    0.770  NH2_A:ARG218
 ZN10  0.17778 -0.17610  0.28304    4.19    1.000    3.146  O_A:HOH816
 ZN11 -0.03473 -0.08544  0.34792    4.15    1.000    1.888  OD2_A:ASP351
 ZN12 -0.19532 -0.08462  0.46678    4.11    1.000    1.144  CB_A:GLN319
 ZN13  0.33612  0.27935  0.34665    4.09    1.000    1.165  NH2_A:ARG218
 ZN14  0.07757  0.15593  0.28974    4.07    1.000    0.483  SD_A:MET369
 ZN15  0.38938  0.20639  0.37188    4.05    1.000    1.223  CA_A:ARG218

   15 Peaks output to file 4WKI_fa.res

    16870 Reflections written to file 4WKI.pha

 Listing file 4WKI.lsa"""

        anode_log = tempfile.NamedTemporaryFile("w", delete=False)
        anode_log.write(content)
        anode_log.close()

        ap = anode_parser.AnodeParser(anode_log.name)

        self.assertEqual(ap.x, '-0.14568')
        self.assertEqual(ap.y, '0.10187')
        self.assertEqual(ap.z, '0.31976')
        self.assertEqual(ap.peak_height, '39.59')
        self.assertEqual(ap.nearest_atom, 'ZN_A:ZN501')


if __name__ == "__main__":
    unittest.main()
