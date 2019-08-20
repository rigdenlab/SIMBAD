"""Test functions for simbad.util.pdb_util"""

__author__ = "Adam Simpkin"
__date__ = "19 Jan 2018"

import os
import tempfile
import unittest
import simbad.util


class Test(unittest.TestCase):
    """Unit test"""

    def test_result_by_score_from_csv_1(self):
        """Test case for simbad.util.result_by_score_from_csv"""

        csv_temp_file = tempfile.NamedTemporaryFile("w", delete=False)
        csv_temp_file.write(
            """pdb_code,alt,a,b,c,alpha,beta,gamma,length_penalty,angle_penalty,total_penalty,volume_difference,probability_score
1DTX, ,23.15,39.06,73.53,90.0,90.0,90.0,0.418,0.0,0.418,398.847,0.842"""
        )
        csv_temp_file.close()

        data = simbad.util.result_by_score_from_csv(csv_temp_file.name, "total_penalty")
        reference_data = ["1DTX", 0.41799999999999998]

        self.assertEqual(data, reference_data)

    def test_result_by_score_from_csv_2(self):
        """Test case for simbad.util.result_by_score_from_csv"""

        csv_temp_file = tempfile.NamedTemporaryFile("w", delete=False)
        csv_temp_file.write(
            """pdb_code,ALPHA,BETA,GAMMA,CC_F,RF_F,CC_I,CC_P,Icp,CC_F_Z_score,CC_P_Z_score,Number_of_rotation_searches_producing_peak
2fbb,21.63,81.88,296.6,14.1,56.2,16.5,18.6,1.0,11.6,8.6,5.0
1f10,34.27,90.0,116.04,13.0,57.1,16.4,14.2,1.0,9.0,7.0,5.0
4w94,29.28,85.42,245.3,12.9,57.2,15.2,10.8,1.0,8.9,7.1,5.0
1xei,38.87,78.75,65.8,12.3,58.0,15.4,13.9,1.0,8.1,6.6,5.0
2z18,27.6,87.35,247.57,12.3,57.5,15.3,12.5,1.0,7.8,6.1,5.0
1ps5,33.92,86.37,67.25,12.6,57.3,15.6,14.8,1.0,7.7,7.4,5.0
1v7s,34.18,87.8,66.84,12.5,57.4,15.7,12.6,1.0,7.6,6.7,5.0
2vb1,37.1,85.56,66.78,12.1,57.3,16.2,12.3,1.0,7.6,6.6,5.0
4yeo,35.02,82.52,67.02,11.8,57.2,15.5,13.8,1.0,7.6,6.7,5.0
2b5z,1.4,38.12,229.38,12.4,57.9,15.4,10.4,1.0,7.6,6.5,5.0
1ykz,26.43,88.72,247.05,12.6,57.5,15.4,11.9,1.0,7.6,6.5,5.0
4xjf,26.78,88.44,245.77,12.9,57.8,15.4,12.7,1.0,7.6,6.5,5.0
2d4j,37.18,84.17,66.64,12.4,57.7,16.1,12.8,1.0,7.5,6.0,5.0
4p2e,29.05,83.8,246.58,12.5,56.9,15.4,12.1,1.0,7.5,7.1,5.0
3wvx,35.67,85.1,67.1,12.6,57.1,15.1,13.0,1.0,7.4,6.4,5.0
2x0a,28.59,85.11,245.89,12.3,57.4,14.8,11.4,1.0,7.4,6.5,5.0
2z19,38.05,79.03,64.98,11.8,57.7,15.9,12.8,1.0,7.1,5.9,5.0
1jj1,28.99,82.92,245.93,12.5,57.3,15.6,11.3,1.0,7.0,5.9,5.0
4j7v,28.54,86.74,246.59,12.0,57.4,14.2,10.2,1.0,7.0,5.7,5.0
2pc2,28.71,76.6,257.7,10.8,57.7,13.4,8.0,1.0,6.7,5.3,5.0"""
        )
        csv_temp_file.close()

        data = simbad.util.result_by_score_from_csv(csv_temp_file.name, "CC_F_Z_score")
        reference_data = ["2fbb", 11.6]

        self.assertEqual(data, reference_data)
