"""Test functions for rotsearch.AmoreRotationSearch and rotsearch.PhaserRotationSearch"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import unittest
import simbad.rotsearch.amore_search
import simbad.rotsearch.phaser_search


class Test(unittest.TestCase):
    """Unit test"""

    @classmethod
    def setUpClass(cls):
        cls.AS = simbad.rotsearch.amore_search.AmoreRotationSearch("mtz", "molrep", "tmp_dir", "work_dir")
        cls.PS = simbad.rotsearch.phaser_search.PhaserRotationSearch("mtz", "molrep", "tmp_dir", "work_dir")

    def test_sortfun(self):
        """Test case for AmoreRotationSearch.sortfun_stdin_template"""

        f = "f"
        sigf = "sigf"

        data = self.AS.sortfun_stdin_template.format(f=f, sigf=sigf)

        reference_data = """TITLE   ** spmi  packing h k l F for crystal**
SORTFUN RESOL 100.  2.5
LABI FP=f  SIGFP=sigf"""

        self.assertEqual(data, reference_data)

    def test_rotfun(self):
        """Test case for AmoreRotationSearch.rotfun_stdin_template"""
        shres = "shres"
        intrad = "intrad"
        pklim = "pklim"
        npic = "npic"
        step = "step"

        data = self.AS.rotfun_stdin_template.format(shres=shres, intrad=intrad, pklim=pklim, npic=npic, step=step)

        reference_data = """TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1
ROTFUN
GENE 1   RESO 100.0 shres  CELL_MODEL 80 75 65
CLMN CRYSTAL ORTH  1 RESO  20.0  shres  SPHERE   intrad
CLMN MODEL 1     RESO  20.0  shres SPHERE   intrad
ROTA  CROSS  MODEL 1  PKLIM pklim  NPIC npic STEP step"""

        self.assertEqual(data, reference_data)

    def test_tabfun(self):
        """Test case for AmoreRotationSearch.tabfun_stdin_template"""
        x = 100
        y = 200
        z = 300
        a = 10
        b = 20
        c = 30

        data = self.AS.tabfun_stdin_template.format(x=x, y=y, z=z, a=a, b=b, c=c)

        reference_data = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL 100 200 300 10 20 30 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0"""

        self.assertEqual(data, reference_data)

    def test_rot_job_succeeded_1(self):
        """Test case for AmoreRotationSearch._rot_job_succeeded"""
        amore_z_score = 15
        data = self.AS._rot_job_succeeded(amore_z_score)

        self.assertTrue(data)

    def test_rot_job_succeeded_2(self):
        """Test case for AmoreRotationSearch._rot_job_succeeded"""
        amore_z_score = 8
        data = self.AS._rot_job_succeeded(amore_z_score)

        self.assertFalse(data)

    def test_rot_job_succeeded_3(self):
        """Test case for PhaserRotationSearch._rot_job_succeeded"""
        phaser_z_score = 15
        data = self.PS._rot_job_succeeded(phaser_z_score)

        self.assertTrue(data)

    def test_rot_job_succeeded_4(self):
        """Test case for PhaserRotationSearch._rot_job_succeeded"""
        phaser_z_score = 6
        data = self.PS._rot_job_succeeded(phaser_z_score)

        self.assertFalse(data)


if __name__ == "__main__":
    unittest.main()
