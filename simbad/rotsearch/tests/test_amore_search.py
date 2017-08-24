"""Test functions for rotsearch.amore_search.AmoreRotationSearch"""

__author__ = "Adam Simpkin"
__date__ = "16 Aug 2017"

import os
import unittest
from simbad.command_line import ccp4_root
from simbad.rotsearch import amore_search


class Test(unittest.TestCase):
    """Unit test"""
    
    def test_calculate_integration_box(self):
        """Test case for AmoreRotationSearch.calculate_integration_box"""
        
        input_model = os.path.join(ccp4_root(), "examples", "toxd", "toxd.pdb")
        data = amore_search.AmoreRotationSearch.calculate_integration_box(input_model)
        
        reference_data = (48.306749999999994, 56.73474999999999, 48.589749999999995, 19.84575)
        
        self.assertEqual(data, reference_data)
        
    def test_rotfun_1(self):
        """Test case for AmoreRotationSearch.rotfun"""
        amore_exe = "amore" 
        table1 = "table1"
        hklpck1 = "hklpck1"
        clmn1 = "clmn1"
        shres = "shres"
        intrad = "intrad"
        
        data_1, data_2 = amore_search.AmoreRotationSearch.rotfun(amore_exe, table1, hklpck1, clmn1, shres, intrad)
        
        reference_data_1 = ['amore', 'table1', 'table1', 'HKLPCK1', 'hklpck1', 'clmn1', 'clmn1']
        reference_data_2 = """ROTFUN
VERB
TITLE : Generate HKLPCK1 from MODEL FRAGMENT   1
GENE 1   RESO 100.0 shres  CELL_MODEL 80 75 65
CLMN MODEL 1     RESO  20.0  shres SPHERE   intrad"""

        self.assertEqual(data_1, reference_data_1)
        self.assertEqual(data_2, reference_data_2)
        
    def test_rotfun_2(self):
        """Test case for AmoreRotationSearch.rotfun"""
        
        amore_exe = "amore"
        table1 = "table1"
        hklpck1 = "hklpck1"
        clmn1 = "clmn1"
        hklpck0 = 'hklpck2'
        clmn0 = "clmn0"
        shres = "shres"
        intrad = "intrad"
        mapout = "mapout"
        pklim = "pklim"
        npic = "npic"
        step = "step"
            
        data_1, data_2 = amore_search.AmoreRotationSearch.rotfun(amore_exe, table1, hklpck1, clmn1, shres, intrad, 
                                                                 clmn0=clmn0, hklpck0=hklpck0, mapout=mapout, pklim=pklim,
                                                                 npic=npic, rotastep=step)
        
        reference_data_1 = ['amore', 'table1', 'table1', 'HKLPCK1', 'hklpck1', 'hklpck0', 'hklpck2', 'clmn1', 'clmn1', 'clmn0', 'clmn0', 'MAPOUT', 'mapout']
        reference_data_2 = """ROTFUN
TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1
GENE 1   RESO 100.0 shres  CELL_MODEL 80 75 65
CLMN CRYSTAL ORTH  1 RESO  20.0  shres  SPHERE   intrad
CLMN MODEL 1     RESO  20.0  shres SPHERE   intrad
ROTA  CROSS  MODEL 1  PKLIM pklim  NPIC npic STEP step"""

        self.assertEqual(data_1, reference_data_1)
        self.assertEqual(data_2, reference_data_2)
        
    def test_solvent_content(self):
        """Test case for AmoreRotationSearch.solvent_content"""
        
        input_model = os.path.join(ccp4_root(), "examples", "toxd", "toxd.pdb")
        unit_cell = '73.58 38.73 23.19 90.00 90.00 90.00'
        space_group = 'P212121'
        data = amore_search.AmoreRotationSearch.solvent_content(input_model, unit_cell, space_group)
        
        reference_data = 46.82229046138755
        
        self.assertEqual(data, reference_data)
        
    def test_tabfun_1(self):
        """Test case for AmoreRotationSearch.tabfun"""
        
        amore_exe = "amore"
        xyzin1 = "xyzin1"
        xyzout1 = "xyzout1"
        table1 = "table1"
        
        data_1, data_2 = amore_search.AmoreRotationSearch.tabfun(amore_exe, xyzin1, xyzout1, table1)
        
        reference_data_1 = ['amore', 'xyzin1', 'xyzin1', 'xyzout1', 'xyzout1', 'table1', 'table1']
        reference_data_2 = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL 200 200 200 90 90 120 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0"""
        
        self.assertEqual(data_1, reference_data_1)
        self.assertEqual(data_2, reference_data_2)
        
    def test_tabfun_2(self):
        """Test case for AmoreRotationSearch.tabfun"""
        
        amore_exe = "amore"
        xyzin1 = "xyzin1"
        xyzout1 = "xyzout1"
        table1 = "table1"
        x=100 
        y=200
        z=300
        a=10
        b=20
        c=30
        
        data_1, data_2 = amore_search.AmoreRotationSearch.tabfun(amore_exe, xyzin1, xyzout1, table1, x, y, z, a, b, c)
        
        reference_data_1 = ['amore', 'xyzin1', 'xyzin1', 'xyzout1', 'xyzout1', 'table1', 'table1']
        reference_data_2 = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL 100 200 300 10 20 30 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0"""
        
        self.assertEqual(data_1, reference_data_1)
        self.assertEqual(data_2, reference_data_2)
        
if __name__ == "__main__":
    unittest.main()