
##############################################################################
#
# Copyright (c) 2012-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2012-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
import esys.ripley
from esys.downunder import *
from esys.escript import *
import numpy as np

                
class TestRegularizaton2D(unittest.TestCase):
    def setUp(self):
        self.domain = esys.ripley.Rectangle(20*getMPISizeWorld()-1,20,d0=getMPISizeWorld()) # expected dimension 1mx1m
        self.COORDINATES=CartesianReferenceSystem()
        
    def tearDown(self):
        del self.domain
        
    # standart inversion case   
    def test_ConstantLevelSet1(self): # doesn't test the regularization

        reg=Regularization(dom, numLevelSets=1,
                               w1=[1,1.], # consider gradient terms
                               #wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        cf=InversionCostFunction(reg, pm1, sm1)
        
        
        m=Scalar(1., Solution(dom))
        args=cf.getArguments(m)
        df=cf.getValue(m, *args)
        gf=cf.getGradient(m, *args)[0]
        
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df, sm1.getRealValue(s0))
        
        self.assertIsInstance(gf, Data)
        self.assertEqual(gf.getFunctionSpace(), Function(dom))
        self.assertEqual(gf.getShape(), ())
        
        
        STEP=0.001
        m2=m+STEP
        df2=cf.getValue(m2, *(cf.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( Lsup(gf2-gf*STEP) < 1e-3 * Lsup(gf2) )
    
    def test_ConstantLevelSet2(self): # doesn't test the regularization
        reg=Regularization(dom, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,0],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        
        
        m=Data([1.,2], Solution(dom))
        args=cf.getArguments(m)
        
        df=cf.getValue(m, *args)
        df_real=sm1.getRealValue(s[0]*1.)+sm2.getRealValue(s[1]*2.)
        gf=cf.getGradient(m, *args)[0]
    
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df, df_real)
        
        self.assertIsInstance(gf, Data)
        self.assertEqual(gf.getFunctionSpace(), Function(dom))
        self.assertEqual(gf.getShape(), (2,))


        STEP=0.001

        m2=m+[STEP,0]
        df2=cf.getValue(m2, *(cf.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( Lsup(gf2-gf[0]*STEP) < 1e-3 * Lsup(gf2) )

        m2=m+[0, STEP]
        df2=cf.getValue(m2, *(cf.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( Lsup(gf2-gf[1]*STEP) < 1e-3 * Lsup(gf2) )

        
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

