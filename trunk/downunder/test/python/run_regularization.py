##############################################################################
#
# Copyright (c) 2012-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2012-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
from esys.downunder import *
from esys.escript import *
import numpy as np

try:
    import esys.ripley
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class Test_Regularizaton2D(unittest.TestCase):
    def setUp(self):
        self.domain = esys.ripley.Rectangle(20*getMPISizeWorld()-1,20*getMPISizeWorld()-1, d0=getMPISizeWorld()) # expected dimension 1mx1m
        self.COORDINATES=CartesianReferenceSystem()
        
    def tearDown(self):
        del self.domain
         
    def test_ConstantLevelSet1(self): 
        a=[1,2]
        
        reg=Regularization(self.domain, numLevelSets=1,
                               w1=[1,1.], # consider gradient terms
                               #wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        
        
        m=Scalar(1., Solution(self.domain))
        args=reg.getArguments(m)
        df=reg.getValue(m, *args)
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df,0.)
        
        gf=reg.getGradient(m, *args)
        self.assertIsInstance(gf[0], Data)
        self.assertIsInstance(gf[1], Data)
        self.assertEqual(gf[0].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[1].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[0].getShape(), ())
        self.assertEqual(gf[1].getShape(), (2,))
        self.assertAlmostEqual(Lsup(gf[0]),0.)
        self.assertAlmostEqual(Lsup(gf[1]),0.)
        
        # lets try another m
        x=Solution(self.domain).getX()
        m=a[0]*x[0]+a[1]*x[1]

        args=reg.getArguments(m)
        df=reg.getValue(m, *args)
        self.assertAlmostEqual(df,(a[0]**2+a[1]**2)/2./2.)
        gf=reg.getGradient(m, *args)
        # and now the derivatives:                
        STEP=0.001
        dm=STEP*x[0]**2
        m2=m+dm
        df2=reg.getValue(m2, *(reg.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( abs(gf2-reg.getDualProduct(dm, gf)) < 1e-3 * abs(gf2) )

    def test_ConstantLevelSet1Hessian(self): 
        x=Solution(self.domain).getX()
        
        
        
        reg=Regularization(self.domain, numLevelSets=1,
                               w1=[1,1.], # consider gradient terms
                               #wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES,
                               location_of_set_m=whereZero(x[0]-inf(x[0])) )
        
        x=Solution(self.domain).getX()
        m=x[0]
        
        gf=reg.getGradient(m, *(reg.getArguments(m)))
        
        STEP=0.001
        dm=STEP*x[0]*(x[0]-1)
        m2=m+dm
        gf2=reg.getGradient(m2, *(reg.getArguments(m2)))

        p=reg.getInverseHessianApproximation(m, gf2-gf, *(reg.getArguments(m)))
        self.assertAlmostEqual(Lsup(p-dm)/Lsup(dm),0.)

    
    def test_ConstantLevelSet2_noCrossGradient(self): 
        a=[1,2]
                
        reg=Regularization(self.domain, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,0],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        
        
        m=Data([1.,2], Solution(self.domain))
        args=reg.getArguments(m)
        
        df=reg.getValue(m, *args)
        
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df,0.)
        
        gf=reg.getGradient(m, *args)
        self.assertIsInstance(gf[0], Data)
        self.assertIsInstance(gf[1], Data)
        self.assertEqual(gf[0].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[1].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[0].getShape(), (2,))
        self.assertEqual(gf[1].getShape(), (2,2))
        self.assertAlmostEqual(Lsup(gf[0]),0.)
        self.assertAlmostEqual(Lsup(gf[1]),0.)

        # lets try another m
        x=Solution(self.domain).getX()
        m=x*a

        args=reg.getArguments(m)
        df=reg.getValue(m, *args)
        self.assertAlmostEqual(df,(a[0]**2+a[1]**2)/2./2.)
        gf=reg.getGradient(m, *args)
        
        # and now the derivatives:                
        STEP=0.0002
        dm=STEP*x[0]*x[1]
        
        m2=m+dm*[1,0]
        df2=reg.getValue(m2, *(reg.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( abs(gf2-reg.getDualProduct(dm*[1,0], gf)) < 1e-3 * abs(gf2) )

        m2=m+dm*[0,1]
        df2=reg.getValue(m2, *(reg.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( abs(gf2-reg.getDualProduct(dm*[0,1], gf)) < 1e-3 * abs(gf2) )

    def test_ConstantLevelSet2_noCrossGradientHessian(self):
        x=Solution(self.domain).getX()
        
        reg=Regularization(self.domain, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,0],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES,
                               location_of_set_m=whereZero(x[0]-inf(x[0]))*[1,1] )
        
        m=x[0]*[1,0]+x[1]*x[1]*[0,1]
        
        gf=reg.getGradient(m, *(reg.getArguments(m)))
        STEP=0.001
        dm=STEP*x[0]*(x[0]-1)
        
        m2=m+dm*[1,0]
        gf2=reg.getGradient(m2, *(reg.getArguments(m2)))
        p=reg.getInverseHessianApproximation(m, gf2-gf, *(reg.getArguments(m)))
        self.assertAlmostEqual(Lsup(p[0]-dm)/Lsup(dm),0.)
        self.assertAlmostEqual(Lsup(p[1])/Lsup(dm),0.)

        m2=m+dm*[0,1]
        gf2=reg.getGradient(m2, *(reg.getArguments(m2)))
        p=reg.getInverseHessianApproximation(m, gf2-gf, *(reg.getArguments(m)))
        self.assertAlmostEqual(Lsup(p[1]-dm)/Lsup(dm),0.)
        self.assertAlmostEqual(Lsup(p[0])/Lsup(dm),0.)

        
    def test_ConstantLevelSet2_WithCrossGradient(self): # doesn't test the regularization
        reg=Regularization(self.domain, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        reg.setTradeOffFactorsForVariation([1.e-8,1.e-8])
        
        m=Data([1.,2], Solution(self.domain))
        args=reg.getArguments(m)
        
        df=reg.getValue(m, *args)
        
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df,0.)
        
        gf=reg.getGradient(m, *args)
        self.assertIsInstance(gf[0], Data)
        self.assertIsInstance(gf[1], Data)
        self.assertEqual(gf[0].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[1].getFunctionSpace(), Function(self.domain))
        self.assertEqual(gf[0].getShape(), (2,))
        self.assertEqual(gf[1].getShape(), (2,2))
        self.assertAlmostEqual(Lsup(gf[0]),0.)
        self.assertAlmostEqual(Lsup(gf[1]),0.)
        
        # lets try another m: 
        x=Solution(self.domain).getX()
        a1=1.
        a2=2.
        # for this one cross gradient should not impact on cost function
        f=0.1
        m=(a1*x[0]+a2*x[1]) * [1,0] + (f*a1*x[0]+f*a2*x[1]) * [0,1] # cross gradient term is zero!
        args=reg.getArguments(m)
        df=reg.getValue(m, *args)
        #self.assertAlmostEqual(df,(a1**2+a2**2)*(1+f**2)/2./2.)
        self.assertAlmostEqual(df,0.)
        gf=reg.getGradient(m, *args)
        self.assertAlmostEqual(Lsup(gf[0]),0.)
        self.assertAlmostEqual(Lsup(gf[1]),0.)
        
        # for this gives maximum impact on  on cost function 
        m=(a1*x[0]+a2*x[1]) * [1,0] + (a2*x[0]-a1*x[1]) * [0,1] # cross gradient term is zero!
        args=reg.getArguments(m)
        df=reg.getValue(m, *args)
        self.assertAlmostEqual(df,((a1**2+a2**2)**2/4.)/2.)
        gf=reg.getGradient(m, *args)
        # and now the derivatives:                
        STEP=0.002
        dm=STEP*x[0]*x[1]
        
        m2=m+dm*[1,0]
        df2=reg.getValue(m2, *(reg.getArguments(m2)))
        gf2=df2-df      
        self.assertTrue( abs(gf2-reg.getDualProduct(dm*[1,0], gf)) < 1e-3 * abs(gf2) )

        m2=m+dm*[0,1]
        df2=reg.getValue(m2, *(reg.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( abs(gf2-reg.getDualProduct(dm*[0,1], gf)) < 1e-3 * abs(gf2) )
        
    def test_ConstantLevelSet2_WithCrossGradientHessian(self):
        x=Solution(self.domain).getX()
        
        reg=Regularization(self.domain, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES,
                               location_of_set_m=whereZero(x[0]-inf(x[0]))*[1,1] )
        reg.setTradeOffFactorsForVariation([1e-10,1e-10])
        m=x[0]*[1,0]+x[1]*x[1]*[0,1]
        
        gf=reg.getGradient(m, *(reg.getArguments(m)))
        STEP=0.001
        dm=STEP*x[0]*x[1]

        for vector in [[1,0], [0,1]]:
            m2=m+dm*vector
            gf2=reg.getGradient(m2, *(reg.getArguments(m2)))

            p=reg.getInverseHessianApproximation(m, gf2-gf, *(reg.getArguments(m)),
                    solve=False)
            X = (gf2-gf)[1]
            A = p.getCoefficient("A")
            res = generalTensorProduct(A, grad(dm*vector), axis_offset=2)
            self.assertLessEqual(Lsup(res-X), 5e-3 * Lsup(X))
        
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

