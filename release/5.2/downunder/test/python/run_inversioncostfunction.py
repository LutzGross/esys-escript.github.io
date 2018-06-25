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

class LinearMappingX(Mapping):
    def __init__(self, mat):
        self.matrix=mat
        self.det=mat[0,0]*mat[1,1]-mat[0,1]*mat[0,1]
        self.inv=np.array([ [mat[1,1]/self.det, -mat[1,0]/self.det ], [-mat[0,1]/self.det, mat[0,0]/self.det]   ])
    
    def getDerivative(self, m):
        return Data(self.matrix, m.getFunctionSpace())
    def getValue(self, p):
        return matrix_mult(self.matrix, p)
    def getInverse(self, p):
        return matrix_mult(self.inv, p)
    def getTypicalDerivative(self):
        return self.matrix[1,1]

class SimpleModel(ForwardModel):
    def __init__(self, domain, coordinates, numComps=1):
        self.domain = domain
        self.trafo=makeTransformation(domain, coordinates)
        self.numComps=numComps

    def getCoordinateTransformation(self):
        return self.trafo 
        
    def getArguments(self, x):
        return tuple([ n+1 for n in range(self.numComps)] )

    def getDefect(self, x, *args):
        if self.numComps == 1:
            out=(1.*x-2.)*(args[0]*x-2.)
        else:
            out=0
            for n in range(self.numComps):
                out=out+((n+1)*x[n]-(n+2))*(args[n]*x[n]-(n+2))  
        return integrate(out, Function(self.domain))/2

    def getGradient(self, x, *args):
        if self.numComps == 1:
            Y=(1.*x-2.)*args[0]
        else:
            Y=Data(0.,(self.numComps,), Function(self.domain))
            for n in range(self.numComps):
                Y[n]=((n+1)*x[n]-(n+2))*args[n]
        return Y
        
    def getRealValue(self, s):
        if self.numComps == 1:
            return (s-2)**2 * 0.5
        else:
            return sum(((n+1)*s[n]-(n+2))**2 for n in range(self.numComps) )*.5    

@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class TestInversionCostfunction(unittest.TestCase):
    def setUp(self):
        self.domain = esys.ripley.Rectangle(20*getMPISizeWorld()-1,20,d0=getMPISizeWorld()) # expected dimension 1mx1m
        self.COORDINATES=CartesianReferenceSystem()
        
    def tearDown(self):
        del self.domain
        
    # standart inversion case   
    def test_ConstantLevelSet1_Mapping1_Model1(self): # doesn't test the regularization
        dom=self.domain
        s0=10.
        
        pm1=LinearMapping(s0, 0.) #p = a * m + p0
        sm1=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)
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
    
    # strong magnetic-gravity coupling    
    def test_ConstantLevelSet1_Mapping2_Model1(self): # doesn't test the regularization
        dom=self.domain
        s=[10., 20]
        
        pm0=LinearMapping(s[0], 0.) #p = a * m + p0
        pm1=LinearMapping(s[1], 0.) 
        
        sm1=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)
        sm2=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)

        reg=Regularization(dom, numLevelSets=1,
                               w1=[1,1.], # consider gradient terms
                               #wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        cf=InversionCostFunction(reg, [ pm0, pm1 ], [ (sm1,0), (sm2,1) ] )
        
        
        m=Scalar(1., Solution(dom))
        args=cf.getArguments(m)
        
        df=cf.getValue(m, *args)
        df_real=sm1.getRealValue(s[0])+sm2.getRealValue(s[1])
        gf=cf.getGradient(m, *args)[0]
    
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df, df_real)
        
        self.assertIsInstance(gf, Data)
        self.assertEqual(gf.getFunctionSpace(), Function(dom))
        self.assertEqual(gf.getShape(), ())

        STEP=0.001
        m2=m+STEP
        df2=cf.getValue(m2, *(cf.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( Lsup(gf2-gf*STEP) < 1e-3 * Lsup(gf2) )
    
    # strong magnetic-gravity coupling
    def test_ConstantLevelSet1_Mapping2_Model1_V2(self): # doesn't test the regularization
        dom=self.domain
        s=[10., 20]
        
        pm0=LinearMapping(s[0], 0.) #p = a * m + p0
        pm1=LinearMapping(s[1], 0.)
        
        sm1=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)
        sm2=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)

        reg=Regularization(dom, numLevelSets=1,
                               w1=[1,1.], # consider gradient terms
                               #wc=[[0,1],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        cf=InversionCostFunction(reg, [ (pm0,0), (pm1,0) ], [ (sm1,0), (sm2,1) ] )
        
        
        m=Scalar(1., Solution(dom))
        args=cf.getArguments(m)
        
        df=cf.getValue(m, *args)
        df_real=sm1.getRealValue(s[0])+sm2.getRealValue(s[1])
        gf=cf.getGradient(m, *args)[0]

    
        self.assertIsInstance(df, float)
        self.assertAlmostEqual(df, df_real)
        
        self.assertIsInstance(gf, Data)
        self.assertEqual(gf.getFunctionSpace(), Function(dom))
        self.assertEqual(gf.getShape(), ())


        STEP=0.001
        m2=m+STEP
        df2=cf.getValue(m2, *(cf.getArguments(m2)))
        gf2=df2-df
        self.assertTrue( Lsup(gf2-gf*STEP) < 1e-3 * Lsup(gf2) )
    
    # weak gravity- magnetic coupling case:        
    def test_ConstantLevelSet2_Mapping2_Model2(self): # doesn't test the regularization
        dom=self.domain
        s=[15., 7.]
        
        pm0=LinearMapping(s[0], 0.) #p = a * m + p0
        pm1=LinearMapping(s[1], 0.) 
        
        sm1=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)
        sm2=SimpleModel(dom, numComps=1,  coordinates=self.COORDINATES)

        reg=Regularization(dom, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,0],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        cf=InversionCostFunction(reg, [ (pm0,0), (pm1,1) ], [ (sm1,0), (sm2,1) ] )
        
        
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

    #  acoustic Wave form inversion case 
    def test_ConstantLevelSet2_Mapping1X2_Model1(self): # doesn't test the regularization
        dom=self.domain
        s=np.array([ [15., 7.] , [-4, 5.] ] )
        
        pm0=LinearMappingX(s) 
        
        sm1=SimpleModel(dom, numComps=2,  coordinates=self.COORDINATES)

        reg=Regularization(dom, numLevelSets=2,
                               w1=[[1,1.],[1.,1.]], # consider gradient terms
                               wc=[[0,0],[0,0]],    # and cross-gradient term
                               coordinates=self.COORDINATES)
                               #location_of_set_m=reg_mask)
        cf=InversionCostFunction(reg, [ (pm0, (0,1))], [ (sm1,0)  ] )
        
        
        m=Data([1.,2], Solution(dom))
        args=cf.getArguments(m)
        
        df=cf.getValue(m, *args)
        df_real=sm1.getRealValue(np.dot(pm0.matrix, np.array([1,2])))
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

