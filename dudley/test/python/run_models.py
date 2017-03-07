# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
      
VERBOSE = False

from esys.escript import *
from esys.escript.models import PowerLaw, FaultSystem, DarcyFlow
from esys.dudley import Rectangle, Brick

from math import pi
import numpy, os, sys, tempfile


class Darcy(unittest.TestCase): #subclassing required
    # this is a simple test for the darcy flux problem
    #
    # 
    #  p = 1/k * ( 1/2* (fb-f0)/xb* x **2 + f0 * x - ub*x ) +  p0
    # 
    #  with f = (fb-f0)/xb* x + f0 
    #
    #    u = f - k * p,x = ub
    #
    #  we prescribe pressure at x=x0=0 to p0
    # 
    #  if we prescribe the pressure on the bottom x=xb we set
    # 
    #  pb= 1/k * ( 1/2* (fb-f0)* xb + f0 * xb - ub*xb ) +  p0 = 1/k * ((fb+f0)/2 - ub ) * xb  + p0
    # 
    #  which leads to ub = (fb+f0)/2-k*(pb-p0)/xb
    #
    def rescaleDomain(self):
        x=self.dom.getX().copy()
        for i in range(self.dom.getDim()):
             x_inf=inf(x[i])
             x_sup=sup(x[i])
             if i == self.dom.getDim()-1:
                x[i]=-self.WIDTH*(x[i]-x_sup)/(x_inf-x_sup)
             else:
                x[i]=self.WIDTH*(x[i]-x_inf)/(x_sup-x_inf)
        self.dom.setX(x)

    def getScalarMask(self,include_bottom=True):
        x=self.dom.getX().copy()
        x_inf=inf(x[self.dom.getDim()-1])
        x_sup=sup(x[self.dom.getDim()-1])
        out=whereZero(x[self.dom.getDim()-1]-x_sup)
        if include_bottom: out+=whereZero(x[self.dom.getDim()-1]-x_inf)
        return wherePositive(out)

    def getVectorMask(self,include_bottom=True):
        x=self.dom.getX().copy()
        out=Vector(0.,Solution(self.dom))
        for i in range(self.dom.getDim()):
             x_inf=inf(x[i])
             x_sup=sup(x[i])
             if i != self.dom.getDim()-1: out[i]+=whereZero(x[i]-x_sup)
             if i != self.dom.getDim()-1 or include_bottom: out[i]+=whereZero(x[i]-x_inf)
        return wherePositive(out)

    def setSolutionFixedBottom(self, p0, pb, f0, fb, k):
         d=self.dom.getDim()
         x=self.dom.getX()[d-1]
         xb=inf(x)
         u=Vector(0.,Solution(self.dom))+kronecker(d)[d-1]*((f0+fb)/2.-k*(pb-p0)/xb)
         p=1./k*((fb-f0)/(xb*2.)* x**2 - (fb-f0)/2.*x)+(pb-p0)/xb*x +  p0
         f= ((fb-f0)/xb* x + f0)*kronecker(Function(self.dom))[d-1]
         return u,p,f
        
    def testConstF_FixedBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u_ref,p)
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
    def testConstF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p )
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")

    def testConstF_FixedBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p )
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        #df.getSolverOptionsPressure().setVerbosityOn()
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p )
        
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p )
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p )
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                    location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)

        
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)
        self.assertLess(Lsup(v-u_ref), self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref), self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)
        self.assertLess(Lsup(v-u_ref),self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref),self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)
        self.assertLess(Lsup(v-u_ref),self.TEST_TOL*Lsup(u_ref), "flux error too big.")  
        self.assertLess(Lsup(p-p_ref),self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)
        self.assertLess(Lsup(v-u_ref),self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref),self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.getSolverOptionsPressure().setTolerance(rtol=self.TOL)
        v,p=df.solve(u,p)
        self.assertLess(Lsup(v-u_ref),self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertLess(Lsup(p-p_ref),self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

class Test_Darcy2DOnDudley(Darcy):
    TOL=1e-8 if hasFeature('paso') else 1e-9
    TEST_TOL=2.e-3
    WIDTH=1.
    def setUp(self):
        NE=40  # wrning smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Rectangle(NE/2,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom

class Test_Darcy3DOnDudley(Darcy):
    TOL=1e-8 if hasFeature('paso') else 1e-9
    WIDTH=1.
    TEST_TOL=4.e-3
    def setUp(self):
        NE=25  # wrning smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Brick(NE,NE,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom

class Test_RheologiesOnDudley(unittest.TestCase):
     """
     this is the program used to generate the powerlaw tests:

     TAU_Y=100.
     N=10
     M=5

     def getE(tau):
         if tau<=TAU_Y:
           return 1./(0.5+20*sqrt(tau))
         else:
           raise ValueError,"out of range."
     tau=[ i* (TAU_Y/N) for i in range(N+1)] + [TAU_Y for i in range(M)]
     e=[ tau[i]/getE(tau[i]) for i in range(N+1)]
     e+= [ (i+1+j)* (max(e)/N) for j in range(M) ]

     print tau
     print e
     """
     TOL=1.e-8
     def test_PowerLaw_Init(self):
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)

         self.assertEqual(pl.getNumMaterials(),3,"num materials is wrong")
         self.assertEqual(pl.validMaterialId(0),True,"material id 0 not found")
         self.assertEqual(pl.validMaterialId(1),True,"material id 1 not found")
         self.assertEqual(pl.validMaterialId(2),True,"material id 2 not found")
         self.assertEqual(pl.validMaterialId(3),False,"material id 3 not found")

         self.assertEqual(pl.getFriction(),None,"initial friction wrong.")
         self.assertEqual(pl.getTauY(),None,"initial tau y wrong.")
         pl.setDruckerPragerLaw(tau_Y=10,friction=3)
         self.assertEqual(pl.getFriction(),3,"friction wrong.")
         self.assertEqual(pl.getTauY(),10,"tau y wrong.")

         self.assertEqual(pl.getElasticShearModulus(),None,"initial shear modulus is wrong")
         pl.setElasticShearModulus(1000)
         self.assertEqual(pl.getElasticShearModulus(),1000.,"shear modulus is wrong")

         e=pl.getEtaN()
         self.assertEqual(len(e),3,"initial length of etaN is wrong.")
         self.assertEqual(e,[None, None, None],"initial etaN are wrong.")
         self.assertEqual(pl.getEtaN(0),None,"initial material id 0 is wrong")
         self.assertEqual(pl.getEtaN(1),None,"initial material id 1 is wrong")
         self.assertEqual(pl.getEtaN(2),None,"initial material id 2 is wrong")
         self.assertRaises(ValueError, pl.getEtaN, 3)

         self.assertRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30,40],tau_t=[2], power=[5,6,7])
         self.assertRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30],tau_t=[2,4,5,5], power=[5,6,7])
         self.assertRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30],tau_t=[2,1,2], power=[5,6,7,7])

         pl.setPowerLaws(eta_N=[10,20,30],tau_t=[100,200,300], power=[1,2,3])
         self.assertEqual(pl.getPower(),[1,2,3],"powers are wrong.")
         self.assertEqual(pl.getTauT(),[100,200,300],"tau t are wrong.")
         self.assertEqual(pl.getEtaN(),[10,20,30],"etaN are wrong.")

         pl.setPowerLaw(id=0,eta_N=40,tau_t=400, power=4)
         self.assertEqual(pl.getPower(),[4,2,3],"powers are wrong.")
         self.assertEqual(pl.getTauT(),[400,200,300],"tau t are wrong.")
         self.assertEqual(pl.getEtaN(),[40,20,30],"etaN are wrong.")

         self.assertRaises(ValueError,pl.getPower,-1)
         self.assertEqual(pl.getPower(0),4,"power 0 is wrong.")
         self.assertEqual(pl.getPower(1),2,"power 1 is wrong.")
         self.assertEqual(pl.getPower(2),3,"power 2 is wrong.")
         self.assertRaises(ValueError,pl.getPower,3)

         self.assertRaises(ValueError,pl.getTauT,-1)
         self.assertEqual(pl.getTauT(0),400,"tau t 0 is wrong.")
         self.assertEqual(pl.getTauT(1),200,"tau t 1 is wrong.")
         self.assertEqual(pl.getTauT(2),300,"tau t 2 is wrong.")
         self.assertRaises(ValueError,pl.getTauT,3)

         self.assertRaises(ValueError,pl.getEtaN,-1)
         self.assertEqual(pl.getEtaN(0),40,"eta n 0 is wrong.")
         self.assertEqual(pl.getEtaN(1),20,"eta n 1 is wrong.")
         self.assertEqual(pl.getEtaN(2),30,"eta n 2 is wrong.")
         self.assertRaises(ValueError,pl.getEtaN,3)

     def checkResult(self,id,gamma_dot_, eta, tau_ref):
         self.assertTrue(eta>=0,"eta needs to be positive (test %s)"%id)
         error=abs(gamma_dot_*eta-tau_ref)
         self.assertTrue(error<=self.TOL*tau_ref,"eta is wrong: error = gamma_dot_*eta-tau_ref = %s * %s - %s = %s (test %s)"%(gamma_dot_,eta,tau_ref,error,id))
        
     def test_PowerLaw_Linear(self):
         taus= [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0]
         pl=PowerLaw(numMaterials=1,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaw(eta_N=2.)
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])
        
     def test_PowerLaw_QuadLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 405.0, 1610.0, 3615.0, 6420.0, 10025.0, 14430.0, 19635.0, 25640.0, 32445.0, 40050.0, 44055.0, 48060.0, 52065.0, 56070.0, 60075.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.4, 11.6, 18.6, 26.4, 35.0, 44.4, 54.6, 65.6, 77.4, 90.0, 99.0, 108.0, 117.0, 126.0, 135.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 8.90625, 41.25, 120.46875, 270.0, 513.28125, 873.75, 1374.84375, 2040.0, 2892.65625, 3956.25, 4351.875, 4747.5, 5143.125, 5538.75, 5934.375]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,1./16.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.0390625, 10.3125, 16.0546875, 22.5, 29.8828125, 38.4375, 48.3984375, 60.0, 73.4765625, 89.0625, 97.96875, 106.875, 115.78125, 124.6875, 133.59375]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,25./4.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 408.90625, 1641.25, 3720.46875, 6670.0, 10513.28125, 15273.75, 20974.84375, 27640.000000000004, 35292.65625, 43956.25, 48351.875, 52747.5, 57143.125, 61538.75, 65934.375]

         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 405.0390625, 1610.3125, 3616.0546875, 6422.5, 10029.8828125, 14438.4375, 19648.3984375, 25660.0, 32473.4765625, 40089.0625, 44097.96875, 48106.875, 52115.78125, 56124.6875, 60133.59375]

         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,25./4.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 9.30625, 42.85, 124.06875, 276.4, 523.28125, 888.15, 1394.44375, 2065.6, 2925.05625, 3996.25, 4395.875, 4795.5, 5195.125, 5594.75, 5994.375]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.4390625, 11.9125, 19.6546875, 28.9, 39.8828125, 52.8375, 67.9984375, 85.6, 105.8765625, 129.0625, 141.96875, 154.875, 167.78125, 180.6875, 193.59375]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.,25./4.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_withShear(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0]
         pl=PowerLaw(numMaterials=1,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaw(eta_N=2.)
         pl.setElasticShearModulus(3.)
         dt=1./3.
         pl.setEtaTolerance(self.TOL)
         self.assertRaises(ValueError, pl.getEtaEff,gamma_dot_s[0])
         for i in range(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i],dt=dt),taus[i])


class Test_FaultSystemOnDudley(unittest.TestCase):
   EPS=1.e-8
   NE=10
   def test_Fault_MaxValue(self):
      dom=Rectangle(2*self.NE,2*self.NE)
      x=dom.getX()
      f=FaultSystem(dim=2)
      f.addFault(V0=[0.5,0.], strikes=[3.*pi/4], ls=[0.70710678118654757], tag=1)
      f.addFault(V0=[1.,0.5], strikes=[pi, pi/2], ls=[0.5,0.5], tag=2)

      u=x[0]*(1.-x[0])*(1-x[1])
      t, loc=f.getMaxValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertEqual(m, 0.25, "wrong max value")
      self.assertEqual(t, 1, "wrong max tag")
      self.assertEqual(l, 0., "wrong max location")

      u=x[1]*(1.-x[1])*(1-x[0])*x[0]
      t, loc=f.getMaxValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertEqual(m, 0.0625, "wrong max value")
      self.assertEqual(t, 2, "wrong max tag")
      self.assertEqual(l, 0.5, "wrong max location")

      u=x[0]*(1.-x[0])*x[1]
      t, loc=f.getMaxValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertEqual(m, 0.25, "wrong max value")
      self.assertEqual(t, 2, "wrong max tag")
      self.assertEqual(l, 1.0, "wrong max location")

      u=x[1]*(1.-x[1])*x[0]
      t, loc=f.getMaxValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertEqual(m, 0.25, "wrong max value")
      self.assertEqual(t, 2, "wrong max tag")
      self.assertEqual(l, 0., "wrong max location")

      u=x[1]*(1.-x[1])*(1.-x[0])
      t, loc=f.getMaxValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertEqual(m, 0.25, "wrong max value")
      self.assertEqual(t, 1, "wrong max tag")
      self.assertLess(abs(l-0.70710678118654), self.EPS,  "wrong max location")

   def test_Fault_MinValue(self):
      dom=Rectangle(2*self.NE,2*self.NE)
      x=dom.getX()
      f=FaultSystem(dim=2)
      f.addFault(V0=[0.5,0.], strikes=[3.*pi/4], ls=[0.70710678118654757], tag=1)
      f.addFault(V0=[1.,0.5], strikes=[pi, pi/2], ls=[0.5,0.5], tag=2)

      u=-x[0]*(1.-x[0])*(1-x[1])
      t, loc=f.getMinValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertTrue(  m == -0.25, "wrong min value")
      self.assertTrue(  t == 1, "wrong min tag")
      self.assertTrue(  l == 0., "wrong min location")
      u=-x[1]*(1.-x[1])*(1-x[0])*x[0]
      t, loc=f.getMinValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertTrue(  m == -0.0625, "wrong min value")
      self.assertTrue(  t == 2, "wrong min tag")
      self.assertTrue(  l == 0.5, "wrong min location")
      u=-x[0]*(1.-x[0])*x[1]
      t, loc=f.getMinValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertTrue(  m == -0.25, "wrong min value")
      self.assertTrue(  t == 2, "wrong min tag")
      self.assertTrue(  l == 1.0, "wrong min location")
      u=-x[1]*(1.-x[1])*x[0]
      t, loc=f.getMinValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertTrue(  m == -0.25, "wrong min value")
      self.assertTrue(  t == 2, "wrong min tag")
      self.assertTrue(  l == 0., "wrong min location")
      u=-x[1]*(1.-x[1])*(1.-x[0])
      t, loc=f.getMinValue(u)
      p=f.getParametrization(x,t)[0]
      m, l=loc(u), loc(p)
      self.assertTrue(  m == -0.25, "wrong min value")
      self.assertTrue(  t == 1, "wrong min tag")
      self.assertTrue(  abs(l-0.70710678118654) <= self.EPS,  "wrong min location")

      
   def test_Fault2D(self):
      f=FaultSystem(dim=2)
      top1=[ [1.,0.], [1.,1.], [0.,1.] ]
      self.assertRaises(ValueError,f.addFault,V0=[1.,0],strikes=[pi/2, pi/2],ls=[1.,1.],tag=1,dips=top1)
      f.addFault(V0=[1.,0],strikes=[pi/2, pi],ls=[1.,1.],tag=1)
      self.assertTrue(f.getDim() == 2, "wrong dimension")
      self.assertTrue( [ 1 ] == f.getTags(), "tags wrong")
      self.assertTrue(  2. == f.getTotalLength(1), "length wrong")
      self.assertTrue(  0. == f.getMediumDepth(1), "depth wrong")
      self.assertTrue( (0., 2.) ==  f.getW0Range(1)," wrong W0 range")
      self.assertTrue( (0., 0.) ==  f.getW1Range(1)," wrong W1 range")
      self.assertTrue( [0., 1., 2.] ==  f.getW0Offsets(1)," wrong W0 offsets")
      segs=f.getTopPolyline(1)
      self.assertTrue( len(segs) == 3, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 2, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[1.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 2, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[1.,1.]) < self.EPS, "wrong vertex. 1 ")
      self.assertTrue( isinstance(segs[2], numpy.ndarray), "wrong class of vertex 2")
      self.assertTrue( segs[2].size == 2, "seg 2 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[2]-[0.,1.]) < self.EPS, "wrong vertex. 2 ")
      c=f.getCenterOnSurface()
      self.assertTrue( isinstance(c, numpy.ndarray), "center has wrong class")
      self.assertTrue( c.size == 2, "center size is wrong")
      self.assertTrue( numpy.linalg.norm(c-[2./3.,2./3.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.assertTrue( abs(o+45.) < self.EPS, "wrong orientation.")

      top2=[ [10.,0.], [0.,10.] ]
      f.addFault(V0=[10.,0],strikes=[3.*pi/4],ls=[14.142135623730951], tag=2, w0_offsets=[0,20], w1_max=20)
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags wrong")
      self.assertTrue(  abs(f.getTotalLength(2)-14.1421356237) < self.EPS * 14.1421356237, "wrong length")
      self.assertTrue(  0. == f.getMediumDepth(2), "depth wrong")
      self.assertTrue( (0., 20.) ==  f.getW0Range(2)," wrong W0 range")
      self.assertTrue( (0., 0.) ==  f.getW1Range(2)," wrong W1 range")
      self.assertTrue( [0., 20.] ==  f.getW0Offsets(2)," wrong W0 offsets")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( numpy.linalg.norm(segs[0]-[10.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0.,10.]) < self.EPS, "wrong vertex. 1 ")
      c=f.getCenterOnSurface()
      self.assertTrue( isinstance(c, numpy.ndarray), "center has wrong class")
      self.assertTrue( c.size == 2, "center size is wrong")
      self.assertTrue( numpy.linalg.norm(c-[12./5.,12./5.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.assertTrue( abs(o+45.) < self.EPS, "wrong orientation.")

      s,d=f.getSideAndDistance([0.,0.], tag=1)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([0.,2.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([1.,2.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([2.,1.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([2.,0.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([0.,-1.], tag=1)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-1.41421356237)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([-1.,0], tag=1)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-1.41421356237)<self.EPS, "wrong distance.")


      f.transform(rot=-pi/2., shift=[-1.,-1.])
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags after transformation wrong")
      self.assertTrue(  2. == f.getTotalLength(1), "length after transformation wrong")
      self.assertTrue(  0. == f.getMediumDepth(1), "depth after transformation wrong")
      self.assertTrue( (0., 2.) ==  f.getW0Range(1)," wrong W0 after transformation range")
      self.assertTrue( (0., 0.) ==  f.getW1Range(1)," wrong W1 rangeafter transformation ")
      self.assertTrue( [0., 1., 2.] ==  f.getW0Offsets(1)," wrong W0 offsetsafter transformation ")
      segs=f.getTopPolyline(1)
      self.assertTrue( len(segs) == 3, "wrong number of segmentsafter transformation ")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.assertTrue( segs[0].size == 2, "seg 0 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[-1.,0.]) < self.EPS, "wrong vertex. 0  after transformation")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex  after transformation1")
      self.assertTrue( segs[1].size == 2, "seg 1 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0.,0.]) < self.EPS, "wrong vertex.  after transformation1 ")
      self.assertTrue( isinstance(segs[2], numpy.ndarray), "wrong class of vertex  after transformation2")
      self.assertTrue( segs[2].size == 2, "seg 2 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[2]-[0., 1.]) < self.EPS, "wrong vertex after transformation. 2 ")
      self.assertTrue(  abs(f.getTotalLength(2)-14.1421356237) < self.EPS * 14.1421356237, "wrong length after transformation")
      self.assertTrue(  0. == f.getMediumDepth(2), "depth wrong after transformation")
      self.assertTrue( (0., 20.) ==  f.getW0Range(2)," wrong W0 range after transformation")
      self.assertTrue( (0., 0.) ==  f.getW1Range(2)," wrong W1 range after transformation")
      self.assertTrue( [0., 20.] ==  f.getW0Offsets(2)," wrong W0 offsets after transformation")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments after transformation")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.assertTrue( numpy.linalg.norm(segs[0]-[-1.,-9]) < self.EPS, "wrong vertex. 0  after transformation")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1 after transformation")
      self.assertTrue( numpy.linalg.norm(segs[1]-[9.,1.]) < self.EPS, "wrong vertex. 1  after transformation")

      c=f.getCenterOnSurface()
      self.assertTrue( isinstance(c, numpy.ndarray), "center has wrong class")
      self.assertTrue( c.size == 2, "center size is wrong")
      self.assertTrue( numpy.linalg.norm(c-[7./5.,-7./5.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.assertTrue( abs(o-45.) < self.EPS, "wrong orientation.")

      p=f.getParametrization([-1.,0.],1)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0])<self.EPS, "wrong value.")
      p=f.getParametrization([-0.5,0.],1)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-0.5)<self.EPS* 0.5, "wrong value.")
      p=f.getParametrization([0.,0.],1)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-1.)<self.EPS, "wrong value.")
      p=f.getParametrization([0.0000001,0.0000001],1, tol=1.e-8)
      self.assertTrue(p[1]==0., "wrong value.")
      p=f.getParametrization([0.0000001,0.0000001],1, tol=1.e-6)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-1.0000001)<self.EPS, "wrong value.")
      p=f.getParametrization([0.,0.5],1)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-1.5)<self.EPS, "wrong value.")
      p=f.getParametrization([0,1.],1)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-2.)<self.EPS, "wrong value.")
      p=f.getParametrization([1.,1.],1)
      self.assertTrue(p[1]==0., "wrong value.")
      p=f.getParametrization([0,1.11],1)
      self.assertTrue(p[1]==0., "wrong value.")
      p=f.getParametrization([-1,-9.],2)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0])<self.EPS, "wrong value.")
      p=f.getParametrization([9,1],2)
      self.assertTrue(p[1]==1., "wrong value.")
      self.assertTrue(abs(p[0]-20.)<self.EPS, "wrong value.")

   def test_Fault3D(self):
      f=FaultSystem(dim=3)
      self.assertTrue(f.getDim() == 3, "wrong dimension")

      top1=[ [0.,0.,0.], [1., 0., 0.] ]
      f.addFault(V0=[0.,0,0],strikes=[0.],ls=[1.,], dips=pi/2, depths=20.,tag=1)
      self.assertTrue( [ 1 ] == f.getTags(), "tags wrong")
      self.assertTrue(  1. == f.getTotalLength(1), "length wrong")
      self.assertTrue(  20. == f.getMediumDepth(1), "depth wrong")
      self.assertTrue( (0., 1.) ==  f.getW0Range(1)," wrong W0 range")
      self.assertTrue( (-20., 0.) ==  f.getW1Range(1)," wrong W1 range")
      self.assertTrue( [0., 1.] ==  f.getW0Offsets(1)," wrong W0 offsets")
      segs=f.getTopPolyline(1)
      self.assertTrue( len(segs) == 2, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[0.,0.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[1.,0.,0]) < self.EPS, "wrong vertex. 1 ")
      c=f.getCenterOnSurface()
      self.assertTrue( isinstance(c, numpy.ndarray), "center has wrong class")
      self.assertTrue( c.size == 3, "center size is wrong")
      self.assertTrue( numpy.linalg.norm(c-[0.5,0.,0.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.assertTrue( abs(o) < self.EPS, "wrong orientation.")
      d=f.getDips(1)
      self.assertTrue( len(d) == 1, "wrong number of dips")
      self.assertTrue(  abs(d[0]-1.5707963267948966) < self.EPS, "wrong dip 0")
      sn=f.getSegmentNormals(1)
      self.assertTrue( len(sn) == 1, "wrong number of normals")
      self.assertTrue( isinstance(sn[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[0]-[0, -1., 0.]) < self.EPS, "wrong bottom vertex 1 ")
      dv=f.getDepthVectors(1)
      self.assertTrue( len(dv) == 2, "wrong number of depth vectors.")
      self.assertTrue( isinstance(dv[0], numpy.ndarray), "wrong class of depth vector 0")
      self.assertTrue( numpy.linalg.norm(dv[0]-[0., 0., -20.]) < self.EPS, "wrong depth vector 0 ")
      self.assertTrue( isinstance(dv[1], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[1]-[0., 0., -20.]) < self.EPS, "wrong depth vector 1 ")
      b=f.getBottomPolyline(1)
      self.assertTrue( len(b) == 2, "wrong number of bottom vertices")
      self.assertTrue( isinstance(b[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(b[0]-[0., 0., -20.]) < self.EPS, "wrong bottom vertex 0 ")
      self.assertTrue( isinstance(b[1], numpy.ndarray), "wrong class of bottom vertex 1")
      self.assertTrue( numpy.linalg.norm(b[1]-[1., 0., -20.]) < self.EPS, "wrong bottom vertex 1 ")
      ds=f.getDepths(1)
      self.assertTrue( len(ds) == 2, "wrong number of depth")
      self.assertTrue( abs(ds[0]-20.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-20.) < self.EPS, "wrong depth at vertex 1 ")

      top2=[ [0.,0.,0.], [0., 10., 0.] ]
      f.addFault(V0=[0.,0,0],strikes=[pi/2],ls=[10.,], dips=pi/2, depths=20.,tag=2)
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags wrong")
      self.assertTrue(  10. == f.getTotalLength(2), "length wrong")
      self.assertTrue(  20. == f.getMediumDepth(2), "depth wrong")
      self.assertTrue( (0., 10.) ==  f.getW0Range(2)," wrong W0 range")
      self.assertTrue( (-20., 0.) ==  f.getW1Range(2)," wrong W1 range")
      self.assertTrue( [0., 10.] ==  f.getW0Offsets(2)," wrong W0 offsets")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[0.,0.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0., 10., 0]) < self.EPS, "wrong vertex. 1 ")
      d=f.getDips(2)
      self.assertTrue( len(d) == 1, "wrong number of dips")
      self.assertTrue(  abs(d[0]-1.5707963267948966) < self.EPS, "wrong dip 0")
      sn=f.getSegmentNormals(2)
      self.assertTrue( len(sn) == 1, "wrong number of normals")
      self.assertTrue( isinstance(sn[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[0]-[1, 0., 0.]) < self.EPS, "wrong bottom vertex 1 ")
      dv=f.getDepthVectors(2)
      self.assertTrue( len(dv) == 2, "wrong number of depth vectors.")
      self.assertTrue( isinstance(dv[0], numpy.ndarray), "wrong class of depth vector 0")
      self.assertTrue( numpy.linalg.norm(dv[0]-[0., 0., -20.]) < self.EPS, "wrong depth vector 0 ")
      self.assertTrue( isinstance(dv[1], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[1]-[0., 0., -20.]) < self.EPS, "wrong depth vector 1 ")
      b=f.getBottomPolyline(2)
      self.assertTrue( len(b) == 2, "wrong number of bottom vertices")
      self.assertTrue( isinstance(b[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(b[0]-[0., 0., -20.]) < self.EPS, "wrong bottom vertex 0 ")
      self.assertTrue( isinstance(b[1], numpy.ndarray), "wrong class of bottom vertex 1")
      self.assertTrue( numpy.linalg.norm(b[1]-[0., 10., -20.]) < self.EPS, "wrong bottom vertex 1 ")
      ds=f.getDepths(2)
      self.assertTrue( len(ds) == 2, "wrong number of depth")
      self.assertTrue( abs(ds[0]-20.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-20.) < self.EPS, "wrong depth at vertex 1 ")

      top2=[ [10.,0.,0.], [0., 10., 0.] ]
      f.addFault(V0=[10.,0,0],strikes=3*pi/4,ls=14.142135623730951, dips=pi/2, depths=30.,tag=2)
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags wrong")
      self.assertTrue(  abs(14.142135623730951 - f.getTotalLength(2)) <self.EPS, "length wrong")
      self.assertTrue(  30. == f.getMediumDepth(2), "depth wrong")
      self.assertTrue( (-30., 0.) ==  f.getW1Range(2)," wrong W1 range")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[10.,0.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0., 10., 0]) < self.EPS, "wrong vertex. 1 ")
      d=f.getDips(2)
      self.assertTrue( len(d) == 1, "wrong number of dips")
      self.assertTrue(  abs(d[0]-1.5707963267948966) < self.EPS, "wrong dip 0")
      sn=f.getSegmentNormals(2)
      self.assertTrue( len(sn) == 1, "wrong number of normals")
      self.assertTrue( isinstance(sn[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[0]-[0.70710678118654746, 0.70710678118654746, 0.]) < self.EPS, "wrong bottom vertex 1 ")
      dv=f.getDepthVectors(2)
      self.assertTrue( len(dv) == 2, "wrong number of depth vectors.")
      self.assertTrue( isinstance(dv[0], numpy.ndarray), "wrong class of depth vector 0")
      self.assertTrue( numpy.linalg.norm(dv[0]-[0., 0., -30.]) < self.EPS, "wrong depth vector 0 ")
      self.assertTrue( isinstance(dv[1], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[1]-[0., 0., -30.]) < self.EPS, "wrong depth vector 1 ")
      b=f.getBottomPolyline(2)
      self.assertTrue( len(b) == 2, "wrong number of bottom vertices")
      self.assertTrue( isinstance(b[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(b[0]-[10., 0., -30.]) < self.EPS, "wrong bottom vertex 0 ")
      self.assertTrue( isinstance(b[1], numpy.ndarray), "wrong class of bottom vertex 1")
      self.assertTrue( numpy.linalg.norm(b[1]-[0., 10., -30.]) < self.EPS, "wrong bottom vertex 1 ")
      ds=f.getDepths(2)
      self.assertTrue( len(ds) == 2, "wrong number of depth")
      self.assertTrue( abs(ds[0]-30.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-30.) < self.EPS, "wrong depth at vertex 1 ")

      top2=[ [10.,0.,0.], [0., 10., 0.] ]
      f.addFault(V0=[10.,0,0],strikes=3*pi/4,ls=14.142135623730951, dips=pi/4, depths=50.,tag=2)
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags wrong")
      self.assertTrue(  abs(14.142135623730951 - f.getTotalLength(2)) <self.EPS, "length wrong")
      self.assertTrue(  50. == f.getMediumDepth(2), "depth wrong")
      self.assertTrue( (-50., 0.) ==  f.getW1Range(2)," wrong W1 range")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[10.,0.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0., 10., 0]) < self.EPS, "wrong vertex. 1 ")
      d=f.getDips(2)
      self.assertTrue( len(d) == 1, "wrong number of dips")
      self.assertTrue(  abs(d[0]-0.78539816339744828) < self.EPS, "wrong dip 0")
      sn=f.getSegmentNormals(2)
      self.assertTrue( len(sn) == 1, "wrong number of normals")
      self.assertTrue( isinstance(sn[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[0]-[0.5,0.5,0.70710678118654746]) < self.EPS, "wrong bottom vertex 1 ")
      dv=f.getDepthVectors(2)
      self.assertTrue( len(dv) == 2, "wrong number of depth vectors.")
      self.assertTrue( isinstance(dv[0], numpy.ndarray), "wrong class of depth vector 0")
      self.assertTrue( numpy.linalg.norm(dv[0]-[25., 25., -35.355339059327378]) < self.EPS, "wrong depth vector 0 ")
      self.assertTrue( isinstance(dv[1], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[1]-[25.,25., -35.355339059327378]) < self.EPS, "wrong depth vector 1 ")
      b=f.getBottomPolyline(2)
      self.assertTrue( len(b) == 2, "wrong number of bottom vertices")
      self.assertTrue( isinstance(b[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(b[0]-[35., 25., -35.355339059327378]) < self.EPS, "wrong bottom vertex 0 ")
      self.assertTrue( isinstance(b[1], numpy.ndarray), "wrong class of bottom vertex 1")
      self.assertTrue( numpy.linalg.norm(b[1]-[25, 35., -35.355339059327378]) < self.EPS, "wrong bottom vertex 1 ")
      ds=f.getDepths(2)
      self.assertTrue( len(ds) == 2, "wrong number of depth")
      self.assertTrue( abs(ds[0]-50.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-50.) < self.EPS, "wrong depth at vertex 1 ")

      top1=[ [10.,0.,0], [10.,10.,0], [0.,10.,0] ]
      f.addFault(V0=[10.,0.,0.],strikes=[pi/2, pi],ls=[10.,10.],tag=1, dips=pi/4, depths=20.)
      self.assertTrue(  20. == f.getTotalLength(1), "length wrong")
      self.assertTrue(  20. == f.getMediumDepth(1), "depth wrong")
      segs=f.getTopPolyline(1)
      self.assertTrue( len(segs) == 3, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[10.,0.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[10.,10.,0.]) < self.EPS, "wrong vertex. 1 ")
      self.assertTrue( isinstance(segs[2], numpy.ndarray), "wrong class of vertex 2")
      self.assertTrue( segs[2].size == 3, "seg 2 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[2]-[0.,10.,0.]) < self.EPS, "wrong vertex. 2 ")
      d=f.getDips(1)
      self.assertTrue( len(d) == 2, "wrong number of dips")
      self.assertTrue(  abs(d[0]-0.78539816339744828) < self.EPS, "wrong dip 0")
      self.assertTrue(  abs(d[1]-0.78539816339744828) < self.EPS, "wrong dip 0")
      ds=f.getDepths(1)
      self.assertTrue( len(ds) == 3, "wrong number of depth")
      self.assertTrue( abs(ds[0]-20.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-20.) < self.EPS, "wrong depth at vertex 1 ")
      sn=f.getSegmentNormals(1)
      self.assertTrue( len(sn) == 2, "wrong number of normals")
      self.assertTrue( isinstance(sn[0], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[0]-[0.70710678118654746,0.,0.70710678118654746]) < self.EPS, "wrong bottom vertex 1 ")
      self.assertTrue( isinstance(sn[1], numpy.ndarray), "wrong class of bottom vertex 0")
      self.assertTrue( numpy.linalg.norm(sn[1]-[0.,0.70710678118654746,0.70710678118654746]) < self.EPS, "wrong bottom vertex 1 ")
      dv=f.getDepthVectors(1)
      self.assertTrue( len(dv) == 3, "wrong number of depth vectors.") 
      self.assertTrue( isinstance(dv[0], numpy.ndarray), "wrong class of depth vector 0")
      self.assertTrue( numpy.linalg.norm(dv[0]-[14.142135623730951, 0., -14.142135623730951]) < self.EPS, "wrong depth vector 0 ")
      self.assertTrue( isinstance(dv[1], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[1]-[11.547005383792515,11.547005383792515, -11.547005383792515]) < self.EPS, "wrong depth vector 2 ")
      self.assertTrue( isinstance(dv[2], numpy.ndarray), "wrong class of depth vector 1")
      self.assertTrue( numpy.linalg.norm(dv[2]-[0.,14.142135623730951, -14.142135623730951]) < self.EPS, "wrong depth vector 2 ")
      segs=f.getBottomPolyline(1)
      self.assertTrue( len(segs) == 3, "wrong number of segments")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[24.142135623730951,0.,-14.142135623730951]) < self.EPS, "wrong vertex. 0 ")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[21.547005383792515,21.547005383792515, -11.547005383792515]) < self.EPS, "wrong vertex. 1 ")
      self.assertTrue( isinstance(segs[2], numpy.ndarray), "wrong class of vertex 2")
      self.assertTrue( segs[2].size == 3, "seg 2 has wrong size.")
      self.assertTrue( numpy.linalg.norm(segs[2]-[0., 24.142135623730951, -14.142135623730951]) < self.EPS, "wrong vertex. 2 ")
      self.assertTrue( abs(0.-f.getW0Range(1)[0]) <=self.EPS," wrong W0 range (0)")
      self.assertTrue( abs(31.857329272664341-f.getW0Range(1)[1]) <=self.EPS," wrong W0 range (1)")
      self.assertTrue( abs(-20. - f.getW1Range(1)[0]) <=self.EPS," wrong W1 range (0)")
      self.assertTrue( abs(0. - f.getW1Range(1)[1]) <=self.EPS," wrong W1 range (1)")
      self.assertTrue( abs(0.0-f.getW0Offsets(1)[0])<=self.EPS," wrong W0 offsets (0)")
      self.assertTrue( abs(15.92866463633217-f.getW0Offsets(1)[1])<=self.EPS," wrong W0 offsets (1)")
      self.assertTrue( abs(31.857329272664341-f.getW0Offsets(1)[2])<=self.EPS," wrong W0 offsets(2)")
      #
      #    ============ fresh start ====================
      #
      f.addFault(V0=[1.,0,0.],strikes=[pi/2, pi],ls=[1.,1.], dips=pi/2,depths=20,tag=1)
      f.addFault(V0=[10.,0,0],strikes=[3.*pi/4],ls=[14.142135623730951], tag=2, w0_offsets=[0,20], w1_max=20, dips=pi/2,depths=20)
      c=f.getCenterOnSurface()
      self.assertTrue( isinstance(c, numpy.ndarray), "center has wrong class")
      self.assertTrue( c.size == 3, "center size is wrong")
      self.assertTrue( numpy.linalg.norm(c-[12./5.,12./5.,0.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.assertTrue( abs(o+45.) < self.EPS, "wrong orientation.")

      f.transform(rot=-pi/2., shift=[-1.,-1.,0.])
      self.assertTrue( [ 1, 2 ] == f.getTags(), "tags after transformation wrong")
      self.assertTrue(  2. == f.getTotalLength(1), "length after transformation wrong")
      self.assertTrue(  20. == f.getMediumDepth(1), "depth after transformation wrong")
      rw0=f.getW0Range(1)
      self.assertTrue( len(rw0) ==2, "wo range has wrong length")
      self.assertTrue( abs(rw0[0]) < self.EPS,"W0 0 wrong.")
      self.assertTrue( abs(rw0[1]-2.) < self.EPS,"W0 1 wrong.")
      self.assertTrue( (-20., 0.) ==  f.getW1Range(1)," wrong W1 rangeafter transformation ")
      dips=f.getDips(1)
      self.assertTrue(len(dips) == 2, "wrong number of dips.")
      self.assertTrue( abs(dips[0]-1.5707963267948966) <= self.EPS, "wrong dip")
      self.assertTrue( abs(dips[1]-1.5707963267948966) <= self.EPS, "wrong dip")
      ds=f.getDepths(1)
      self.assertTrue( len(ds) == 3, "wrong number of depth")
      self.assertTrue( abs(ds[0]-20.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-20.) < self.EPS, "wrong depth at vertex 1 ")
      self.assertTrue( abs(ds[2]-20.) < self.EPS, "wrong depth at vertex 1 ")
      segs=f.getTopPolyline(1)
      self.assertTrue( len(segs) == 3, "wrong number of segmentsafter transformation ")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.assertTrue( segs[0].size == 3, "seg 0 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[0]-[-1.,0.,0.]) < self.EPS, "wrong vertex. 0  after transformation")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex  after transformation1")
      self.assertTrue( segs[1].size == 3, "seg 1 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[1]-[0.,0.,0.]) < self.EPS, "wrong vertex.  after transformation1 ")
      self.assertTrue( isinstance(segs[2], numpy.ndarray), "wrong class of vertex  after transformation2")
      self.assertTrue( segs[2].size == 3, "seg 2 has wrong size after transformation.")
      self.assertTrue( numpy.linalg.norm(segs[2]-[0., 1.,0.]) < self.EPS, "wrong vertex after transformation. 2 ")
      self.assertTrue(  abs(f.getTotalLength(2)-14.1421356237) < self.EPS * 14.1421356237, "wrong length after transformation")
      self.assertTrue(  20. == f.getMediumDepth(2), "depth wrong after transformation")
      rw0=f.getW0Range(2)
      self.assertTrue( len(rw0) ==2, "wo range has wrong length")
      self.assertTrue( abs(rw0[0]) < self.EPS,"W0 0 wrong.")
      self.assertTrue( abs(rw0[1]-20.) < self.EPS,"W0 1 wrong.")
      self.assertTrue( (-20., 0.) ==  f.getW1Range(2)," wrong W1 range after transformation")
      self.assertTrue( [0., 20.] ==  f.getW0Offsets(2)," wrong W0 offsets after transformation")
      dips=f.getDips(2)
      self.assertTrue(len(dips) == 1, "wrong number of dips.")
      self.assertTrue( abs(dips[0]-1.5707963267948966) <= self.EPS, "wrong dip")
      ds=f.getDepths(2)
      self.assertTrue( len(ds) == 2, "wrong number of depth")
      self.assertTrue( abs(ds[0]-20.) < self.EPS, "wrong depth at vertex 0 ")
      self.assertTrue( abs(ds[1]-20.) < self.EPS, "wrong depth at vertex 1 ")
      segs=f.getTopPolyline(2)
      self.assertTrue( len(segs) == 2, "wrong number of segments after transformation")
      self.assertTrue( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.assertTrue( numpy.linalg.norm(segs[0]-[-1.,-9,0.]) < self.EPS, "wrong vertex. 0  after transformation")
      self.assertTrue( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1 after transformation")
      self.assertTrue( numpy.linalg.norm(segs[1]-[9.,1.,0.]) < self.EPS, "wrong vertex. 1  after transformation")
      #
      #    ============ fresh start ====================
      #
      f=FaultSystem(dim=3)

      top1=[ [0.,0.,0.], [1., 0., 0.] ]
      f.addFault(V0=[0.,0,0],strikes=[0.],ls=[1.,], dips=pi/2, depths=1.,tag=1)
      top1=[ [10.,0.,0], [10.,10.,0], [0.,10.,0] ]
      f.addFault(V0=[10.,0.,0.],strikes=[pi/2, pi],ls=[10.,10.],tag=2, dips=pi/4, depths=20.)

      p,m=f.getParametrization([0.3,0.,-0.5],1)
      self.assertTrue(length(p-[0.3,-0.5]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([0.5,0.,-0.5],1)
      self.assertTrue(length(p-[0.5,-0.5]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([0.25,0.,-0.5],1)
      self.assertTrue(length(p-[0.25,-0.5]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([0.5,0.,-0.25],1)
      self.assertTrue(length(p-[0.5,-0.25]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([0.001,0.,-0.001],1)
      self.assertTrue(length(p-[0.001, -0.001]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([0.001,0.,0.001],1)
      self.assertTrue(m==0., "wrong value.")

      p,m=f.getParametrization([0.999,0.,0.001],1)
      self.assertTrue(m==0., "wrong value.")

      p,m=f.getParametrization([1.001,0.,-0.001],1)
      self.assertTrue(m==0., "wrong value.")
      p,m=f.getParametrization([1.001,0.,-0.1],1)
      self.assertTrue(m==0., "wrong value.")
      p,m=f.getParametrization([1.001,0.,-0.000000001],1)
      self.assertTrue(m==0., "wrong value.")

      p,m=f.getParametrization([0.999,0.,-0.001],1)
      self.assertTrue(length(p-[0.999, -0.001]) <= self.EPS, "wrong value.")
      self.assertTrue(m==1., "wrong value.")

      p,m=f.getParametrization([ 16.29252873 , 6.46410161 ,-6.29252873], 2, tol=1.e-7)
      self.assertTrue(m==1., "wrong value.")
      self.assertTrue(length(p-[4.7785993908996511, -10]) <= self.EPS*10., "wrong value.")
      p,m=f.getParametrization([15.77350269, 12.77350269, -5.77350269], 2, tol=1.e-7)
      self.assertTrue(m==1., "wrong value.")
      self.assertTrue(length(p-[11.150065245432518, -10]) <= self.EPS*10., "wrong value.")

      p,m=f.getParametrization([  3., 17.0710678, -7.0710678], 2, tol=1.e-7)
      self.assertTrue(m==1., "wrong value.")
      self.assertTrue(length(p-[27.078729881764687, -10]) <= self.EPS*10., "wrong value.")
      p,m=f.getParametrization([9.30940108, 16.55204176, -6.55204176], 2, tol=1.e-7)
      self.assertTrue(m==1., "wrong value.")
      self.assertTrue(length(p-[20.707264027231822, -10]) <= self.EPS*10., "wrong value.")

      p,m=f.getParametrization([ 21.54700538,  21.54700538, -11.54700538], 2, tol=1.e-7)
      self.assertTrue(m==1., "wrong value.")
      self.assertTrue(length(p-[15.92866463633217, -20]) <= self.EPS*10., "wrong value.")

      p,m=f.getParametrization([ 0.,0.,0.], 2, tol=1.e-7)
      self.assertTrue(m==0., "wrong value.")

      p,m=f.getParametrization([ 11.,11.,0.], 2, tol=1.e-7)
      self.assertTrue(m==0., "wrong value.")


      s,d=f.getSideAndDistance([0.,-1.,0.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([1.,-1.,0.], tag=1)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([0.,1.,0.], tag=1)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([1.,1.,0.], tag=1)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-1.)<self.EPS, "wrong distance.")

    
      s,d=f.getSideAndDistance([0.,0.,0.], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-10.)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([5.,5.,0.], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-5.)<self.EPS, "wrong distance.")

      s,d=f.getSideAndDistance([10.,10.,-1.], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-0.70710678118654757)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([10.,10.,-2.], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-2.*0.70710678118654757)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([10.,10.,-3.], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-3.*0.70710678118654757)<self.EPS, "wrong distance.")

      s,d=f.getSideAndDistance([5.,12.,0], tag=2)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-2*0.70710678118654757)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([5.,12.,-1], tag=2)
      self.assertTrue( s>0, "wrong side.")
      self.assertTrue( abs(d-0.70710678118654757)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([5.,12.,-2], tag=2)
      # s not checked as it is undefined.
      self.assertTrue( abs(d)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([5.,12.,-3], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-0.70710678118654757)<self.EPS, "wrong distance.")
      s,d=f.getSideAndDistance([5.,12.,-4], tag=2)
      self.assertTrue( s<0, "wrong side.")
      self.assertTrue( abs(d-2.*0.70710678118654757)<self.EPS, "wrong distance.")

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

