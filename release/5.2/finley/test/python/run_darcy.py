# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import *
from esys.escript.models import DarcyFlow
from esys.finley import Rectangle, Brick

import os

VERBOSE=False

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

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
        k=1.e-8
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u_ref,p)

        
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p )

        
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")

    def testConstF_FixedBottom_largeK(self):
        k=1.e8
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")

    def testVarioF_FixedBottom_smallK(self):
        k=1.e-8
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        #df.getSolverOptionsPressure().setVerbosityOn()
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)

        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")

    def testVarioF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref  *mp
        u=u_ref *mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p) 

        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
 
    def testVarioF_FixedBottom_largeK(self):
        k=1.e8
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)

        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_smallK(self):
        k=1.e-8
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref # *mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                    location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
                     
        v,p=df.solve(u,p)
  
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref     *mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)

        
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_largeK(self):
        k=1.e8
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        #df.getSolverOptionsPressure().setVerbosityOn()
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")

    def testVarioF_FreeBottom_smallK(self):
        k=1.e-8
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)

        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")  

    def testVarioF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_largeK(self):
        k=1.e8
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom, verbose=VERBOSE, solver=self.SOLVER)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        v,p=df.solve(u,p)
        self.assertTrue(Lsup(v-u_ref)<self.TEST_TOL*Lsup(u_ref), "flux error too big.")
        self.assertTrue(Lsup(p-p_ref)<self.TEST_TOL*Lsup(p_ref), "pressure error too big.")

class Darcy2D(Darcy): #subclassing required
    TOL=1e-4
    WIDTH=1.
    SOLVER=DarcyFlow.POST
    def setUp(self):
        NE=40  # warning: smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Rectangle(NE,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom

class Test_Darcy2D_EVAL(Darcy2D):
    TEST_TOL=0.01
    SOLVER=DarcyFlow.EVAL

class Test_Darcy2D_POST(Darcy2D):
    TEST_TOL=1.e-3
    SOLVER=DarcyFlow.POST

class Test_Darcy2D_SMOOTH(Darcy2D):
    TEST_TOL=0.01
    SOLVER=DarcyFlow.SMOOTH

class Darcy3D(Darcy): #subclassing required
    WIDTH=1.
    SOLVER=DarcyFlow.POST
    def setUp(self):
        NE=29  # warning: smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Brick(NE,NE,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom

class Test_Darcy3D_EVAL(Darcy3D):
    TEST_TOL=0.01
    SOLVER=DarcyFlow.EVAL

class Test_Darcy3D_POST(Darcy3D):
    TEST_TOL=1.e-3
    SOLVER=DarcyFlow.POST

class Test_Darcy3D_SMOOTH(Darcy3D):
    TEST_TOL=0.01
    SOLVER=DarcyFlow.SMOOTH
    
    
if __name__ == '__main__':
    run_tests(__name__, classes=[
            Test_Darcy3D_EVAL, Test_Darcy3D_POST, Test_Darcy3D_SMOOTH,
            Test_Darcy2D_EVAL, Test_Darcy2D_POST, Test_Darcy2D_SMOOTH
        ],exit_on_failure=True)


