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

import logging
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import numpy as np
import os
import sys
import cmath
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.escript import *
from esys.weipa import saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE
from esys.escript.pdetools import Locator

try:
    from esys.ripley import Rectangle as ripRectangle, Brick as ripBrick
    HAVE_RIPLEY = True
except ImportError as e:
    HAVE_RIPLEY = False

try:
    from esys.finley import Rectangle as finRectangle
    HAVE_FINLEY = True
except ImportError as e:
    HAVE_FINLEY = False

mpisize = getMPISizeWorld()
# this is mainly to avoid warning messages
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

HAVE_DIRECT = hasFeature("PASO_DIRECT") or hasFeature('trilinos')

@unittest.skipUnless(HAVE_RIPLEY, "Ripley module not available")
@unittest.skipUnless(HAVE_DIRECT, "more than 1 MPI rank or missing direct solver")
class TestAcousticInversion(unittest.TestCase):
    def test_API(self):
        NE=20
        for x in [int(sqrt(mpisize)),2,3,5,7,1]:
            NX=x
            NY=mpisize//x
            if NX*NY == mpisize:
                break
        domain=ripRectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY, diracPoints=[(0.5,1.)], diracTags=['sss'])

        omega=2.

        data=Data([1,2], FunctionOnBoundary(domain))
        F=Data([2,3], Function(domain))
        w=1.
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, w, data, 1.) # F is a scalar
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, w, [1,2], F) # data is not Data
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, w, Data([1,2], Function(domain)), F) # data is not on boundary
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, w, Scalar(1, Function(domain)), F) # data is not of shape (2,)
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, [1,2], data, F) # w is not a scalar
        self.assertRaises(ValueError, AcousticWaveForm, domain, omega, Scalar(1, Function(domain)), data, F) # w is not a scalar

        # now we do a real one
        acw=AcousticWaveForm(domain, omega, w, data, F)
        self.assertEqual(acw.getDomain(),  domain)
        pde=acw.setUpPDE()
        self.assertIsInstance(pde, LinearPDE)
        self.assertEqual(pde.getNumEquations(), 2)
        self.assertEqual(pde.getNumSolutions(), 2)
        self.assertEqual(pde.getDomain(),  domain)

    def test_numeric2DscaleF(self):
        NE=40
        for x in [int(sqrt(mpisize)),2,3,5,7,1]:
            NX=x
            NY=mpisize//x
            if NX*NY == mpisize:
                break
        domain=ripRectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY, diracPoints=[(0.5,1.)], diracTags=['sss'])
        #domain=ripRectangle(100,100, diracPoints=[(0.5,1.)], diracTags=['sss'])
        omega=2.

        # test solution is u = a * z where a is complex
        a=complex(3.45, 0.56)
        sigma=complex(1e-3, 0.056)

        data=Data([a.real, a.imag], FunctionOnBoundary(domain))
        mydata=data.copy()

        z=FunctionOnBoundary(domain).getX()[1]
        w=whereZero(z-1.)
        # source:
        F=Data( [1,0],Function(domain))

        acw=AcousticWaveForm(domain, omega, w, data, F, coordinates=None, fixAtBottom=False, tol=1e-8, saveMemory=True, scaleF=True)
        # check rescaled data
        surv=acw.getSurvey()
        self.assertAlmostEqual( integrate(length(surv[0])**2 * surv[1]), 1.)

        mydata_scale=sqrt( integrate(w*length(mydata)**2) )
        self.assertAlmostEqual( acw.getSourceScaling(z*[1, 0.]) , a/mydata_scale )
        self.assertAlmostEqual( acw.getSourceScaling(mydata) , 1./mydata_scale )

        # this should be zero:
        sigma_comps=[sigma.real, sigma.imag]
        args=acw.getArguments(sigma_comps)
        d=acw.getDefect(sigma_comps, *args)
        self.assertTrue(isinstance(d, float))
        self.assertLess(abs(d), 1e-10)

        dg=acw.getGradient(sigma_comps, *args)
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==(2,))
        self.assertTrue(dg.getFunctionSpace()==Solution(domain))
        self.assertLess(Lsup(dg), 1e-10)

        # this shuld be zero' too
        sigma_comps=[2*sigma.real, sigma.imag/2.]
        args=acw.getArguments(sigma_comps)
        d=acw.getDefect(sigma_comps, *args)
        self.assertTrue(isinstance(d, float))
        self.assertLess(abs(d), 1e-10)

        dg=acw.getGradient(sigma_comps, *args)
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==(2,))
        self.assertTrue(dg.getFunctionSpace()==Solution(domain))
        self.assertLess(Lsup(dg), 1e-10)

        # this shouldn't be zero:
        sigma0=[2*sigma.real, 10*a.imag]*(27*Function(domain).getX()[0]-Function(domain).getX()[1])
        args=acw.getArguments(sigma0)
        d0=acw.getDefect(sigma0, *args)
        self.assertTrue(isinstance(d0, float))
        self.assertTrue(d0 >= 0)
        self.assertTrue(d0 > 1e-10)

        dg0=acw.getGradient(sigma0, *args)
        self.assertTrue(isinstance(dg0, Data))
        self.assertTrue(dg0.getShape()==(2,))
        self.assertTrue(dg0.getFunctionSpace()==Solution(domain))
        self.assertTrue(Lsup(dg0) > 1e-10)

        # test the gradient numerrically:
        h=0.002
        X=Function(domain).getX()
        # .. increment:
        p=h*exp(-(length(X-[0.6,0.6])/10)**2)*Lsup(length(sigma0))

        sigma1=sigma0+p*[1,0]
        args=acw.getArguments(sigma1)
        d1=acw.getDefect(sigma1, *args)
        self.assertLess( abs( d1-d0-integrate(dg0[0]*p) ), 1e-2*abs(d1-d0) )

        sigma2=sigma0+p*[0,1]
        args=acw.getArguments(sigma2)
        d2=acw.getDefect(sigma2, *args)
        self.assertLess( abs(d2-d0-integrate(dg0[1]*p)), 1e-2*abs(d2-d0) )

    def test_numeric2DnoscaleF(self):
        domain=ripRectangle(n0=10, n1=20, l0=1., l1=1., diracPoints=[(0.5,1.)], diracTags=['sss'])
        omega=1.5

        # test solution is u = a * z where a is complex
        a=complex(3.45, 0.56)
        sigma=complex(1e-3, 0.056)

        data=Data([a.real, a.imag], FunctionOnBoundary(domain))
        z=FunctionOnBoundary(domain).getX()[1]
        w=whereZero(z-1.)
        # F = - a*omega* sigma
        F=Data( [-(a*omega**2*sigma).real, -(a*omega**2*sigma).imag ],Function(domain))

        acw=AcousticWaveForm(domain, omega, w, data, F, coordinates=None, fixAtBottom=False, tol=1e-8, saveMemory=True, scaleF=False)
        # this should be zero:
        sigma_comps=[sigma.real, sigma.imag]
        args=acw.getArguments(sigma_comps)
        d=acw.getDefect(sigma_comps, *args)
        self.assertTrue(isinstance(d, float))
        self.assertLess(Lsup(d), 1e-10)
        #self.assertTrue(d >= 0)

        dg=acw.getGradient(sigma_comps, *args)

        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==(2,))
        self.assertTrue(dg.getFunctionSpace()==Solution(domain))
        self.assertLess(Lsup(dg), 1e-8)
        # this shouldn't be zero:
        sigma0=Data([2*sigma.real, sigma.imag/2], Function(domain) )
        args=acw.getArguments(sigma0)
        d0=acw.getDefect(sigma0, *args)
        self.assertTrue(isinstance(d0, float))
        self.assertGreaterEqual(d0, 0)
        self.assertGreater(d0, 1e-10)

        dg0=acw.getGradient(sigma0, *args)
        self.assertTrue(isinstance(dg0, Data))
        self.assertTrue(dg0.getShape()==(2,))
        self.assertTrue(dg0.getFunctionSpace()==Solution(domain))
        self.assertTrue(Lsup(dg0) > 1e-10)
        # test the gradient numerically:
        h=0.001
        X=Function(domain).getX()
        p=h*sin(length(X)*np.pi)*Lsup(length(sigma0))

        sigma1=sigma0+p*[1,0]
        args=acw.getArguments(sigma1)
        d1=acw.getDefect(sigma1, *args)
        self.assertLess( abs( d1-d0-integrate(dg0[0]*p) ), 1e-2*abs(d1-d0) )

        sigma2=sigma0+p*[0,1]
        args=acw.getArguments(sigma2)
        d2=acw.getDefect(sigma2, *args)
        self.assertLess( abs(d2-d0-integrate(dg0[1]*p)), 1e-2*abs(d2-d0) )


@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class TestSubsidence(unittest.TestCase):
    def test_PDE(self):
        lam=2.
        mu=1.

        domain=ripBrick(20,20,max(19,2*mpisize-1), d2=mpisize)

        xb=FunctionOnBoundary(domain).getX()
        m=whereZero(xb[2]-1)
        w=m*[0,0,1]
        d=m*2.5
        acw=Subsidence(domain, w,d, lam, mu )

        P0=10.
        args0=acw.getArguments(P0)
        u=args0[0]
        self.assertLess(Lsup(u[0]), 1.e-8)
        self.assertLess(Lsup(u[1]), 1.e-8)
        self.assertLess(Lsup(u[2]-2.5*domain.getX()[2]), 1.e-8)

        dd=acw.getDefect(P0, *args0)

        self.assertTrue( dd >= 0.)
        self.assertTrue( dd <= 1e-7 * 2.5 )

    def test_Differential(self):
        lam=2.
        mu=1.
        INC=0.01
        domain=ripBrick(20,20,min(99,20*mpisize-1) , d2=mpisize)

        xb=FunctionOnBoundary(domain).getX()
        m=whereZero(xb[2]-1)
        w=m*[0,0,1]
        d=m*2.5
        acw=Subsidence(domain, w,d, lam, mu )

        x=Function(domain).getX()
        P0=x[0]*x[1]
        args0=acw.getArguments(P0)
        d0=acw.getDefect(P0, *args0)
        grad_d=acw.getGradient(P0, *args0)

        dP=exp(-(length(x-[0.5,0.5,0.5])/0.06)**2)
        P1=P0+INC*dP
        args1=acw.getArguments(P1)
        d1=acw.getDefect(P1, *args1)
        ref=abs((d1-d0)/INC)
        self.assertLess(abs((d1-d0)/INC-integrate(grad_d* dP)), ref * 1.e-5)

        dP=exp(-(length(x-[0.3,0.3,0.5])/0.06)**2)
        P2=P0-INC*dP
        args2=acw.getArguments(P2)
        d2=acw.getDefect(P2, *args2)
        ref=abs((d2-d0)/INC)
        self.assertLess(abs((d2-d0)/INC+integrate(grad_d* dP)), ref * 1.e-5)

@unittest.skipUnless(HAVE_FINLEY, "Finley module not available")
class TestDCResistivity(unittest.TestCase):

    def test_PDE2D(self):
        dx_tests=0.1
        sigma0=1.
        electrodes=[(0.5-2*dx_tests,1.), (0.5-dx_tests,1.), (0.5+dx_tests,1.), (0.5+2*dx_tests,1.)]
        domain=finRectangle(20,20, d1=mpisize,  diracPoints=electrodes, diracTags=["sl0", "sl1", "sr0", "sr1"] )
        loc=Locator(domain,electrodes[2:])

        # this creates some reference Data:
        x=domain.getX()
        q=whereZero(x[0]-inf(x[0]))+whereZero(x[0]-sup(x[0]))+whereZero(x[1]-inf(x[1]))
        ppde=LinearPDE(domain, numEquations=1)
        s=Scalar(0.,DiracDeltaFunctions(domain))
        s.setTaggedValue("sl0" ,1.)
        s.setTaggedValue("sl1",-1.)
        ppde.setValue(A=kronecker(2)*sigma0, q=q, y_dirac=s)
        pp=ppde.getSolution()
        uu=loc(pp)

        # arguments for DcRes
        current = 10.
        sourceInfo = [ "sl0",  "sl1" ]
        sampleTags = [ ("sr0", "sr1") ]

        sigmaPrimary=7.
        phiPrimary=pp*current*sigma0/sigmaPrimary

        uuscale=1-current*sigma0/sigmaPrimary
        delphi_in = [ (uu[1]-uu[0]) * uuscale]

        acw=DcRes(domain, loc, delphi_in, sampleTags,  phiPrimary, sigmaPrimary)

        self.assertLess(Lsup(phiPrimary-acw.getPrimaryPotential()), 1.e-10 * Lsup(acw.getPrimaryPotential()))

        SIGMA=10. # matches current
        args0=acw.getArguments(SIGMA)
        p=args0[0]
        u=args0[1]

        # true secondary potential
        pps=pp-phiPrimary
        self.assertLess(Lsup(p-pps), 1.e-6*Lsup(pps))

        # test return values at electrodes:
        self.assertLess(abs(u[0]-uu[0]*uuscale), 1.e-6 * abs(uu[0]*uuscale))
        self.assertLess(abs(u[1]-uu[1]*uuscale), 1.e-6 * abs(uu[1]*uuscale))

        # this sould be zero
        dd=acw.getDefect(SIGMA, *args0)
        self.assertTrue( dd >= 0.)
        self.assertTrue( dd <= 1e-7 )

    def test_Differential2D(self):
        INC=0.001
        sigma0=1.
        dx_tests=0.1
        electrodes=[(0.5-2*dx_tests,1.), (0.5-dx_tests,1.), (0.5+dx_tests,1.), (0.5+2*dx_tests,1.)]
        domain=finRectangle(20,20, d1=mpisize,  diracPoints=electrodes, diracTags=["sl0", "sl1", "sr0", "sr1"] )
        loc=Locator(domain,electrodes[2:])

        # arguments for DcRes
        #current = 10.
        sampleTags = [ ("sr0", "sr1") ]

        delphi_in = [ 0.05 ]

        sigmaPrimary=1
        x=domain.getX()
        phiPrimary=(x[0]-inf(x[0]))*(x[1]-inf(x[1]))*(x[0]-sup(x[0]))

        acw=DcRes(domain, loc, delphi_in, sampleTags,  phiPrimary, sigmaPrimary)

        #=====================================================================
        x=Function(domain).getX()
        SIGMA0=x[0]*x[1]+1
        args0=acw.getArguments(SIGMA0)
        d0=acw.getDefect(SIGMA0, *args0)
        grad_d=acw.getGradient(SIGMA0, *args0)

        dS=exp(-(length(x-[0.5,0.5])/0.2)**2)
        SIGMA1=SIGMA0+INC*dS
        args1=acw.getArguments(SIGMA1)
        d1=acw.getDefect(SIGMA1, *args1)
        ref=abs((d1-d0)/INC)
        self.assertLess(abs((d1-d0)/INC-integrate(grad_d* dS)), ref * 1.e-3)

        dS=-exp(-(length(x-[0.5,0.5])/0.2)**2)
        SIGMA2=SIGMA0+INC*dS
        args2=acw.getArguments(SIGMA2)
        d2=acw.getDefect(SIGMA2, *args2)
        ref=abs((d2-d0)/INC)
        self.assertLess(abs((d2-d0)/INC-integrate(grad_d* dS)), ref * 1.e-3)

        dS=-1
        SIGMA3=SIGMA0+INC*dS
        args3=acw.getArguments(SIGMA3)
        d3=acw.getDefect(SIGMA3, *args3)
        ref=abs((d3-d0)/INC)
        self.assertLess(abs((d3-d0)/INC-integrate(grad_d* dS)), ref * 1.e-3)

        dS=1
        SIGMA4=SIGMA0+INC*dS
        args4=acw.getArguments(SIGMA4)
        d4=acw.getDefect(SIGMA4, *args4)
        ref=abs((d4-d0)/INC)
        self.assertLess(abs((d4-d0)/INC-integrate(grad_d* dS)), ref * 1.e-3)

class TestIsostaticPressure(unittest.TestCase):
    @unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
    def test_all(self):
        domain=ripBrick(50,50,20*mpisize-1, d2=mpisize)

        ps=IsostaticPressure(domain, level0=1., coordinates=None)
        g=Vector(0., Function(domain))
        rho=Scalar(100, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*981.
        self.assertLess(Lsup(p0-p_ref), 1e-6 * Lsup(p_ref))

        g=Vector([0,0,-10], Function(domain))
        rho=Scalar(0, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*26700
        self.assertLess(Lsup(p0-p_ref), 1e-6 * Lsup(p_ref))

        g=Vector([0,0,-10], Function(domain))
        rho=Scalar(100, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*(981.+26700+1000)
        self.assertLess(Lsup(p0-p_ref), 1e-6 * Lsup(p_ref))

@unittest.skipUnless(HAVE_RIPLEY, "Ripley module not available")
@unittest.skipUnless(HAVE_DIRECT, "more than 1 MPI rank or missing direct solver")
class TestMT2DModelTEMode(unittest.TestCase):
    def test_API(self):
        domain=ripRectangle(25, 25, d1=mpisize)
        omega=2.
        SIGMA=15.

        x=[ [0.2,0.5], [0.3,0.5] ]
        Z_XY=[ complex(1.2,1.5), complex(1.3,2.5) ]

        eta=1.
        w0=1.
        # now we do a real one
        model=MT2DModelTEMode(domain, omega, x, Z_XY, eta, sigma0=SIGMA, w0=w0)
        self.assertEqual(model.getDomain(),  domain)
        pde=model.setUpPDE()
        self.assertIsInstance(pde, LinearPDE)
        self.assertEqual(pde.getNumEquations(), 2)
        self.assertEqual(pde.getNumSolutions(), 2)
        self.assertEqual(pde.getDomain(),  domain)

        # other things that should work
        model=MT2DModelTEMode(domain, omega, x, Z_XY, eta=None, w0=[2.,3.] , sigma0=SIGMA, airLayerLevel=1.)

        # these shouldn't work
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, x, [3.], eta=[1.,1.], w0=[2.,3.],  sigma0=SIGMA)
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, x, Z_XY, eta=[1.], w0=[2.,3.] , sigma0=SIGMA)
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, [(6.7,5)], Z_XY, eta=[1.,1.], w0=[2.,3.] , sigma0=SIGMA)

    def test_PDE(self):
        omega=1.
        mu0=0.123
        SIGMA=15.
        k=cmath.sqrt(1j*omega*mu0*SIGMA)  # Ex=exp(k*z)
        NE=101
        domain=ripRectangle(max(NE,30*mpisize-1),max(NE,30*mpisize-1), d1=mpisize)

        Z0=0.5
        H=1./NE
        X1=(int(0.3/H)+0.5)*H
        X2=(int(0.6/H)+0.5)*H

        
        IMP=-(1j*omega*mu0)/k*cmath.exp(k*Z0)/cmath.exp(k*Z0)
        Z_XY=[ IMP, IMP ]

        x=[ [X1,Z0], [X2,Z0] ]
        eta=None

        z=domain.getX()[1]
        Ex0_ex=cos(k.imag*(z-1))*exp(k.real*(z-1))
        Ex1_ex=sin(k.imag*(z-1))*exp(k.real*(z-1))
        Ex0_ex_z=-sin(k.imag*(z-1))*k.imag*exp(k.real*(z-1))+cos(k.imag*(z-1))*exp(k.real*(z-1))*k.real
        Ex1_ex_z=cos(k.imag*(z-1))*k.imag*exp(k.real*(z-1))+sin(k.imag*(z-1))*exp(k.real*(z-1))*k.real

        model=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, tol=1e-9,  directSolver=True, sigma0=SIGMA)

        args=model.getArguments(SIGMA)
        Ex=args[0]
        Exz=args[1]
        self.assertLess(Lsup(Ex[0]-Ex0_ex), 1e-4 * Lsup(Ex0_ex))
        self.assertLess(Lsup(Ex[1]-Ex1_ex), 1e-4 * Lsup(Ex1_ex))
        self.assertLess(Lsup(Exz[0]-Ex0_ex_z), 1e-2 * Lsup(Ex0_ex_z))
        self.assertLess(Lsup(Exz[1]-Ex1_ex_z), 1e-2 * Lsup(Ex1_ex_z))

        argsr=model.getArguments(0.)
        ref=model.getDefect(0., *argsr)

        # this should be almost zero:
        args=model.getArguments(SIGMA)
        d=model.getDefect(SIGMA, *args)
        self.assertTrue( d > 0.)
        self.assertTrue( ref > 0.)
        self.assertLess( d, 3e-3 * ref ) # d should be zero (some sort of)

        z=ReducedFunction(domain).getX()[1]
        Ex0_ex=cos(k.imag*(z-1))*exp(k.real*(z-1))
        Ex1_ex=sin(k.imag*(z-1))*exp(k.real*(z-1))
        Ex0_ex_z=-sin(k.imag*(z-1))*k.imag*exp(k.real*(z-1))+cos(k.imag*(z-1))*exp(k.real*(z-1))*k.real
        Ex1_ex_z=cos(k.imag*(z-1))*k.imag*exp(k.real*(z-1))+sin(k.imag*(z-1))*exp(k.real*(z-1))*k.real
        # and this should be zero
        d0=model.getDefect(SIGMA, Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], Ex0_ex_z*[1.,0]+ Ex1_ex_z*[0,1.])
        self.assertTrue( d0 <= 1e-8 * ref ) # d should be zero (some sort of)

        # and this too
        dg=model.getGradient(SIGMA, Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], Ex0_ex_z*[1.,0]+ Ex1_ex_z*[0,1.])
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==())
        self.assertLess(Lsup(dg), 1e-10)

    def test_Differential(self):
        INC=0.001
        omega=5.
        mu0=0.123
        SIGMA=15.
        k=cmath.sqrt(1j*omega*mu0*SIGMA)  # Ex=exp(k*z)

        NE=101
        domain=ripRectangle(max(NE,50*mpisize-1), max(NE,50*mpisize-1), d1=mpisize)

        Z0=0.5
        IMP=-(1j*omega*mu0)/k*(cmath.exp(k*Z0)-cmath.exp(-k*Z0))/(cmath.exp(k*Z0)+cmath.exp(-k*Z0))
        Z_XY=[ IMP, IMP ]
        H=1./NE
        X1=(int(0.3/H)+0.5)*H
        X2=(int(0.6/H)+0.5)*H
        x=[ [X1,Z0], [X2,Z0] ]
        eta=None

        z=domain.getX()[1]
        Ex0_ex=cos(k.imag*z)*(exp(k.real*z)-exp(-k.real*z))
        Ex0_ex_z=-sin(k.imag*z)*k.imag*(exp(k.real*z)-exp(-k.real*z))+cos(k.imag*z)*(exp(k.real*z)+exp(-k.real*z))*k.real
        Ex1_ex=sin(k.imag*z)*(exp(k.real*z)+exp(-k.real*z))
        Ex1_ex_z=cos(k.imag*z)*k.imag*(exp(k.real*z)+exp(-k.real*z))+sin(k.imag*z)*(exp(k.real*z)-exp(-k.real*z))*k.real

        model=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, sigma0=SIGMA, tol=1e-9,  directSolver=True,  airLayerLevel=1.)

        # this is the base line:
        xx=domain.getX()[0]
        SIGMA0=3.*(xx+0.3)
        args0=model.getArguments(SIGMA0)
        d0=model.getDefect(SIGMA0, *args0)
        dg0=model.getGradient(SIGMA0, *args0)
        self.assertTrue(isinstance(dg0, Data))
        self.assertTrue(dg0.getShape()==())

        X=Function(domain).getX()

        # test 1
        p=INC
        SIGMA1=SIGMA0+p
        args1=model.getArguments(SIGMA1)
        d1=model.getDefect(SIGMA1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )

        # test 2
        p=exp(-length(X-(0.2,0.2))**2/10)*INC
        SIGMA1=SIGMA0+p
        args1=model.getArguments(SIGMA1)
        d1=model.getDefect(SIGMA1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )

        # test 3
        p=sin(length(X)*3*3.14)*INC
        SIGMA1=SIGMA0+p
        args1=model.getArguments(SIGMA1)
        d1=model.getDefect(SIGMA1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )


@unittest.skipUnless(HAVE_RIPLEY, "Ripley module not available")
@unittest.skipUnless(HAVE_DIRECT, "more than 1 MPI rank or missing direct solver")
class TestMT2DModelTMMode(unittest.TestCase):
    def test_API(self):
        domain=ripRectangle(25, 25, d0=mpisize)
        omega=2.
        SIGMA=15.
        x=[ [0.2,0.5], [0.3,0.5] ]
        Z_XY=[ complex(1.2,1.5), complex(1.3,2.5) ]
        eta=1.
        w0=1.
        Hx0=1.
        # now we do a real one
        model=MT2DModelTMMode(domain, omega, x, Z_XY, eta, w0=w0, sigma0=SIGMA)
        self.assertEqual(model.getDomain(),  domain)
        pde=model.setUpPDE()
        self.assertIsInstance(pde, LinearPDE)
        self.assertEqual(pde.getNumEquations(), 2)
        self.assertEqual(pde.getNumSolutions(), 2)
        self.assertEqual(pde.getDomain(),  domain)

        # other things that should work
        model=MT2DModelTMMode(domain, omega, x, Z_XY, eta=None, w0=[2.,3.], sigma0=SIGMA,  airLayerLevel=1.)

        # these shouldn't work
        self.assertRaises(ValueError, MT2DModelTMMode, domain, omega, x, [3.], eta=[1.,1.], w0=[2.,3.], sigma0=SIGMA)
        self.assertRaises(ValueError, MT2DModelTMMode, domain, omega, x, Z_XY, eta=[1.], w0=[2.,3.], sigma0=SIGMA)
        self.assertRaises(ValueError, MT2DModelTMMode, domain, omega, [(6.7,5)], Z_XY, eta=[1.,1.], w0=[2.,3.], sigma0=SIGMA)

    def test_PDE(self):
        omega=10.
        mu0=0.123
        RHO=0.15
        SIGMA=1/RHO
        k=cmath.sqrt(1j*omega*mu0/RHO)  # Hx=exp(k*z)
        NE=151
        L=1
        domain=ripRectangle(NE,NE, d0=mpisize)

        Z0=0.5
        H=1./NE
        X1=(int(0.3/H)+0.5)*H
        X2=(int(0.6/H)+0.5)*H

        IMP=RHO*k*cmath.exp(k*(Z0-L))/cmath.exp(k*(Z0-L))
        Z_XY=[ IMP, IMP ]

        x=[ [X1,Z0], [X2,Z0] ]
        eta=None

        z=domain.getX()[1]
        Hx0_ex=cos(k.imag*(z-L))*exp(k.real*(z-L))
        Hx1_ex=sin(k.imag*(z-L))*exp(k.real*(z-L))
        Hx0_ex_z=-sin(k.imag*(z-L))*k.imag*exp(k.real*(z-L))+cos(k.imag*(z-L))*exp(k.real*(z-L))*k.real
        Hx1_ex_z=cos(k.imag*(z-L))*k.imag*exp(k.real*(z-L))+sin(k.imag*(z-L))*exp(k.real*(z-L))*k.real

        model=MT2DModelTMMode(domain, omega, x, Z_XY, eta, mu=mu0, sigma0=SIGMA, tol=1e-9,  directSolver=True,  airLayerLevel=1.)

        args=model.getArguments(RHO)
        Hx=args[0]
        g_Hx=args[1]
        Hxz=g_Hx[:,1]
        self.assertLess(Lsup(Hx[0]-Hx0_ex), 1e-4 * Lsup(Hx0_ex))
        self.assertLess(Lsup(Hx[1]-Hx1_ex), 1e-4 * Lsup(Hx1_ex))
        self.assertLess(Lsup(Hxz[0]-Hx0_ex_z), 1e-2 * Lsup(Hx0_ex_z))
        self.assertLess(Lsup(Hxz[1]-Hx1_ex_z), 1e-2 * Lsup(Hx1_ex_z))

        argsr=model.getArguments(1.)
        ref=model.getDefect(1., *argsr)

        # this should be almost zero:
        args=model.getArguments(RHO)
        d=model.getDefect(RHO, *args)
        self.assertTrue( d > 0.)
        self.assertTrue( ref > 0.)
        self.assertTrue( d <= 3e-3 * ref ) # d should be zero (some sort of)

        z=ReducedFunction(domain).getX()[1]
        Hx0_ex=cos(k.imag*(z-L))*exp(k.real*(z-L))
        Hx1_ex=sin(k.imag*(z-L))*exp(k.real*(z-L))
        Hx0_ex_z=-sin(k.imag*(z-L))*k.imag*exp(k.real*(z-L))+cos(k.imag*(z-L))*exp(k.real*(z-L))*k.real
        Hx1_ex_z=cos(k.imag*(z-L))*k.imag*exp(k.real*(z-L))+sin(k.imag*(z-L))*exp(k.real*(z-L))*k.real
        g_Hx = Data(0, (2,2), Hx0_ex_z.getFunctionSpace())
        g_Hx[0,1] = Hx0_ex_z
        g_Hx[1,1] = Hx1_ex_z
        # and this should be zero
        d0=model.getDefect(RHO, Hx0_ex*[1.,0]+ Hx1_ex*[0,1.], g_Hx)
        self.assertLess( d0, 1e-8 * ref ) # d should be zero (some sort of)

        # and this too
        dg=model.getGradient(RHO, Hx0_ex*[1.,0]+Hx1_ex*[0,1.], g_Hx)
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==())
        self.assertLess(Lsup(dg), 1e-10)

    def test_Differential(self):
        INC=0.001
        omega=5.
        mu0=0.123
        RHO=0.15
        SIGMA=1/RHO
        k=cmath.sqrt(1j*omega*mu0*1/RHO)  # Hx=exp(k*z)

        L=1
        NE=101
        domain=ripRectangle(max(NE,50*mpisize-1), max(NE,50*mpisize-1), d1=mpisize)

        Z0=0.5
        IMP=RHO*k*(cmath.exp(k*(Z0-L))-cmath.exp(-k*(Z0-L)))/(cmath.exp(k*(Z0-L))+cmath.exp(-k*(Z0-L)))
        Z_XY=[ IMP, IMP ]
        H=1./NE
        X1=(int(0.3/H)+0.5)*H
        X2=(int(0.6/H)+0.5)*H
        x=[ [X1,Z0], [X2,Z0] ]
        eta=None

        model=MT2DModelTMMode(domain, omega, x, Z_XY, eta, mu=mu0, tol=1e-9, sigma0=SIGMA, directSolver=True)

        # this is the base line:
        xx=domain.getX()[0]
        RHO0=3.*(xx+0.3)
        args0=model.getArguments(RHO0)
        d0=model.getDefect(RHO0, *args0)
        dg0=model.getGradient(RHO0, *args0)
        self.assertTrue(isinstance(dg0, Data))
        self.assertTrue(dg0.getShape()==())

        X=Function(domain).getX()

        # test 1
        p=INC
        RHO1=RHO0+p
        args1=model.getArguments(RHO1)
        d1=model.getDefect(RHO1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )

        # test 2
        p=exp(-length(X-(0.2,0.2))**2/10)*INC
        RHO1=RHO0+p
        args1=model.getArguments(RHO1)
        d1=model.getDefect(RHO1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )

        # test 3
        p=sin(length(X)*3*3.14)*INC
        RHO1=RHO0+p
        args1=model.getArguments(RHO1)
        d1=model.getDefect(RHO1, *args1)
        self.assertLess( abs( d1-d0-integrate(dg0*p) ), 1e-2*abs(d1-d0) )

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

