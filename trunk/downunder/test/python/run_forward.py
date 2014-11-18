
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
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
from esys.escript import getEscriptParamInt
from esys.escript.pdetools import Locator

mpisize = getMPISizeWorld()
# this is mainly to avoid warning messages
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='ref_data'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

    
have_direct=getEscriptParamInt("PASO_DIRECT")   


@unittest.skipIf(mpisize>1 or have_direct!=1, "more than 1 MPI rank or missing direct solver")
class TestAcousticInversion(unittest.TestCase):
    def test_API(self):
        from esys.ripley import Rectangle
        domain=Rectangle(20,20, diracPoints=[(0.5,1.)], diracTags=['sss'])
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
         
        from esys.ripley import Rectangle
        domain=Rectangle(100,100, diracPoints=[(0.5,1.)], diracTags=['sss'])
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
        # 
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
        self.assertTrue(Lsup(dg) < 1e-10)

        # this shuld be zero' too
        sigma_comps=[2*sigma.real, sigma.imag/2.]
        args=acw.getArguments(sigma_comps)
        d=acw.getDefect(sigma_comps, *args)
        self.assertTrue(isinstance(d, float))
        self.assertTrue(abs(d)< 1e-10)
        
        dg=acw.getGradient(sigma_comps, *args)
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==(2,))
        self.assertTrue(dg.getFunctionSpace()==Solution(domain))                 
        self.assertTrue(Lsup(dg) < 1e-10)
        
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
        self.assertTrue( abs( d1-d0-integrate(dg0[0]*p) ) < 1e-2  * abs(d1-d0) )

        sigma2=sigma0+p*[0,1]
        args=acw.getArguments(sigma2)
        d2=acw.getDefect(sigma2, *args)
        self.assertTrue( abs(d2-d0-integrate(dg0[1]*p))  < 1e-2  * abs(d2-d0) )
    
    def test_numeric2DnoscaleF(self):
         
        from esys.ripley import Rectangle
        domain=Rectangle(10,20, diracPoints=[(0.5,1.)], diracTags=['sss'])
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
        self.assertTrue(Lsup(d) < 1e-10)
        #self.assertTrue(d >= 0)
        #self.assertTrue(d < 1e-10)
        
        dg=acw.getGradient(sigma_comps, *args)

        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==(2,))
        self.assertTrue(dg.getFunctionSpace()==Solution(domain))                 
        self.assertTrue(Lsup(dg) < 5e-10)
        # this shouldn't be zero:
        sigma0=Data([2*sigma.real, sigma.imag/2], Function(domain) )
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
        h=0.001
        X=Function(domain).getX()
        p=h*sin(length(X)*np.pi)*Lsup(length(sigma0))
        
        sigma1=sigma0+p*[1,0]
        args=acw.getArguments(sigma1)
        d1=acw.getDefect(sigma1, *args)

        self.assertTrue( abs( d1-d0-integrate(dg0[0]*p) ) < 1e-2  * abs(d1-d0) )

        sigma2=sigma0+p*[0,1]
        args=acw.getArguments(sigma2)
        d2=acw.getDefect(sigma2, *args)
        self.assertTrue( abs(d2-d0-integrate(dg0[1]*p))  < 1e-2  * abs(d2-d0) )

class TestMT2DModelTEMode(unittest.TestCase):
    def test_API(self):
        from esys.ripley import Rectangle
        domain=Rectangle(19, 19, d1=mpisize)
        omega=2.
        x=[ [0.2,0.5], [0.3,0.5] ]
        Z_XY=[ complex(1.2,1.5), complex(1.3,2.5) ]
        eta=1.
        w0=1.
        E_x0=1.
        # now we do a real one
        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta, w0=w0, E_x0=E_x0)
        self.assertEqual(acw.getDomain(),  domain)
        pde=acw.setUpPDE()
        self.assertIsInstance(pde, LinearPDE)
        self.assertEqual(pde.getNumEquations(), 2)
        self.assertEqual(pde.getNumSolutions(), 2)
        self.assertEqual(pde.getDomain(),  domain)

        # other things that should work 
        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta=[1.,1.], w0=[2.,3.], E_x0=complex(4.5,6) )

        # these shouldn't work
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, x, [3.], eta=[1.,1.], w0=[2.,3.], E_x0=complex(4.5,6) )
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, x, Z_XY, eta=[1.], w0=[2.,3.], E_x0=complex(4.5,6) )
        self.assertRaises(ValueError, MT2DModelTEMode, domain, omega, [(6.7,5)], Z_XY, eta=[1.,1.], w0=[2.,3.], E_x0=complex(4.5,6) )

    def test_PDE(self):
        omega=2.
        mu0=0.123
        SIGMA=15.
        k=cmath.sqrt(1j*omega*mu0*SIGMA)  # Ex=exp(k*z)

        from esys.ripley import Rectangle
        domain=Rectangle(199,199, d1=mpisize)

        
        IMP=cmath.sqrt(1j*omega*mu0/SIGMA)
        Z_XY=[ IMP, IMP ]
        x=[ [0.3,0.5], [0.6,0.5] ]
        eta=0.005
        z=domain.getX()[1]
        Ex0_ex=exp(-k.real*z)*cos(-k.imag*z)
        Ex0_ex_z=(-k.real*cos(-k.imag*z)+k.imag*sin(-k.imag*z)) * exp(-k.real*z)
        Ex1_ex=exp(-k.real*z)*sin(-k.imag*z)
        Ex1_ex_z=(-k.real*sin(-k.imag*z)-k.imag*cos(-k.imag*z)) * exp(-k.real*z)

        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, fixAtBottom=True, E_x0=Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], tol=1e-9)

        args=acw.getArguments(SIGMA)
        Ex=args[0]
        Exz=args[1]        
        self.assertTrue(Lsup(Ex[0]-Ex0_ex) <= 1e-4 * Lsup(Ex0_ex))
        self.assertTrue(Lsup(Ex[1]-Ex1_ex) <= 1e-4 * Lsup(Ex1_ex))
        self.assertTrue(Lsup(Exz[0]-Ex0_ex_z) <= 1e-2 * Lsup(Ex0_ex_z))
        self.assertTrue(Lsup(Exz[1]-Ex1_ex_z) <= 1e-2 * Lsup(Ex1_ex_z))

        argsr=acw.getArguments(0.)
        ref=acw.getDefect(0., *argsr)
        
        # this should be almost zero:
        args=acw.getArguments(SIGMA)
        d=acw.getDefect(SIGMA, *args)
        
        self.assertTrue( d > 0.)
        self.assertTrue( ref > 0.)
        self.assertTrue( d <= 1e-4 * ref ) # d should be zero (some sort of)  

        # and this should be zero        
        d0=acw.getDefect(SIGMA, Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], Ex0_ex_z*[1.,0]+ Ex1_ex_z*[0,1.])
        self.assertTrue( d0 <= 1e-8 * ref ) # d should be zero (some sort of)  
        

        # and this too
        dg=acw.getGradient(SIGMA, Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], Ex0_ex_z*[1.,0]+ Ex1_ex_z*[0,1.])
        self.assertTrue(isinstance(dg, Data))
        self.assertTrue(dg.getShape()==())              
        self.assertTrue(Lsup(dg) < 1e-10)
        
    def test_Differential(self):
    
        INC=0.01 
        
        omega=2.
        mu0=0.123
        SIGMA=15.
        k=cmath.sqrt(1j*omega*mu0*SIGMA)  # Ex=exp(k*z)

        from esys.ripley import Rectangle
        domain=Rectangle(99,99, d1=mpisize)

        
        IMP=-cmath.sqrt(1j*omega*mu0/SIGMA)
        Z_XY=[ IMP, IMP ]
        x=[ [0.3,0.5], [0.6,0.5] ]
        eta=0.005
        z=domain.getX()[1]
        Ex0_ex=exp(-k.real*z)*cos(-k.imag*z)
        Ex0_ex_z=(-k.real*cos(-k.imag*z)+k.imag*sin(-k.imag*z)) * exp(-k.real*z)
        Ex1_ex=exp(-k.real*z)*sin(-k.imag*z)
        Ex1_ex_z=(-k.real*sin(-k.imag*z)-k.imag*cos(-k.imag*z)) * exp(-k.real*z)

        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, fixAtBottom=True, E_x0=Ex0_ex*[1.,0]+ Ex1_ex*[0,1.], tol=1e-9 )
    
        # this is the base line:
        SIGMA0=10.        
        args0=acw.getArguments(SIGMA0)
        d0=acw.getDefect(SIGMA0, *args0)
        
        dg0=acw.getGradient(SIGMA0, *args0)
        self.assertTrue(isinstance(dg0, Data))
        self.assertTrue(dg0.getShape()==())
        
        
        X=Function(domain).getX()

        # test 1
        p=INC
        SIGMA1=SIGMA0+p
        args1=acw.getArguments(SIGMA1)
        d1=acw.getDefect(SIGMA1, *args1)
        self.assertTrue( abs( d1-d0-integrate(dg0*p) ) < 1e-2  * abs(d1-d0) )
                
        # test 2
        p=exp(-length(X-(0.2,0.2))**2/10)*INC
        SIGMA1=SIGMA0+p
        args1=acw.getArguments(SIGMA1)
        d1=acw.getDefect(SIGMA1, *args1)
        self.assertTrue( abs( d1-d0-integrate(dg0*p) ) < 1e-2  * abs(d1-d0) )

        # test 3
        p=sin(length(X)*3*3.14)*INC
        SIGMA1=SIGMA0+p
        args1=acw.getArguments(SIGMA1)
        d1=acw.getDefect(SIGMA1, *args1)
        self.assertTrue( abs( d1-d0-integrate(dg0*p) ) < 1e-2  * abs(d1-d0) )

class TestSubsidence(unittest.TestCase):
    def test_PDE(self):
         
        lam=2.
        mu=1.

        from esys.ripley import Brick
        domain=Brick(20,20,19, d2=mpisize)
        
        xb=FunctionOnBoundary(domain).getX()
        m=whereZero(xb[2]-1)
        w=m*[0,0,1]
        d=m*2.5
        acw=Subsidence(domain, w,d, lam, mu )
        
        P0=10.        
        args0=acw.getArguments(P0)
        u=args0[0]
        self.assertTrue(Lsup(u[0]) < 1.e-8)
        self.assertTrue(Lsup(u[1]) < 1.e-8)
        self.assertTrue(Lsup(u[2]-2.5*domain.getX()[2]) < 1.e-8)
        
        dd=acw.getDefect(P0, *args0)
        
        self.assertTrue( dd >= 0.)
        self.assertTrue( dd <= 1e-7 * 2.5 )
    def test_Differential(self):
    
        lam=2.
        mu=1.
        
        INC=0.01 
        from esys.ripley import Brick
        domain=Brick(20,20,20*mpisize-1 , d2=mpisize)
        
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
        self.assertTrue(abs((d1-d0)/INC-integrate(grad_d* dP)) < ref * 1.e-5) 

        dP=exp(-(length(x-[0.3,0.3,0.5])/0.06)**2)
        P2=P0-INC*dP
        args2=acw.getArguments(P2)
        d2=acw.getDefect(P2, *args2)
        ref=abs((d2-d0)/INC)
        self.assertTrue(abs((d2-d0)/INC+integrate(grad_d* dP)) < ref * 1.e-5) 

class TestDCResistivity(unittest.TestCase):

    def test_PDE2D(self):
         
        dx_tests=0.1 

        sigma0=1.
        electrodes=[(0.5-2*dx_tests,1.), (0.5-dx_tests,1.), (0.5+dx_tests,1.), (0.5+2*dx_tests,1.)]
        from esys.finley import Rectangle
        domain=Rectangle(20,20, d1=mpisize,  diracPoints=electrodes, diracTags=["sl0", "sl1", "sr0", "sr1"] )
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

        self.assertTrue(Lsup(phiPrimary-acw.getPrimaryPotential()) < 1.e-10 * Lsup(acw.getPrimaryPotential()))

        SIGMA=10. # matches current        
        args0=acw.getArguments(SIGMA)
        p=args0[0]
        u=args0[1]

        # true secondary potential
        pps=pp-phiPrimary
        self.assertTrue(Lsup(p-pps) < 1.e-6 * Lsup(pps))


        # test return values at electrodes:
        self.assertTrue(abs(u[0]-uu[0]*uuscale) < 1.e-6 * abs(uu[0]*uuscale))
        self.assertTrue(abs(u[1]-uu[1]*uuscale) < 1.e-6 * abs(uu[1]*uuscale))

        # this sould be zero
        dd=acw.getDefect(SIGMA, *args0)
        self.assertTrue( dd >= 0.)
        self.assertTrue( dd <= 1e-7 )

    def test_Differential2D(self):

        INC=0.001

        sigma0=1.
        dx_tests=0.1 
        electrodes=[(0.5-2*dx_tests,1.), (0.5-dx_tests,1.), (0.5+dx_tests,1.), (0.5+2*dx_tests,1.)]
        from esys.finley import Rectangle
        domain=Rectangle(20,20, d1=mpisize,  diracPoints=electrodes, diracTags=["sl0", "sl1", "sr0", "sr1"] )
        loc=Locator(domain,electrodes[2:])

        # arguments for DcRes
        #current = 10.
        sampleTags = [ ("sr0", "sr1") ]

        delphi_in = [ 0.05 ]

        sigmaPrimary=1
        x=domain.getX()
        phiPrimary=(x[0]-inf(x[0]))*(x[1]-inf(x[1]))*(x[0]-sup(x[0]))

        acw=DcRes(domain, loc, delphi_in, sampleTags,  phiPrimary, sigmaPrimary)

        #===========================================================================
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
        self.assertTrue(abs((d1-d0)/INC-integrate(grad_d* dS)) < ref * 1.e-3) 

        dS=-exp(-(length(x-[0.5,0.5])/0.2)**2)
        SIGMA2=SIGMA0+INC*dS
        args2=acw.getArguments(SIGMA2)
        d2=acw.getDefect(SIGMA2, *args2)
        ref=abs((d2-d0)/INC)
        self.assertTrue(abs((d2-d0)/INC-integrate(grad_d* dS)) < ref * 1.e-3) 

        dS=-1
        SIGMA3=SIGMA0+INC*dS
        args3=acw.getArguments(SIGMA3)
        d3=acw.getDefect(SIGMA3, *args3)
        ref=abs((d3-d0)/INC)
        self.assertTrue(abs((d3-d0)/INC-integrate(grad_d* dS)) < ref * 1.e-3) 

        dS=1
        SIGMA4=SIGMA0+INC*dS
        args4=acw.getArguments(SIGMA4)
        d4=acw.getDefect(SIGMA4, *args4)
        ref=abs((d4-d0)/INC)
        self.assertTrue(abs((d4-d0)/INC-integrate(grad_d* dS)) < ref * 1.e-3) 

class TestIsostaticPressure(unittest.TestCase):
    def test_all(self):
        from esys.ripley import Brick
        domain=Brick(50,50,20*mpisize-1, d2=mpisize)

        ps=IsostaticPressure(domain, level0=1., coordinates=None)
    
        g=Vector(0., Function(domain))
        rho=Scalar(100, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*981.
        self.assertTrue(Lsup(p0-p_ref) < 1e-6 * Lsup(p_ref))

        g=Vector([0,0,-10], Function(domain))
        rho=Scalar(0, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*26700
        self.assertTrue(Lsup(p0-p_ref) < 1e-6 * Lsup(p_ref))

        g=Vector([0,0,-10], Function(domain))
        rho=Scalar(100, Function(domain))
        p0=ps.getPressure(g, rho)
        p_ref=-(1.-domain.getX()[2])*(981.+26700+1000)
        self.assertTrue(Lsup(p0-p_ref) < 1e-6 * Lsup(p_ref))
                               
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
    