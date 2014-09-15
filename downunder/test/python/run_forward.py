
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
        domain=Rectangle(20,20)
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
        domain=Rectangle(200,200)

        
        IMP=cmath.sqrt(1j*omega*mu0/SIGMA)
        Z_XY=[ IMP, IMP ]
        x=[ [0.3,0.5], [0.6,0.5] ]
        eta=0.005
        z=domain.getX()[1]
        Ex0_ex=exp(-k.real*z)*cos(-k.imag*z)
        Ex0_ex_z=(-k.real*cos(-k.imag*z)+k.imag*sin(-k.imag*z)) * exp(-k.real*z)
        Ex1_ex=exp(-k.real*z)*sin(-k.imag*z)
        Ex1_ex_z=(-k.real*sin(-k.imag*z)-k.imag*cos(-k.imag*z)) * exp(-k.real*z)

        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, fixAtBottom=True, E_x0=Ex0_ex*[1.,0]+ Ex1_ex*[0,1.] )

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
        domain=Rectangle(200,200)

        
        IMP=-cmath.sqrt(1j*omega*mu0/SIGMA)
        Z_XY=[ IMP, IMP ]
        x=[ [0.3,0.5], [0.6,0.5] ]
        eta=0.005
        z=domain.getX()[1]
        Ex0_ex=exp(-k.real*z)*cos(-k.imag*z)
        Ex0_ex_z=(-k.real*cos(-k.imag*z)+k.imag*sin(-k.imag*z)) * exp(-k.real*z)
        Ex1_ex=exp(-k.real*z)*sin(-k.imag*z)
        Ex1_ex_z=(-k.real*sin(-k.imag*z)-k.imag*cos(-k.imag*z)) * exp(-k.real*z)

        acw=MT2DModelTEMode(domain, omega, x, Z_XY, eta, mu=mu0, fixAtBottom=True, E_x0=Ex0_ex*[1.,0]+ Ex1_ex*[0,1.] )
    
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
        
                                  
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
    
