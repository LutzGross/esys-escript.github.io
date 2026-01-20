
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Generic base class for PDE solving tests
"""

from esys.escript import Data, Function, Lsup, Solution, Tensor4, Vector, \
                         grad, inner, kronecker, matrixmult, whereZero, hasFeature
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import esys.escriptcore.utestselect as unittest
import numpy

HAVE_DIRECT_PASO = hasFeature('paso') and (hasFeature('umfpack') or hasFeature("mkl") or hasFeature("mumps"))
HAVE_MUMPS = hasFeature("mumps")
HAVE_TRILINOS = hasFeature('trilinos')
HAVE_SOLVER = HAVE_DIRECT_PASO or HAVE_TRILINOS
HAVE_SOLVER_COMPLEX = HAVE_TRILINOS or HAVE_MUMPS

class SolveTestCaseTemplate(unittest.TestCase):
    """
    this is the template class for testing solvers:
    """
    REL_TOL = 1.e-6
    SOLVER_VERBOSE = False
    SOLVER_TOL = 1.e-8

    # the following members must be set by the test methods in subclasses
    domain = None
    package = None
    method = None
    preconditioner = SolverOptions.NO_PRECONDITIONER


    def getPDE(self, system, iscomplex=False):
        dim = self.domain.getDim()
        if system:
            pde=LinearPDE(self.domain, numEquations=dim, isComplex=iscomplex)
        else:
            pde=LinearPDE(self.domain, numEquations=1, isComplex=iscomplex)

        self.setCoefficients(pde, system)
        so = pde.getSolverOptions()
        so.setPackage(self.package)
        so.setSolverMethod(self.method)
        so.setPreconditioner(self.preconditioner)
        so.setTolerance(self.SOLVER_TOL)
        so.setVerbosity(self.SOLVER_VERBOSE)
        pde.setSolverOptions(so)
        return pde, self.getSolution(system), self.getGrad(system)

    

class SolveTestCaseOrder1(SolveTestCaseTemplate):
    """
    this is the class for testing solvers for order 1 meshes:
    """

    def getGrad(self, system):
        """returns exact gradient"""
        dim = self.domain.getDim()
        if system:
            g_ex = Data(0., (dim,dim), Solution(self.domain))
            if dim == 2:
                g_ex[0,0] = 2.
                g_ex[0,1] = 3.
                g_ex[1,0] = 3.
                g_ex[1,1] = 2.
            else:
                g_ex[0,0] = 2.
                g_ex[0,1] = 3.
                g_ex[0,2] = 4.
                g_ex[1,0] = 4.
                g_ex[1,1] = 1.
                g_ex[1,2] = -2.
                g_ex[2,0] = 8.
                g_ex[2,1] = 4.
                g_ex[2,2] = 5.
        else:
            g_ex = Data(0., (dim,), Solution(self.domain))
            if dim == 2:
                g_ex[0] = 2.
                g_ex[1] = 3.
            else:
                g_ex[0] = 2.
                g_ex[1] = 3.
                g_ex[2] = 4.
        return g_ex

    def getSolution(self, system):
        """returns exact solution"""
        dim = self.domain.getDim()
        x = Solution(self.domain).getX()
        if system:
            u_ex = Vector(0., Solution(self.domain))
            if dim == 2:
                u_ex[0] =  1.+2.*x[0]+3.*x[1]
                u_ex[1] = -1.+3.*x[0]+2.*x[1]
            else:
                u_ex[0] =  1.+2.*x[0]+3.*x[1]+4.*x[2]
                u_ex[1] = -1.+4.*x[0]+1.*x[1]-2.*x[2]
                u_ex[2] =  5.+8.*x[0]+4.*x[1]+5.*x[2]
        else:
            if dim == 2:
                u_ex = 1.+2.*x[0]+3.*x[1]
            else:
                u_ex = 1.+2.*x[0]+3.*x[1]+4.*x[2]
        return u_ex
    
    def setCoefficients(self, pde, system):
        """sets PDE coefficients"""
        FAC_DIAG = self.FAC_DIAG
        FAC_OFFDIAG =self.FAC_OFFDIAG
        x = Solution(self.domain).getX()
        mask = whereZero(x[0])
        dim = self.domain.getDim()
        u_ex = self.getSolution(system)
        g_ex = self.getGrad(system)

        if system:
            A = Tensor4(0., Function(self.domain))
            for i in range(dim):
                A[i,:,i,:] = kronecker(dim)

            Y = Vector(0., Function(self.domain))
            if dim == 2:
                Y[0] = u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
                Y[1] = u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
            else:
                Y[0] = u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
                Y[1] = u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
                Y[2] = u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
            pde.setValue(r=u_ex, q=mask*numpy.ones(dim,),
                         A=A,
                         D=kronecker(dim)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((dim,dim))*FAC_OFFDIAG,
                         Y=Y,
                         y=matrixmult(g_ex,self.domain.getNormal()))
        else:
            pde.setValue(r=u_ex, q=mask, A=kronecker(dim),
                         y=inner(g_ex, self.domain.getNormal()))

    
class SolveTestCaseOrder2(SolveTestCaseTemplate):
    """
    this is the class for testing solvers for order 2 meshes:
    """
    
    def getGrad(self, system):
        """returns exact gradient"""
        dim = self.domain.getDim()
        x = Solution(self.domain).getX()
        if system:
            g_ex = Data(0., (dim,dim), Solution(self.domain))
            if dim == 2:
                g_ex[0,0] = 2.+8.*x[0]+ 5.*x[1]
                g_ex[0,1] = 3.+5.*x[0]+12.*x[1]
                g_ex[1,0] = 4.+2.*x[0]+ 6.*x[1]
                g_ex[1,1] = 2.+6.*x[0]+ 8.*x[1]
            else:
                g_ex[0,0] =  2.+6.*x[1]+8.*x[2]+18.*x[0]
                g_ex[0,1] =  3.+6.*x[0]+7.*x[2]+20.*x[1]
                g_ex[0,2] =  4.+7.*x[1]+8.*x[0]+22.*x[2]
                g_ex[1,0] =  4.+3.*x[1]-8.*x[2]- 4.*x[0]
                g_ex[1,1] =  1.+3.*x[0]+2.*x[2]+14.*x[1]
                g_ex[1,2] = -6.+2.*x[1]-8.*x[0]+10.*x[2]
                g_ex[2,0] =  7.-6.*x[1]+2.*x[2]+ 4.*x[0]
                g_ex[2,1] =  9.-6.*x[0]+8.*x[2]+16.*x[1]
                g_ex[2,2] =  2.+8.*x[1]+2.*x[0]+ 2.*x[2]
        else:
            g_ex = Data(0., (dim,), Solution(self.domain))
            if dim == 2:
                g_ex[0] = 2.+8.*x[0]+5.*x[1]
                g_ex[1] = 3.+5.*x[0]+12.*x[1]
            else:
                g_ex[0] = 2.+6.*x[1]+8.*x[2]+18.*x[0]
                g_ex[1] = 3.+6.*x[0]+7.*x[2]+20.*x[1]
                g_ex[2] = 4.+7.*x[1]+8.*x[0]+22.*x[2]
        return g_ex

    def getSolution(self, system):
        """returns exact solution"""
        dim = self.domain.getDim()
        x = Solution(self.domain).getX()
        if system:
            u_ex = Vector(0., Solution(self.domain))
            if dim == 2:
                u_ex[0] =  1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
                u_ex[1] = -1.+4.*x[0]+2.*x[1]+1.*x[0]**2+6.*x[1]*x[0]+4.*x[1]**2
            else:
                u_ex[0] = 1.+2.*x[0]+3.*x[1]+4.*x[2]+\
                          6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+\
                          9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
                u_ex[1] = 2.+4.*x[0]+1.*x[1]-6.*x[2]+\
                          3.*x[0]*x[1]+2.*x[1]*x[2]-8.*x[2]*x[0]-\
                          2.*x[0]**2+7.*x[1]**2+5.*x[2]**2
                u_ex[2] = -2.+7.*x[0]+9.*x[1]+2*x[2]-\
                          6.*x[0]*x[1]+8.*x[1]*x[2]+2.*x[2]*x[0]+\
                          2.*x[0]**2+8.*x[1]**2+1.*x[2]**2
        else:
            if dim == 2:
                u_ex = 1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
            else:
                u_ex = 1.+2.*x[0]+3.*x[1]+4.*x[2]+\
                       6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+\
                       9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        return u_ex

    def setCoefficients(self, pde, system):
        """sets PDE coefficients"""        
        FAC_DIAG = self.FAC_DIAG
        FAC_OFFDIAG =self.FAC_OFFDIAG
        x = Solution(self.domain).getX()
        mask = whereZero(x[0])
        dim = self.domain.getDim()
        u_ex = self.getSolution(system)
        g_ex = self.getGrad(system)

        if system:
            A = Tensor4(0., Function(self.domain))
            for i in range(dim):
                A[i,:,i,:] = kronecker(dim)

            Y = Vector(0., Function(self.domain))
            if dim == 2:
                Y[0] = u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG-20
                Y[1] = u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG-10
            else:
                Y[0] = u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG-60
                Y[1] = u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG-20
                Y[2] = u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG-22
            pde.setValue(r=u_ex, q=mask*numpy.ones(dim,),
                         A=A,
                         D=kronecker(dim)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((dim,dim))*FAC_OFFDIAG,
                         Y=Y,
                         y=matrixmult(g_ex,self.domain.getNormal()))
        else:
            pde.setValue(r=u_ex, q=mask, A=kronecker(dim),
                         y=inner(g_ex, self.domain.getNormal()))
            if dim == 2:
                pde.setValue(Y=-20.)
            else:
                pde.setValue(Y=-60.)


class SimpleSolveTestCase(SolveTestCaseOrder1):
    """
    testing the real PDEs 
    """
    FAC_DIAG = 1.
    FAC_OFFDIAG = -0.4
    def test_single(self):
        pde, u_ex, g_ex = self.getPDE(False)
        g=grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))

        u = pde.getSolution()
        self.assertFalse(u.isComplex())
        self.assertEqual(u.getShape(), ( ))
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

    @unittest.skipIf(not HAVE_SOLVER, "No solver available")
    def test_system(self):
        pde, u_ex, g_ex = self.getPDE(True)
        g = grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        self.assertFalse(u.isComplex())
        self.assertEqual(u.getShape(), (pde.getDim(), ))
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolveTestCaseOrder2(SolveTestCaseOrder2):
    """
    testing the real PDEs 
    """
    FAC_DIAG = 1.
    FAC_OFFDIAG = -0.4
    def test_single(self):
        pde, u_ex, g_ex = self.getPDE(False)
        g=grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        self.assertFalse(u.isComplex())
        self.assertEqual(u.getShape(), ( ))
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

    @unittest.skipIf(not HAVE_SOLVER, "No solver available")
    def test_system(self):
        pde, u_ex, g_ex = self.getPDE(True)
        g = grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        self.assertFalse(u.isComplex())
        self.assertEqual(u.getShape(), (pde.getDim(), ))
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class ComplexSolveTestCase(SolveTestCaseOrder1):
    """
    testing the complex PDEs 
    """
    FAC_DIAG = 1.+0.2j
    FAC_OFFDIAG = -0.4

    @unittest.skipIf(not HAVE_SOLVER_COMPLEX, "No solver available")
    def test_singlecomplex(self):
        pde, u_ex, g_ex = self.getPDE(False, iscomplex=True)
        g=grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        self.assertTrue(u.isComplex())
        self.assertEqual(u.getShape(), ())
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

    @unittest.skipIf(not HAVE_SOLVER_COMPLEX, "No solver available")
    def test_systemcomplex(self):
        pde, u_ex, g_ex = self.getPDE(True, iscomplex=True)
        g = grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        
        self.assertTrue(u.isComplex())
        self.assertEqual(u.getShape(), (pde.getDim(),))
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class ComplexSolveTestCaseOrder2(SolveTestCaseOrder2):
    """
    testing the complex  PDEs for order 2 meshes 
    """
    FAC_DIAG = 1.+0.2j
    FAC_OFFDIAG = -0.4

    @unittest.skipIf(not HAVE_SOLVER_COMPLEX, "No solver available")
    def test_singlecomplex(self):
        pde, u_ex, g_ex = self.getPDE(False, iscomplex=True)
        g=grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        self.assertTrue(u.isComplex())
        self.assertEqual(u.getShape(), ( ))
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

    @unittest.skipIf(not HAVE_SOLVER_COMPLEX, "No solver available")
    def test_systemcomplex(self):
        pde, u_ex, g_ex = self.getPDE(True, iscomplex=True)
        g = grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        self.assertTrue(u.isComplex())
        self.assertEqual(u.getShape(), (pde.getDim(), ))
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
