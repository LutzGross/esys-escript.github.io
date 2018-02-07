
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

"""
Generic base class for PDE solving tests
"""

from esys.escript import Data, Function, Lsup, Solution, Tensor4, Vector, \
                         grad, inner, kronecker, matrixmult, whereZero
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import esys.escriptcore.utestselect as unittest
import numpy


class SimpleSolveTestCase(unittest.TestCase):
    REL_TOL = 1.e-6
    SOLVER_VERBOSE = False
    SOLVER_TOL = 1.e-8
    FAC_DIAG = 1.
    FAC_OFFDIAG = -0.4
    # the following members must be set by the test methods in subclasses
    domain = None
    package = None
    method = None
    preconditioner = None

    def _getGrad(self, system):
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

    def _getSolution(self, system):
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

    def _setCoefficients(self, pde, system):
        """sets PDE coefficients"""
        FAC_DIAG = 1.
        FAC_OFFDIAG = -0.4
        x = Solution(self.domain).getX()
        mask = whereZero(x[0])
        dim = self.domain.getDim()
        u_ex = self._getSolution(system)
        g_ex = self._getGrad(system)

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

    def _setSolverOptions(self, so):
        """override this to modify solver options prior to solving"""
        pass

    def getPDE(self, system):
        dim = self.domain.getDim()
        if system:
            pde=LinearPDE(self.domain, numEquations=dim)
        else:
            pde=LinearPDE(self.domain, numEquations=1)

        self._setCoefficients(pde, system)
        so = pde.getSolverOptions()
        so.setPackage(self.package)
        so.setSolverMethod(self.method)
        so.setPreconditioner(self.preconditioner)
        so.setTolerance(self.SOLVER_TOL)
        so.setVerbosity(self.SOLVER_VERBOSE)
        self._setSolverOptions(so)
        return pde, self._getSolution(system), self._getGrad(system)

    def test_single(self):
        pde, u_ex, g_ex = self.getPDE(False)
        g=grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

    def test_system(self):
        pde, u_ex, g_ex = self.getPDE(True)
        g = grad(u_ex)
        self.assertLess(Lsup(g_ex-g), self.REL_TOL*Lsup(g_ex))
        u = pde.getSolution()
        error = Lsup(u-u_ex)
        self.assertLess(error, self.REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)


