
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['ForwardModel','GravityModel']

from esys.escript import unitsSI as U
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE

PI = 3.14159265358979323846
G = 6.6742e-11*U.m**3/(U.kg*U.sec**2)

class ForwardModel(object):
    """
    An abstract forward model that can be plugged into a cost function.
    Subclasses should implement getValue() and getGradient().
    """
    def __init__(self):
        pass

    def getArguments(self, x):
        return None

    def getValue(self, x, *args):
        pass

    def getGradient(self, x, *args):
        pass

class GravityModel(ForwardModel):
    """
    Forward Model for gravity inversion as described in the inversion cookbook.
    """
    def __init__(self, domain, chi, g, constrain_top=True, gravity_constant=G, tol=1e-8):
        """
        Creates a new gravity model on the given domain with one or more
        surveys (chi, g).
        """
        self.__domain = domain
        try:
            n=len(chi)
            m=len(g)
            if m != n:
                raise ValueError("Length of chi and g must be the same.")
            self.__chi = chi
            self.__g = g
        except TypeError:
            self.__chi = [chi]
            self.__g = [g]

        A=0
        for s in xrange(len(self.__chi)):
            A2 = integrate(inner(self.__chi[s], self.__g[s]**2))
            if A2 < 0:
                raise ValueError("Negative weighting factor for survey %s"%s)
            A=max(A2, A)
        if not A > 0:
            raise ValueError("No reference data set.")
        self.__G = gravity_constant
        BX = boundingBox(domain)
        DIM = domain.getDim()
        x = domain.getX()
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        if constrain_top:
            constraint=whereZero(x[DIM-1]-BX[DIM-1][1])+whereZero(x[DIM-1]-BX[DIM-1][0])
        else:
            constraint=whereZero(x[DIM-1]-BX[DIM-1][0])
        for i in xrange(DIM-1):
            constraint=constraint+whereZero(x[i]-BX[i][1])+whereZero(x[i]-BX[i][0])
        self.__pde.setValue(A=kronecker(DIM), q=constraint)

    def getSurvey(self, index=None):
        """
        Returns the pair (g_index, chi_index), where g_i is the gravity
        anomaly of survey i, chi_i is the weighting factor for survey i.
        If index is None, all surveys will be returned in a pair of lists.
        """
        if index is None:
            return self.__g, self.__chi
        if index>=len(self.__g):
            raise IndexError("Forward model only has %d surveys"%len(self.__g))
        return self.__g[index], self.__chi[index]

    def getArguments(self, rho):
        """
        Returns precomputed values shared by getValue() and getGradient().
        """
        phi = self.getPotential(rho)
        return phi, -grad(phi)

    def getPotential(self, rho):
        self.__pde.setValue(Y=(-4*PI*self.__G) * rho, X=Data())
        phi=self.__pde.getSolution()
        return phi

    def getValue(self, rho, phi, gravity_force):
        A=0.
        for s in xrange(len(self.__chi)):
            A = inner(self.__chi[s], (gravity_force-self.__g[s])**2) + A
        return 0.5*integrate(A)

    def getGradient(self, rho, phi, gravity_force):
        Z=0.
        for s in xrange(len(self.__chi)):
            Z = self.__chi[s] * (-gravity_force+self.__g[s]) + Z

        self.__pde.setValue(Y=Data(), X=Z)
        ZT=self.__pde.getSolution()
        return ZT*(-4*PI*self.__G)

