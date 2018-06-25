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

"""Forward model for Subsidence modelling"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['Subsidence']

from .base import ForwardModel
from esys.escript import Data, FunctionOnBoundary
from esys.escript.linearPDEs import LinearPDESystem
from esys.escript.util import *


class Subsidence(ForwardModel):
    """
    Forward Model for subsidence inversion minimizing
    integrate( (inner(w,u)-d)**2)
    where u is the surface displacement due to a pressure change P
    """
    def __init__(self, domain, w, d, lam, mu, coordinates=None, tol=1e-8):
        """
        Creates a new subsidence on the given domain

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors and direction
        :type w: ``Vector`` with ``FunctionOnBoundary``
        :param d: displacement measured at surface
        :type d: ``Scalar`` with ``FunctionOnBoundary``
        :param lam: 1st Lame coefficient
        :type lam: ``Scalar`` with ``Function``
        :param lam: 2st Lame coefficient/Shear modulus
        :type lam: ``Scalar`` with ``Function``
        :param coordinates: defines coordinate system to be used (not supported yet))
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        """
        super(Subsidence, self).__init__()
        DIM=domain.getDim()

        self.__pde=LinearPDESystem(domain)
        self.__pde.setSymmetryOn()
        self.__pde.getSolverOptions().setTolerance(tol)
        #... set coefficients ...
        C=self.__pde.createCoefficient('A')
        for i in range(DIM):
            for j in range(DIM):
                C[i,i,j,j]+=lam
                C[i,j,i,j]+=mu
                C[i,j,j,i]+=mu
        x=domain.getX()
        msk=whereZero(x[DIM-1])*kronecker(DIM)[DIM-1]
        for i in range(DIM-1):
            xi=x[i]
            msk+=(whereZero(xi-inf(xi))+whereZero(xi-sup(xi))) *kronecker(DIM)[i]
        self.__pde.setValue(A=C,q=msk)

        self.__w=interpolate(w, FunctionOnBoundary(domain))
        self.__d=interpolate(d, FunctionOnBoundary(domain))

    def rescaleWeights(self, scale=1., P_scale=1.):
        """
        rescales the weights
        
        :param scale: scale of data weighting factors
        :type scale: positive ``float``
        :param P_scale: scale of pressure increment
        :type P_scale: ``Scalar``
        """
        pass

    def getArguments(self, P):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param P: pressure
        :type P: ``Scalar``
        :return: displacement u
        :rtype: ``Vector``
        """
        DIM=self.__pde.getDim()
        self.__pde.setValue(y=Data(),X=P*kronecker(DIM))
        u= self.__pde.getSolution()
        return u,

    def getDefect(self, P,u):
        """
        Returns the value of the defect.

        :param P: pressure
        :type P: ``Scalar``
        :param u: corresponding displacement
        :type u: ``Vector``
        :rtype: ``float``
        """
        return 0.5*integrate((inner(u,self.__w)-self.__d)**2)

    def getGradient(self, P, u):
        """
        Returns the gradient of the defect with respect to susceptibility.

        :param P: pressure
        :type P: ``Scalar``
        :param u: corresponding displacement
        :type u: ``Vector``
        :rtype: ``Scalar``
        """
        d=inner(u,self.__w)-self.__d
        self.__pde.setValue(y=d*self.__w,X=Data())
        ustar=self.__pde.getSolution()

        return div(ustar)

