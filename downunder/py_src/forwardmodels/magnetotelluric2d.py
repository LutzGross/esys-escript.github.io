from __future__ import print_function
from __future__ import division
##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
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

"""Forward models for 2D MT (TE and TM mode)"""

__copyright__="""Copyright (c) 2003-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['MT2DModelTEMode']

from .base import ForwardModel
from esys.downunder.coordinates import makeTranformation
from esys.escript import Data, Scalar, Vector, Function, FunctionOnBoundary, Solution
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.util import *
from math import pi as PI


class MT2DModelTEMode(ForwardModel):
    """
    Forward Model for two dimensional MT model in the TE mode for a given
    frequency omega.
    It defines a cost function:

      *  defect = 1/2 integrate( sum_s w^s * ( E_x/Hy - Z_XY^s ) ) ** 2  *

    where E_x is the horizontal electric field perpendicular to the YZ-domain,
    horizontal magnetic field H_y=1/(i*omega*mu) * E_{x,z} with complex unit
    i and permeability mu. The weighting factor w^s is set to

        * w^s(X) = w_0^s/(2pi*eta**2) * exp(-length(X-X^s)**2/(2*eta**2))  *

    where X^s is the location of impedance measurement Z_XY^s, w_0 is the level
    of confidence (eg. 1/measurement error) and eta is level of spatial
    confidence.

    E_x is given as solution of the PDE

        * -E_{x,ii} - i omega * mu * sigma * E_x = 0

    where E_x is set to E_0 on the top of the domain. Homogeneous Neuman
    conditions are assumed elsewhere.
    """
    def __init__(self, domain, omega, x, Z_XY, eta=None, w0=1., mu=4*PI*1e-7, Ex_top=1,
        coordinates=None, fixAtTop=False, tol=1e-8, saveMemory=False, directSolver=True):
        """
        initializes a new forward model.

        :param domain: domain of the model
        :type domain: `Domain`
        :param omega: frequency
        :type omega: positive ``float``
        :param x: coordinates of measurements
        :type x: ``list`` of ``tuple`` with ``float``
        :param Z_XY: measured impedance
        :type Z_XY: ``list`` of ``complex``
        :param eta: spatial confidence radius
        :type eta:  positive ``float`` or ``list`` of  positive ``float``
        :param w0: confidence factors for meassurements. If not set one is assumed.
        :type w0: ``None`` or a list of positive ``float``
        :param Ex_top: value of Ex_ at top the domain and if `fixAtTop` is set at the bottom of the domain.
        :type Ex_top: ``float``, ``complex`` or ``Data`` of shape (2,)
        :param mu: permeability
        :type mu: ``float``
        :param coordinates: defines coordinate system to be used (not supported yet)
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param fixAtTop: if true Ex is set Ex_top at the top of the domain
        :type fixAtTop: ``bool``
        :param saveMemory: if true stiffness matrix is deleted after solution
                           of the PDE to minimize memory requests. This will
                           require more compute time as the matrix needs to be
                           reallocated.
        :type saveMemory: ``bool``
        """
        super(MT2DModelTEMode, self).__init__()
        self.__trafo = makeTranformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian coordinates are not supported yet.")

        if len(x) != len(Z_XY):
            raise ValueError("Number of data points and number of impedance values don't match.")

        if eta is None:
            eta = sup(domain.getSize())*0.45

        if isinstance(eta, float) or isinstance(eta, int):
            eta = [float(eta)]*len(Z_XY)
        elif not len(eta) == len(Z_XY):
            raise ValueError("Number of confidence radii and number of impedance values don't match.")

        if isinstance(w0, float) or isinstance(w0, int):
            w0 =[float(w0)]*len(Z_XY)
        elif not len(w0) == len(Z_XY):
            raise ValueError("Number of confidence factors and number of impedance values don't match.")

        self.__domain = domain
        self._omega_mu = omega * mu

        xx=Function(domain).getX()
        f = -1./(complex(0,1)*self._omega_mu)
        scaledZxy = [ z*f for z in Z_XY ]

        totalS=0
        self._scaledZ = Vector([0.,0.], Function(domain))
        self._weight = Scalar(0., Function(domain))

        for s in range(len(scaledZxy)):
            chi = self.getWeightingFactor(xx, 1., x[s], eta[s])
            f = integrate(chi)
            if f < eta[s]**2 * 0.01 :
                raise ValueError("Zero weight (almost) for data point %s. Change eta or refine mesh."%(s,))
            w02 = w0[s]/f
            totalS += w02
            self._scaledZ[0] += chi*scaledZxy[s].real
            self._scaledZ[1] += chi*scaledZxy[s].imag
            self._weight += chi*w02

        if not totalS > 0:
            raise ValueError("Scaling of weight factors failed as sum is zero.")

        DIM = domain.getDim()
        z = domain.getX()[DIM-1]
        self.__ztop = sup(z)
        self.__zbottom = inf(z)
        #====================================
        if fixAtTop:
            if isinstance(Ex_top, float) or isinstance(Ex_top, int) :
                self.__Ex_0 = Data((Ex_top,0), Solution(domain))
            elif isinstance(Ex_top, tuple):
                self.__Ex_0 = Data((Ex_top[0],Ex_top[1]), Solution(domain))
            elif isinstance(Ex_top, complex):
                self.__Ex_0 = Data((Ex_top.real,Ex_top.imag), Solution(domain))
            else:
                if not Ex_top.getShape() == (2,):
                    raise ValueError("Expected shape of Ex_0 is (2,)")
                self.__Ex_0 = Ex_top
            self.__Ex_0 = self.__Ex_0 * whereZero(z-self.__ztop)
        else:
            self.__Ex_0 = [0.,0.]
        #====================================
        self.__tol = tol
        self.__fixAtTop = fixAtTop
        self.__directSolver = directSolver
        self.__saveMemory = saveMemory
        self.__pde = None
        if not saveMemory:
            self.__pde = self.setUpPDE()

    def getDomain(self):
        """
        Returns the domain of the forward model.

        :rtype: `Domain`
        """
        return self.__domain

    def getCoordinateTransformation(self):
        """
        returns the coordinate transformation being used

        :rtype: ``CoordinateTransformation``
        """
        return self.__trafo

    def setUpPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        if self.__pde is None:
            DIM=self.__domain.getDim()
            pde=LinearPDE(self.__domain, numEquations=2)
            if self.__directSolver == True:
                pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
            D=pde.createCoefficient('D')
            A=pde.createCoefficient('A')
            Y=pde.createCoefficient('Y')
            X=pde.createCoefficient('X')
            y=pde.createCoefficient('y')
            A[0,:,0,:]=kronecker(DIM)
            A[1,:,1,:]=kronecker(DIM)
            pde.setValue(A=A, D=D, X=X, Y=Y, y=y)
            z = self.__domain.getX()[DIM-1]
            q = whereZero(z-self.__zbottom)
            if self.__fixAtTop:
                q+=whereZero(z-self.__ztop)
            pde.setValue(q=q*[1,1])
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOff()
        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        return pde

    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: conductivity
        :type sigma: ``Data`` of shape (2,)
        :return: Ex_, Ex_,z
        :rtype: ``Data`` of shape (2,)
        """
        DIM = self.__domain.getDim()
        pde = self.setUpPDE()
        D = pde.getCoefficient('D')
        f = self._omega_mu * sigma
        D[0,1] = -f
        D[1,0] =  f
        Z = FunctionOnBoundary(self.__domain).getX()[DIM-1]
        pde.setValue(D=D, y=whereZero(Z-self.__ztop)*[1.,0.], r=self.__Ex_0)
        u = pde.getSolution()
        return u, grad(u)[:,1]

    def getWeightingFactor(self, x, wx0, x0, eta):
        """
        returns the weighting factor
        """
        return wx0 * whereNegative(length(x-x0)-eta)

    def getDefect(self, sigma, Ex, dExdz):
        """
        Returns the defect value.

        :param sigma: a suggestion for conductivity
        :type sigma: ``Data`` of shape ()
        :param Ex_: electric field
        :type Ex_: ``Data`` of shape (2,)
        :param dExdz: vertical derivative of electric field
        :type dExdz: ``Data`` of shape (2,)

        :rtype: ``float``
        """
        x=dExdz.getFunctionSpace().getX()
        Ex=interpolate(Ex, x.getFunctionSpace())
        u0=Ex[0]
        u1=Ex[1]
        u01=dExdz[0]
        u11=dExdz[1]
        scale = self._weight / ( u01**2 + u11**2 )

        Z = self._scaledZ
        A = integrate(scale * ((Z[0]**2+Z[1]**2)*(u01**2+u11**2) \
                + 2*Z[1]*(u0*u11-u01*u1)
                - 2*Z[0]*(u0*u01+u11*u1) + u0**2 + u1**2))

        return A/2

    def getGradient(self, sigma, Ex, dExdz):
        """
        Returns the gradient of the defect with respect to density.

        :param sigma: a suggestion for conductivity
        :type sigma: ``Data`` of shape ()
        :param Ex: electric field
        :type Ex: ``Data`` of shape (2,)
        :param dExdz: vertical derivative of electric field
        :type dExdz: ``Data`` of shape (2,)
        """
        pde=self.setUpPDE()

        x=dExdz.getFunctionSpace().getX()
        Ex=interpolate(Ex, x.getFunctionSpace())
        u0 = Ex[0]
        u1 = Ex[1]
        u01 = dExdz[0]
        u11 = dExdz[1]

        D=pde.getCoefficient('D')
        Y=pde.getCoefficient('Y')
        if Y.isEmpty():
            Y=pde.createCoefficient('Y')
        X=pde.getCoefficient('X')
        if X.isEmpty():
            X=pde.createCoefficient('X')

        f = self._omega_mu * sigma
        D[0,1] =  f
        D[1,0] = -f

        Z = self._scaledZ
        scale = 1./( u01**2 + u11**2 )
        scale2 = scale**2
        scale *= self._weight
        scale2 *= self._weight

        Y[0] = scale * (u0 - u01*Z[0] + u11*Z[1])
        Y[1] = scale * (u1 - u01*Z[1] - u11*Z[0])
        X[0,1] = scale2 * (2*u01*u11*(Z[0]*u1-Z[1]*u0) \
                + (Z[0]*u0+Z[1]*u1)*(u01**2-u11**2)
                - u01*(u0**2 + u1**2))
        X[1,1] = scale2 * (2*u01*u11*(Z[1]*u1+Z[0]*u0) \
                + (Z[1]*u0-Z[0]*u1)*(u01**2-u11**2)
                - u11*(u0**2 + u1**2))

        pde.setValue(D=D, X=X, Y=Y)
        Zstar=pde.getSolution()
        if self.__saveMemory: self.__pde=None
        return (-self._omega_mu)* (Zstar[1]*u0-Zstar[0]*u1)

