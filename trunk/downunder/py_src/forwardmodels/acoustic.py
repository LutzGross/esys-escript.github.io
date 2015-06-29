##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
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

"""Forward model for acoustic wave forms"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['AcousticWaveForm']

from .base import ForwardModel
from esys.downunder.coordinates import makeTranformation
from esys.escript import Data, DiracDeltaFunctions, FunctionOnBoundary
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.util import *
import numpy as np


class AcousticWaveForm(ForwardModel):
    """
    Forward Model for acoustic waveform inversion in the frequency domain.
    It defines a cost function:

    :math: `defect = 1/2 integrate( ( w * ( a * u - data ) ) ** 2 )`

    where w are weighting factors, data are the measured data (as a 2-comp
    vector of real and imaginary part) for real frequency omega, and u is
    the corresponding result produced by the forward model.
    u (as a 2-comp vector) is the solution of the complex Helmholtz equation
    for frequency omega, source F and complex, inverse, squared p-velocity
    sigma:

    :math: `-u_{ii} - omega**2 * sigma * u = F`

    It is assumed that the exact scale of source F is unknown and the scaling
    factor a of F is calculated by minimizing the defect.
    """
    def __init__(self, domain, omega, w, data, F, coordinates=None,
                 fixAtBottom=False, tol=1e-8, saveMemory=True, scaleF=True):
        """
        initializes a new forward model with acoustic wave form inversion.

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: weighting factors
        :type w: ``Scalar``
        :param data: real and imaginary part of data
        :type data: ``Data`` of shape (2,)
        :param F: real and imaginary part of source given at Dirac points,
                  on surface or at volume.
        :type F: ``Data`` of shape (2,)
        :param coordinates: defines coordinate system to be used (not supported yet)
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param saveMemory: if true stiffness matrix is deleted after solution
                           of PDE to minimize memory requests. This will
                           require more compute time as the matrix needs to be
                           reallocated.
        :type saveMemory: ``bool``
        :param scaleF: if true source F is scaled to minimize defect.
        :type scaleF: ``bool``
        :param fixAtBottom: if true pressure is fixed to zero at the bottom of
                            the domain
        :type fixAtBottom: ``bool``
        """
        super(AcousticWaveForm, self).__init__()
        self.__trafo = makeTranformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian Coordinates are not supported yet.")
        if not isinstance(data, Data):
            raise ValueError("data must be an escript.Data object.")
        if not data.getFunctionSpace() == FunctionOnBoundary(domain):
            raise ValueError("data must be defined on boundary")
        if not data.getShape() == (2,):
            raise ValueError("data must have shape (2,) (real and imaginary part).")
        if w is None:
            w = 1.
        if not isinstance(w, Data):
            w = Data(w, FunctionOnBoundary(domain))
        else:
            if not w.getFunctionSpace() == FunctionOnBoundary(domain):
                raise ValueError("Weights must be defined on boundary.")
            if not w.getShape() == ():
                raise ValueError("Weights must be scalar.")

        self.__domain = domain
        self.__omega = omega
        self.__weight = w
        self.__data = data
        self.scaleF = scaleF
        if scaleF:
            A = integrate(self.__weight*length(self.__data)**2)
            if A > 0:
                self.__data*=1./sqrt(A)

        self.__BX = boundingBox(domain)
        self.edge_lengths = np.asarray(boundingBoxEdgeLengths(domain))

        if not isinstance(F, Data):
            F=interpolate(F,  DiracDeltaFunctions(domain))
        if not F.getShape() == (2,):
            raise ValueError("Source must have shape (2,) (real and imaginary part).")

        self.__F=Data()
        self.__f=Data()
        self.__f_dirac=Data()

        if F.getFunctionSpace() == DiracDeltaFunctions(domain):
            self.__f_dirac=F
        elif F.getFunctionSpace() == FunctionOnBoundary(domain):
            self.__f=F
        else:
            self.__F=F
        self.__tol=tol
        self.__fixAtBottom=fixAtBottom
        self.__pde=None
        if not saveMemory:
            self.__pde=self.setUpPDE()

    def getSurvey(self, index=None):
        """
        Returns the pair (data, weight)

        If argument index is ignored.
        """
        return self.__data, self.__weight

    def rescaleWeights(self, scale=1., sigma_scale=1.):
        """
        rescales the weights such that

        :math: integrate( ( w omega**2 * sigma_scale * data * ((1/L_j)**2)**-1) +1 )/(data*omega**2 * ((1/L_j)**2)**-1) * sigma_scale )=scale

        :param scale: scale of data weighting factors
        :type scale: positive ``float``
        :param sigma_scale: scale of 1/vp**2 velocity.
        :type sigma_scale: ``Scalar``
        """
        raise Warning("rescaleWeights is not tested yet.")
        if not scale > 0:
            raise ValueError("Value for scale must be positive.")
        if not sigma_scale*omega**2*d > 0:
            raise ValueError("Rescaling of weights failed due to zero denominator.")
        # copy back original weights before rescaling
        #self.__weight=[1.*ow for ow in self.__origweight]
        L2=1/length(1/self.edge_length)**2
        d=Lsup(length(data))
        A=integrate(self.__weight*(sigma_scale*omega**2*d+1)/(sigma_scale*omega**2*d) )
        if A > 0:
            self.__weight*=1./A
            if self.scaleF:
                self.__data*=sqrt(A)
        else:
            raise ValueError("Rescaling of weights failed.")

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
        Creates and returns the underlying PDE.

        :rtype: `LinearPDE`
        """
        from esys.escript import getEscriptParamInt
        if self.__pde is None:
            if getEscriptParamInt("PASO_DIRECT")==0:
                raise ValueError("Either this build of escript or the current MPI configuration does not support direct solvers.")
            pde=LinearPDE(self.__domain, numEquations=2)
            D=pde.createCoefficient('D')
            A=pde.createCoefficient('A')
            A[0,:,0,:]=kronecker(self.__domain.getDim())
            A[1,:,1,:]=kronecker(self.__domain.getDim())
            pde.setValue(A=A, D=D)
            if self.__fixAtBottom:
                DIM=self.__domain.getDim()
                z = self.__domain.getX()[DIM-1]
                pde.setValue(q=whereZero(z-self.__BX[DIM-1][0])*[1,1])

            pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOff()
        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        return pde

    def getSourceScaling(self, u):
        """
        returns the scaling factor s required to rescale source F to minimize defect ``|s * u- data|^2``

        :param u: value of pressure solution (real and imaginary part)
        :type u: ``Data`` of shape (2,)
        :rtype: `complex`
        """
        uTu = integrate(self.__weight * length(u)**2)
        uTar = integrate(self.__weight * ( u[0]*self.__data[0]+u[1]*self.__data[1]) )
        uTai = integrate(self.__weight * ( u[0]*self.__data[1]-u[1]*self.__data[0]) )
        if uTu > 0:
            return complex(uTar/uTu, uTai/uTu)
        else:
            return complex(1.,0)

    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: a suggestion for complex 1/V**2
        :type sigma: ``Data`` of shape (2,)
        :return: solution,  uTar, uTai, uTu
        :rtype: ``Data`` of shape (2,), 3 x `float`
        """
        pde=self.setUpPDE()
        D=pde.getCoefficient('D')
        D[0,0]=-self.__omega**2 * sigma[0]
        D[0,1]= self.__omega**2 * sigma[1]
        D[1,0]=-self.__omega**2 * sigma[1]
        D[1,1]=-self.__omega**2 * sigma[0]
        pde.setValue(D=D, Y=self.__F, y=self.__f, y_dirac=self.__f_dirac)
        u=pde.getSolution()

        uTar=integrate(self.__weight * ( u[0]*self.__data[0]+u[1]*self.__data[1]) )
        uTai=integrate(self.__weight * ( u[0]*self.__data[1]-u[1]*self.__data[0]) )
        uTu = integrate( self.__weight * length(u)**2 )
        return u, uTar, uTai, uTu

    def getDefect(self, sigma, u, uTar, uTai, uTu):
        """
        Returns the defect value.

        :param sigma: a suggestion for complex 1/V**2
        :type sigma: ``Data`` of shape (2,)
        :param u: a u vector
        :type u: ``Data`` of shape (2,)
        :param uTar: equals `integrate( w  * (data[0]*u[0]+data[1]*u[1]))`
        :type uTar: `float`
        :param uTai: equals `integrate( w  * (data[1]*u[0]-data[0]*u[1]))`
        :type uTa: `float`
        :param uTu: equals `integrate( w  * (u,u))`
        :type uTu: `float`

        :rtype: ``float``
        """
        # assuming integrate(w * length(data)**2) =1
        if self.scaleF and abs(uTu) >0:
           A = 1.-(uTar**2 + uTai**2)/uTu
        else:
           A = integrate(self.__weight*length(self.__data)**2)- 2 * uTar + uTu
        return  A/2

    def getGradient(self, sigma, u, uTar, uTai, uTu):
        """
        Returns the gradient of the defect with respect to density.

        :param sigma: a suggestion for complex 1/V**2
        :type sigma: ``Data`` of shape (2,)
        :param u: a u vector
        :type u: ``Data`` of shape (2,)
        :param uTar: equals `integrate( w  * (data[0]*u[0]+data[1]*u[1]))`
        :type uTar: `float`
        :param uTai: equals `integrate( w  * (data[1]*u[0]-data[0]*u[1]))`
        :type uTa: `float`
        :param uTu: equals `integrate( w  * (u,u))`
        :type uTu: `float`
        """
        pde=self.setUpPDE()

        if self.scaleF and abs(uTu) >0:
            Z=((uTar**2+uTai**2)/uTu**2) *interpolate(u, self.__data.getFunctionSpace())
            Z[0]+= (-uTar/uTu) * self.__data[0]+ (-uTai/uTu) * self.__data[1]
            Z[1]+= (-uTar/uTu) * self.__data[1]+   uTai/uTu  * self.__data[0]

        else:
            Z = u - self.__data
        if Z.getFunctionSpace() == DiracDeltaFunctions(self.getDomain()):
            pde.setValue(y_dirac=self.__weight * Z)
        else:
            pde.setValue(y=self.__weight * Z)
        D=pde.getCoefficient('D')
        D[0,0]=-self.__omega**2 * sigma[0]
        D[0,1]=-self.__omega**2 * sigma[1]
        D[1,0]= self.__omega**2 * sigma[1]
        D[1,1]=-self.__omega**2 * sigma[0]
        pde.setValue(D=D)
        ZTo2=pde.getSolution()*self.__omega**2
        return inner(ZTo2,u)*[1,0]+(ZTo2[1]*u[0]-ZTo2[0]*u[1])*[0,1]

