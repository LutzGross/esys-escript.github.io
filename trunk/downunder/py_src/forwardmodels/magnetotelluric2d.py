from __future__ import print_function
from __future__ import division
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

"""Forward models for 2D MT (TE and TM mode)"""

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['MT2DModelTEMode', 'MT2DModelTMMode']

from .base import ForwardModel
from esys.downunder.coordinates import makeTranformation
from esys.escript import Data, Scalar, Vector, Function, FunctionOnBoundary, Solution
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.util import *
from math import pi as PI


class MT2DBase(ForwardModel):
    """
    Base class for 2D MT forward models. See `MT2DModelTEMode` and
    `MT2DModelTMMode` for actual implementations.
    """
    def __init__(self, domain, omega, x, Z, eta=None, w0=1., mu=4*PI*1e-7,
                 Ftop=1., fixAtTop=False, fixAboveLevel=None, Fbottom=0.,
                 fixAtBottom=False, coordinates=None, tol=1e-8,
                 saveMemory=False, directSolver=True):
        """
        initializes a new forward model.

        :param domain: domain of the model
        :type domain: `Domain`
        :param omega: frequency
        :type omega: positive ``float``
        :param x: coordinates of measurements
        :type x: ``list`` of ``tuple`` with ``float``
        :param Z: measured impedance (possibly scaled)
        :type Z: ``list`` of ``complex``
        :param eta: spatial confidence radius
        :type eta:  positive ``float`` or ``list`` of  positive ``float``
        :param w0: confidence factors for meassurements.
        :type w0: ``None`` or a list of positive ``float``
        :param mu: permeability
        :type mu: ``float``
        :param Ftop: value of field at top of the domain, see `fixAtTop` and
                     `fixAboveLevel`
        :type Ftop: ``float``, ``complex`` or ``Data`` of shape (2,)
        :param fixAtTop: if true F is set to Ftop at the top of the domain.
                         If both `fixAtTop` and `fixAboveLevel` are set, then
                         `fixAboveLevel` takes precedence.
        :type fixAtTop: ``bool``
        :param fixAboveLevel: level above which F is set to Ftop (typically
                              the level of the air layer).
                              Use `fixAtTop` *or* `fixAboveLevel`, not both.
        :type fixAboveLevel : ``float`` or ``None``
        :param Fbottom: value of field at base of the domain
        :type Fbottom: ``float``, ``complex`` or ``Data`` of shape (2,)
        :param fixAtBottom: if true F is set to Fbottom at the bottom of the domain
        :type fixAtBottom: ``bool``
        :param coordinates: defines coordinate system to be used (not supported yet)
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param saveMemory: if true stiffness matrix is deleted after solution
                           of the PDE to minimize memory use. This will
                           require more compute time as the matrix needs to be
                           reallocated at each iteration.
        :type saveMemory: ``bool``
        :param directSolver: if true a direct solver (rather than an iterative
                             solver) will be used to solve the PDE
        :type directSolver: ``bool``
        """
        super(MT2DBase, self).__init__()
        self.__trafo = makeTranformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian coordinates are not supported yet.")
        if len(x) != len(Z):
            raise ValueError("Number of data points and number of impedance values don't match.")

        if eta is None:
            eta = sup(domain.getSize())*0.45

        if isinstance(eta, float) or isinstance(eta, int):
            eta = [float(eta)]*len(Z)
        elif not len(eta) == len(Z):
            raise ValueError("Number of confidence radii and number of impedance values don't match.")

        if isinstance(w0, float) or isinstance(w0, int):
            w0 =[float(w0)]*len(Z)
        elif not len(w0) == len(Z):
            raise ValueError("Number of confidence factors and number of impedance values don't match.")

        self.__domain = domain
        self._omega_mu = omega * mu

        xx=Function(domain).getX()
        totalS=0
        self._Z = [ Scalar(0., Function(domain)),  Scalar(0., Function(domain)) ]
        self._weight = Scalar(0., Function(domain))

        for s in range(len(Z)):
            chi = self.getWeightingFactor(xx, 1., x[s], eta[s])
            f = integrate(chi)
            if f < eta[s]**2 * 0.01 :
                raise ValueError("Zero weight (almost) for data point %s. Change eta or refine mesh."%(s,))
            w02 = w0[s]/f
            totalS += w02
            self._Z[0] += chi*Z[s].real
            self._Z[1] += chi*Z[s].imag
            self._weight += chi*w02/(abs(Z[s])**2)

        if not totalS > 0:
            raise ValueError("Scaling of weight factors failed as sum is zero.")

        DIM = domain.getDim()
        z = domain.getX()[DIM-1]
        self._ztop = sup(z)
        self._zbottom = inf(z)
        self._q=Vector(0.,Solution(domain))
        self._r=Vector(0.,Solution(domain))
        #====================================
        if fixAtTop or fixAboveLevel is not None:
            if fixAboveLevel is not None:
                m=whereNonNegative(z-fixAboveLevel)
            else:
                m=whereZero(z-self._ztop)
            if isinstance(Ftop, float) or isinstance(Ftop, int):
                d = Data((Ftop,0), Solution(domain))
            elif isinstance(Ftop, tuple):
                d = Data((Ftop[0],Ftop[1]), Solution(domain))
            elif isinstance(Ftop, complex):
                d = Data((Ftop.real,Ftop.imag), Solution(domain))
            else:
                if not Ftop.getShape() == (2,):
                    raise ValueError("Expected shape of top value is (2,)")
                d = Ftop
            self._r+=m*d
            self._q+=m*[1.,1.]
        if fixAtBottom:
            m=whereZero(z-self._zbottom)
            if isinstance(Fbottom, float) or isinstance(Fbottom, int) :
                d = Data((Fbottom,0), Solution(domain))
            elif isinstance(Fbottom, tuple):
                d = Data((Fbottom[0],Fbottom[1]), Solution(domain))
            elif isinstance(Fbottom, complex):
                d = Data((Fbottom.real,Fbottom.imag), Solution(domain))
            else:
                if not Fbottom.getShape() == (2,):
                    raise ValueError("Expected shape of top value is (2,)")
                d = Fbottom
            self._r+=m*d
            self._q+=m*[1.,1.]
        #====================================
        self.__tol = tol
        self._directSolver = directSolver
        self._saveMemory = saveMemory
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
            if self._directSolver == True:
                pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
            D=pde.createCoefficient('D')
            A=pde.createCoefficient('A')
            pde.setValue(A=A, D=D, q=self._q)
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOff()
        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        pde.setValue(X=pde.createCoefficient('X'), Y=pde.createCoefficient('Y'))
        return pde

    def getWeightingFactor(self, x, wx0, x0, eta):
        """
        returns the weighting factor
        """
        try:
            origin, spacing, NE = x.getDomain().getGridParameters()
            cell = [int((x0[i]-origin[i])/spacing[i]) for i in range(2)]
            midpoint = [origin[i]+cell[i]*spacing[i]+spacing[i]/2. for i in range(2)]
            return wx0 * whereNegative(length(x-midpoint)-eta)
        except:
            return wx0 * whereNegative(length(x-x0)-eta)

    def getArguments(self, x):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.
        Needs to be implemented in subclasses.
        """
        raise NotImplementedError

    def getDefect(self, x, Ex, dExdz):
        """
        Returns the defect value. Needs to be implemented in subclasses.
        """
        raise NotImplementedError

    def getGradient(self, x, Ex, dExdz):
        """
        Returns the gradient. Needs to be implemented in subclasses.
        """
        raise NotImplementedError


class MT2DModelTEMode(MT2DBase):
    """
    Forward Model for two dimensional MT model in the TE mode for a given
    frequency omega.
    It defines a cost function:

      *  defect = 1/2 integrate( sum_s w^s * ( E_x/H_y - Z_XY^s ) ) ** 2  *

    where E_x is the horizontal electric field perpendicular to the YZ-domain,
    horizontal magnetic field H_y=1/(i*omega*mu) * E_{x,z} with complex unit
    i and permeability mu. The weighting factor w^s is set to 
        *  w^s(X) = w_0^s  *
    if length(X-X^s) <= eta and zero otherwise. X^s is the location of 
    impedance measurement Z_XY^s, w_0^s is the level
    of confidence (eg. 1/measurement error) and eta is level of spatial
    confidence.

    E_x is given as solution of the PDE

        * -E_{x,ii} - i omega * mu * sigma * E_x = 0

    where by default the normal derivative of E_x at the top of the domain 
    is set to Ex_n=1 and E_x is set to zero at the bottom. Homogeneous Neuman
    conditions are assumed elsewhere. If fixAtTop is set E_x is set to Ex_top. 
    """
    def __init__(self, domain, omega, x, Z_XY, eta=None, w0=1., mu=4*PI*1e-7,
                 Ex_n=1, coordinates=None, Ex_top=1, fixAtTop=False, tol=1e-8,
                 saveMemory=False, directSolver=True):
        """
        initializes a new forward model. See base class for a description of
        the arguments.
        """
        
        f = -1./(complex(0,1)*omega*mu)
        scaledZXY = [ z*f for z in Z_XY ]
        self.__Ex_n=complex(Ex_n)
        super(MT2DModelTEMode, self).__init__(domain=domain, omega=omega, x=x, 
                                              Z=scaledZXY, eta=eta, w0=w0, mu=mu, 
                                              Ftop=Ex_top, fixAtTop=fixAtTop, fixAboveLevel= None, Fbottom=0., fixAtBottom=True, 
                                              coordinates=coordinates, tol=tol, saveMemory=saveMemory, directSolver=directSolver)

        
    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: conductivity
        :type sigma: ``Data`` of shape (2,)
        :return: Ex_, Ex_,z
        :rtype: ``Data`` of shape (2,)
        """
        DIM = self.getDomain().getDim()
        pde = self.setUpPDE()
        D = pde.getCoefficient('D')
        f = self._omega_mu * sigma
        D[0,1] = -f
        D[1,0] =  f
        A= pde.getCoefficient('A')
        A[0,:,0,:]=kronecker(DIM)
        A[1,:,1,:]=kronecker(DIM)
        
        Z = FunctionOnBoundary(self.getDomain()).getX()[DIM-1]
        pde.setValue(A=A, D=D, y=whereZero(Z-self._ztop)*[self.__Ex_n.real,self.__Ex_n.imag],  r=self._r)
        u = pde.getSolution()
        return u, grad(u)[:,1]

    def getDefect(self, sigma, Ex, dExdz):
        """
        Returns the defect value.

        :param sigma: a suggestion for conductivity
        :type sigma: ``Data`` of shape ()
        :param Ex: electric field
        :type Ex: ``Data`` of shape (2,)
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

        Z = self._Z
        A = integrate(scale * ( (Z[0]**2+Z[1]**2)*(u01**2+u11**2)
                              + 2*Z[1]*(u0*u11-u01*u1)
                              - 2*Z[0]*(u0*u01+u11*u1)
                              + u0**2 + u1**2 ))

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
        DIM = self.getDomain().getDim()

        x=dExdz.getFunctionSpace().getX()
        Ex=interpolate(Ex, x.getFunctionSpace())
        u0 = Ex[0]
        u1 = Ex[1]
        u01 = dExdz[0]
        u11 = dExdz[1]

        D=pde.getCoefficient('D')
        Y=pde.getCoefficient('Y')
        X=pde.getCoefficient('X')
        A= pde.getCoefficient('A')

        A[0,:,0,:]=kronecker(DIM)
        A[1,:,1,:]=kronecker(DIM)
        
        f = self._omega_mu * sigma
        D[0,1] =  f
        D[1,0] = -f

        Z = self._Z
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

        pde.setValue(A=A, D=D, X=X, Y=Y)
        Zstar=pde.getSolution()
        return (-self._omega_mu)* (Zstar[1]*u0-Zstar[0]*u1)


class MT2DModelTMMode(MT2DBase):
    """
    Forward Model for two-dimensional MT model in the TM mode for a given
    frequency omega.
    It defines a cost function:

      *  defect = 1/2 integrate( sum_s w^s * ( rho*H_x/Hy - Z_YX^s ) ) ** 2  *

    where H_x is the horizontal magnetic field perpendicular to the YZ-domain,
    horizontal magnetic field H_y=1/(i*omega*mu) * E_{x,z} with complex unit
    i and permeability mu. The weighting factor w^s is set to 
        *  w^s(X) = w_0^s  *
    if length(X-X^s) <= eta and zero otherwise. X^s is the location of 
    impedance measurement Z_XY^s, w_0^s is the level
    of confidence (eg. 1/measurement error) and eta is level of spatial
    confidence.

    H_x is given as solution of the PDE

        * -(rho*H_{x,i})_{,i} + i omega * mu * H_x = 0

    where H_x is set to Hx_top on the top of the domain. Homogeneous Neuman
    conditions are assumed elsewhere. If fixAtBottom is set H_x is set to Hx_bottom
    at the bottom of the domain overwrtining the Neuman condition. 
    If fixAboveLevel is set then H_x is set to Hx_top for any location above and including fixAboveLevel
    typically including teh top boundary. 
    """
    def __init__(self, domain, omega, x, Z_YX, eta=None, w0=1., mu=4*PI*1e-7,
                 fixAboveLevel=None, Hx_top=1, coordinates=None, Hx_bottom=1.,
                 fixAtBottom=False, tol=1e-8, saveMemory=False,
                 directSolver=True):
        """
        initializes a new forward model. See base class for a description of
        the arguments.
        """
        super(MT2DModelTMMode, self).__init__(domain=domain, omega=omega, x=x,
                Z=Z_YX, eta=eta, w0=w0, mu=mu, Ftop=Hx_top, fixAtTop=True,
                fixAboveLevel=fixAboveLevel, Fbottom=Hx_bottom,
                fixAtBottom=fixAtBottom, coordinates=coordinates, tol=tol,
                saveMemory=saveMemory, directSolver=directSolver)

    def getArguments(self, rho):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param rho: resistivity
        :type rho: ``Data`` of shape (2,)
        :return: Hx, grad(Hx)
        :rtype: ``tuple`` of ``Data``
        """
        DIM = self.getDomain().getDim()
        pde = self.setUpPDE()
        
        D = pde.getCoefficient('D')
        f = self._omega_mu
        D[0,1] = -f
        D[1,0] =  f
        
        A= pde.getCoefficient('A')
        for i in xrange(DIM):
            A[0,i,0,i]=rho
            A[1,i,1,i]=rho
        
        pde.setValue(A=A, D=D, r=self._r)
        u = pde.getSolution()
        return u, grad(u)

    def getDefect(self, rho, Hx, g_Hx):
        """
        Returns the defect value.

        :param rho: a suggestion for resistivity
        :type rho: ``Data`` of shape ()
        :param Hx: magnetic field
        :type Hx: ``Data`` of shape (2,)
        :param g_Hx: gradient of magnetic field
        :type g_Hx: ``Data`` of shape (2,2)

        :rtype: ``float``
        """
        x = g_Hx.getFunctionSpace().getX()
        Hx = interpolate(Hx, x.getFunctionSpace())
        u0 = Hx[0]
        u1 = Hx[1]
        u01 = g_Hx[0,1]
        u11 = g_Hx[1,1]
        scale = rho / ( u0**2 + u1**2 )

        Z = self._Z
        A = integrate(self._weight * ( Z[0]**2 + Z[1]**2
                                    + scale*(-2*Z[0]*(u0*u01 + u1*u11)
                                             +2*Z[1]*(u1*u01 - u0*u11)
                                             +rho*(u01**2 + u11**2)) ))
        return A/2

    def getGradient(self, rho, Hx, g_Hx):
        """
        Returns the gradient of the defect with respect to resistivity.

        :param rho: a suggestion for resistivity
        :type rho: ``Data`` of shape ()
        :param Hx: magnetic field
        :type Hx: ``Data`` of shape (2,)
        :param g_Hx: gradient of magnetic field
        :type g_Hx: ``Data`` of shape (2,2)
        """
        pde=self.setUpPDE()
        DIM = self.getDomain().getDim()

        x=g_Hx.getFunctionSpace().getX()
        Hx=interpolate(Hx, x.getFunctionSpace())
        u0 = Hx[0]
        u1 = Hx[1]
        u00 = g_Hx[0,0]
        u10 = g_Hx[1,0]
        u01 = g_Hx[0,1]
        u11 = g_Hx[1,1]

        A=pde.getCoefficient('A')
        D=pde.getCoefficient('D')
        Y=pde.getCoefficient('Y')
        X=pde.getCoefficient('X')

        for i in xrange(DIM):
            A[0,i,0,i]=rho
            A[1,i,1,i]=rho

        f = self._omega_mu
        D[0,1] =  f
        D[1,0] = -f

        Z = self._Z
        scale = 1./( u0**2 + u1**2 )
        scale2 = scale**2
        scale *= self._weight
        scale2 *= rho*self._weight
        rho_scale = rho*scale

        gscale = u01**2 + u11**2

        Y[0] = scale2 * ( (Z[0]*u01+Z[1]*u11)*(u0**2-u1**2)
                      + 2*u0*u1*(Z[0]*u11-Z[1]*u01)
                      - rho*u0*gscale )
        Y[1] = scale2 * ( (Z[0]*u11-Z[1]*u01)*(u1**2-u0**2)
                      + 2*u0*u1*(Z[0]*u01+Z[1]*u11)
                      - rho*u1*gscale )
        X[0,1] = rho_scale * (-Z[0]*u0 + Z[1]*u1 + rho*u01)
        X[1,1] = rho_scale * (-Z[0]*u1 - Z[1]*u0 + rho*u11)

        pde.setValue(A=A, D=D, X=X, Y=Y)
        g=grad(pde.getSolution())

        Hstarr_x = g[0,0]
        Hstari_x = g[1,0]
        Hstarr_z = g[0,1]
        Hstari_z = g[1,1]
        return -scale*(u0*(Z[0]*u01+Z[1]*u11)+u1*(Z[0]*u11-Z[1]*u01)-rho*gscale)\
               - Hstarr_x*u00 - Hstarr_z*u01 - Hstari_x*u10 - Hstari_z*u11

