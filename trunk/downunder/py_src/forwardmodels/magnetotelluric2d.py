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

"""Forward models for 2D MT (TE and TM mode)"""

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['MT2DModelTEMode', 'MT2DModelTMMode']

from .base import ForwardModel
import esys.escript as escript
import esys.escript.linearPDEs as lpde
import math
from esys.downunder import coordinates as coords

class MT2DBase(ForwardModel):
    """
    Base class for 2D MT forward models. See `MT2DModelTEMode` and
    `MT2DModelTMMode` for actual implementations.
    
    """
    def __init__(self, domain, omega, x, Z, eta=None, w0=1., mu=4*math.pi*1e-7, sigma0=0.01,
                 airLayerLevel=None, fixAirLayer=False,
                 coordinates=None, tol=1e-8, saveMemory=False, directSolver=True):
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
        :param sigma0: background conductivity
        :type sigma0: ``float``
        :param airLayerLevel: position of the air layer from to bottom of the domain. If
                              not set the air layer starts at the top of the domain
        :type airLayerLevel: ``float`` or ``None``        
        :param fixAirLayer: fix air layer (TM mode)
        :type fixAirLayer: ``bool``
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
        self.__trafo = coords.makeTransformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian coordinates are not supported yet.")
        if len(x) != len(Z):
            raise ValueError("Number of data points and number of impedance values don't match.")

        if eta is None:
            eta = escript.sup(domain.getSize())*0.45

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
        self._ks=escript.sqrt(self._omega_mu * sigma0 /2.)

        xx=escript.Function(domain).getX()
        totalS=0
        self._Z = [ escript.Scalar(0., escript.Function(domain)),  escript.Scalar(0., escript.Function(domain)) ]
        self._weight = escript.Scalar(0., escript.Function(domain))

        for s in range(len(Z)):
            chi = self.getWeightingFactor(xx, 1., x[s], eta[s])
            f = escript.integrate(chi)
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
        self._ztop= escript.sup(z)
        self._zbottom = escript.inf(z)
        if airLayerLevel is None:
            airLayerLevel=self._ztop
        self._airLayerLevel=airLayerLevel
                
        # botton:
        mask0=escript.whereZero(z-self._zbottom)
        r=mask0* [ escript.exp(self._ks*(self._zbottom -airLayerLevel))*escript.cos(self._ks*(self._zbottom -airLayerLevel)), 
                   escript.exp(self._ks*(self._zbottom -airLayerLevel))*escript.sin(self._ks*(self._zbottom -airLayerLevel))]

        #top:
        if  fixAirLayer:
          mask1=escript.whereNonNegative(z-airLayerLevel)
          r+=mask1*[ 1, 0 ]
        else:
          mask1=escript.whereZero(z-self._ztop)
          r+=mask1*[ self._ks*(self._ztop-airLayerLevel)+1, self._ks*(self._ztop-airLayerLevel) ]
          
        self._q=(mask0+mask1)*[1,1]
        self._r=r       
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
            pde=lpde.LinearPDE(self.__domain, numEquations=2)
            if self._directSolver == True:
                pde.getSolverOptions().setSolverMethod(lpde.SolverOptions.DIRECT)
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
            return wx0 * escript.whereNegative(escript.length(x-midpoint)-eta)
        except:
            return wx0 * escript.whereNegative(escript.length(x-x0)-eta)

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

    where E_x at top and bottom is set to solution for background field. Homogeneous Neuman
    conditions are assumed elsewhere. 
    """
    def __init__(self, domain, omega, x, Z_XY, eta=None, w0=1., mu=4*math.pi*1e-7, sigma0=0.01,
                 airLayerLevel=None, coordinates=None, Ex_top=1, fixAtTop=False, 
                 tol=1e-8, saveMemory=False, directSolver=True):
        """
        initializes a new forward model. See base class for a description of
        the arguments.
        """
        
        f = -1./(complex(0,1)*omega*mu)
        scaledZXY = [ z*f for z in Z_XY ]
        super(MT2DModelTEMode, self).__init__(domain=domain, omega=omega, x=x, 
                                              Z=scaledZXY, eta=eta, w0=w0, mu=mu, sigma0=sigma0,
                                              airLayerLevel= airLayerLevel,  fixAirLayer=False,
                                              coordinates=coordinates, tol=tol, saveMemory=saveMemory, directSolver=directSolver)

        
    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: conductivity
        :type sigma: ``Data`` of shape (2,)
        :return: E_x, E_z
        :rtype: ``Data`` of shape (2,)
        """
        DIM = self.getDomain().getDim()
        pde = self.setUpPDE()
        D = pde.getCoefficient('D')
        f = self._omega_mu * sigma
        D[0,1] = -f
        D[1,0] =  f
        A= pde.getCoefficient('A')
        A[0,:,0,:]=escript.kronecker(DIM)
        A[1,:,1,:]=escript.kronecker(DIM)
    
        pde.setValue(A=A, D=D,  r=self._r)
        u = pde.getSolution()
        return u, escript.grad(u)[:,1]

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
        Ex=escript.interpolate(Ex, x.getFunctionSpace())
        u0=Ex[0]
        u1=Ex[1]
        u01=dExdz[0]
        u11=dExdz[1]
        scale = self._weight / ( u01**2 + u11**2 )

        Z = self._Z
        A = escript.integrate(scale * ( (Z[0]**2+Z[1]**2)*(u01**2+u11**2)
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
        Ex=escript.interpolate(Ex, x.getFunctionSpace())
        u0 = Ex[0]
        u1 = Ex[1]
        u01 = dExdz[0]
        u11 = dExdz[1]

        D=pde.getCoefficient('D')
        Y=pde.getCoefficient('Y')
        X=pde.getCoefficient('X')
        A= pde.getCoefficient('A')

        A[0,:,0,:]=escript.kronecker(DIM)
        A[1,:,1,:]=escript.kronecker(DIM)
        
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

    where H_x at top and bottom is set to solution for background field. 
    Homogeneous Neuman conditions are assumed elsewhere.
    """
    def __init__(self, domain, omega, x, Z_YX, eta=None, w0=1., mu=4*math.pi*1e-7, sigma0=0.01,
                 airLayerLevel=None, coordinates=None, tol=1e-8, saveMemory=False,
                 directSolver=True):
        """
        initializes a new forward model. See base class for a description of
        the arguments.
        """
        super(MT2DModelTMMode, self).__init__(domain=domain, omega=omega, x=x,
                Z=Z_YX, eta=eta, w0=w0, mu=mu, sigma0=sigma0, 
                airLayerLevel=airLayerLevel, fixAirLayer=True,
                coordinates=coordinates, tol=tol,saveMemory=saveMemory, directSolver=directSolver)

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
        for i in range(DIM):
            A[0,i,0,i]=rho
            A[1,i,1,i]=rho
        
        pde.setValue(A=A, D=D, r=self._r)
        u = pde.getSolution()
        return u, escript.grad(u)

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
        Hx = escript.interpolate(Hx, x.getFunctionSpace())
        u0 = Hx[0]
        u1 = Hx[1]
        u01 = g_Hx[0,1]
        u11 = g_Hx[1,1]
        scale = rho / ( u0**2 + u1**2 )

        Z = self._Z
        A = escript.integrate(self._weight * ( Z[0]**2 + Z[1]**2
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
        Hx=escript.interpolate(Hx, x.getFunctionSpace())
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

        for i in range(DIM):
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
        g=escript.grad(pde.getSolution())

        Hstarr_x = g[0,0]
        Hstari_x = g[1,0]
        Hstarr_z = g[0,1]
        Hstari_z = g[1,1]
        return -scale*(u0*(Z[0]*u01+Z[1]*u11)+u1*(Z[0]*u11-Z[1]*u01)-rho*gscale)\
               - Hstarr_x*u00 - Hstarr_z*u01 - Hstari_x*u10 - Hstari_z*u11

