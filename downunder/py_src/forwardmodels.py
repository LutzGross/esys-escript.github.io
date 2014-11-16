from __future__ import print_function
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

"""Collection of forward models that define the inversion problem"""

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['ForwardModel','ForwardModelWithPotential','GravityModel','MagneticModel', 'SelfDemagnetizationModel', 'AcousticWaveForm', 'MT2DModelTEMode','DcRes', 'IsostaticPressure', 'Subsidence']


from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
from esys.escript import Data, Vector, Scalar, Function, DiracDeltaFunctions, FunctionOnBoundary, Solution, length, exp
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDESystem, LinearPDE, SolverOptions
from .coordinates import makeTranformation
from esys.escript.util import *
from math import pi as PI
from esys.weipa import saveSilo
import numpy as np


class ForwardModel(object):
    """
    An abstract forward model that can be plugged into a cost function.
    Subclasses need to implement `getDefect()`, `getGradient()`, and possibly
    `getArguments()` and 'getCoordinateTransformation'.
    """
    def __init__(self):
        pass

    def getArguments(self, x):
        return ()

    def getDefect(self, x, *args):
        raise NotImplementedError

    def getGradient(self, x, *args):
        raise NotImplementedError

    def getCoordinateTransformation(self):
        return None


class ForwardModelWithPotential(ForwardModel):
    """
    Base class for a forward model using a potential such as magnetic or
    gravity. It defines a cost function:

        defect = 1/2 sum_s integrate( ( weight_i[s] * ( r_i - data_i[s] ) ) ** 2 )

    where s runs over the survey, weight_i are weighting factors, data_i are
    the data, and r_i are the results produced by the forward model.
    It is assumed that the forward model is produced through postprocessing
    of the solution of a potential PDE.
    """
    def __init__(self, domain, w, data,  coordinates=None,
                                 fixPotentialAtBottom=False,
                                 tol=1e-8):
        """
        initializes a new forward model with potential.

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Vector`` or list of ``Vector``
        :param data: data
        :type data: ``Vector`` or list of ``Vector``
        :param coordinates: defines coordinate system to be used
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param fixPotentialAtBottom: if true potential is fixed to zero at the bottom of the domain
                                     in addition to the top.
        :type fixPotentialAtBottom: ``bool``
        """
        super(ForwardModelWithPotential, self).__init__()
        self.__domain = domain
        self.__trafo=makeTranformation(domain, coordinates)

        try:
            n=len(w)
            m=len(data)
            if not m == n:
                raise ValueError("Length of weight and data must be the same.")
            self.__weight = w
            self.__data = data
        except TypeError:
            self.__weight = [w]
            self.__data = [data]

        BX = boundingBox(domain)
        DIM = domain.getDim()
        x = domain.getX()
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()
        z=x[DIM-1]
        q0=whereZero(z-BX[DIM-1][1])
        if fixPotentialAtBottom: q0+=whereZero(z-BX[DIM-1][0])
        self.__pde.setValue(q=q0)

        self.edge_lengths=np.asarray(boundingBoxEdgeLengths(domain))
        self.diameter=1./sqrt(sum(1./self.edge_lengths**2))

        self.__origweight=[]
        for s in range(len(self.__weight)):
            # save a copy of the original weights in case of rescaling
            self.__origweight.append(1.*self.__weight[s])

        if not self.__trafo.isCartesian():
            fd=1./self.__trafo.getScalingFactors()
            fw=self.__trafo.getScalingFactors()*sqrt(self.__trafo.getVolumeFactor())
            for s in range(len(self.__weight)):
                self.__weight[s] = fw * self.__weight[s]
                self.__data[s]   = fd * self.__data[s]

    def _rescaleWeights(self, scale=1., fetch_factor=1.):
        """
        rescales the weights such that

        *sum_s integrate( ( weight_i[s] *data_i[s]) (weight_j[s]*1/L_j) * L**2 * fetch_factor )=scale*
        """
        if not scale > 0:
             raise ValueError("Value for scale must be positive.")
        A=0
        # copy back original weights before rescaling
        self.__weight=[1.*ow for ow in self.__origweight]

        for s in range(len(self.__weight)):
            A += integrate(abs(inner(self.__weight[s], self.__data[s]) * inner(self.__weight[s], 1/self.edge_lengths) * fetch_factor))
        if A > 0:
            A=sqrt(scale/A)/self.diameter
            if not self.__trafo.isCartesian():
                A*=self.__trafo.getScalingFactors()*sqrt(self.__trafo.getVolumeFactor())
            for s in range(len(self.__weight)):
                self.__weight[s]*=A
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

    def getPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        return self.__pde

    def _getDefect(self, result):
        """
        Returns the defect value.

        :param result: a result vector
        :type result: `Vector`
        :rtype: ``float``
        """
        A=0.
        for s in range(len(self.__weight)):
            A += integrate( inner(self.__weight[s], self.__data[s]-result)**2 )
        return  A/2

    def getDefectGradient(self, result):
        Y=0.
        for s in range(len(self.__weight)):
            Y = inner(self.__weight[s], self.__data[s]-result) * self.__weight[s] + Y
        return Y

    def getSurvey(self, index=None):
        """
        Returns the pair (data_index, weight_index), where data_i is the data
        of survey i, weight_i is the weighting factor for survey i.
        If index is None, all surveys will be returned in a pair of lists.
        """
        if index is None:
            return self.__data, self.__weight
        if index>=len(self.__data):
            raise IndexError("Forward model only has %d surveys"%len(self.__data))
        return self.__data[index], self.__weight[index]



class GravityModel(ForwardModelWithPotential):
    """
    Forward Model for gravity inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, w, g,  gravity_constant=U.Gravitational_Constant,
                  coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
        """
        Creates a new gravity model on the given domain with one or more
        surveys (w, g).

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Vector`` or list of ``Vector``
        :param g: gravity anomaly data
        :type g: ``Vector`` or list of ``Vector``
        :param coordinates: defines coordinate system to be used
        :type coordinates: ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param fixPotentialAtBottom: if true potential is fixed to zero at the bottom of the domain
                                     in addition to the top.
        :type fixPotentialAtBottom: ``bool``


        :note: It is advisable to call rescaleWeights() to rescale weights before start inversion.
        """
        super(GravityModel, self).__init__(domain, w, g, coordinates, fixPotentialAtBottom, tol)

        trafo = self.getCoordinateTransformation()
        if not trafo.isCartesian():
            self.__G = 4*PI*gravity_constant * trafo.getVolumeFactor()
                    #* trafo.getReferenceSystem().getHeightUnit()**(-3)

            fw = trafo.getScalingFactors()**2 * trafo.getVolumeFactor()
            A=self.getPDE().createCoefficient("A")
            for i in range(self.getDomain().getDim()): A[i,i]=fw[i]
            self.getPDE().setValue(A=A)
        else: # cartesian
            self.__G = 4*PI*gravity_constant
            self.getPDE().setValue(A=kronecker(self.getDomain()))

    def rescaleWeights(self, scale=1., rho_scale=1.):
        """
        rescales the weights such that

        *sum_s integrate( ( w_i[s] *g_i[s]) (w_j[s]*1/L_j) * L**2 * 4*pi*G*rho_scale )=scale*

        :param scale: scale of data weighting factors
        :type scale: positive ``float``
        :param rho_scale: scale of density.
        :type rho_scale: ``Scalar``
        """
        self._rescaleWeights(scale, self.__G*rho_scale)

    def getArguments(self, rho):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param rho: a suggestion for the density distribution
        :type rho: ``Scalar``
        :return: gravity potential and corresponding gravity field.
        :rtype: ``Scalar``, ``Vector``
        """
        phi = self.getPotential(rho)
        gravity_force = -grad(phi)
        return phi, gravity_force

    def getPotential(self, rho):
        """
        Calculates the gravity potential for a given density distribution.

        :param rho: a suggestion for the density distribution
        :type rho: ``Scalar``
        :return: gravity potential
        :rtype: ``Scalar``
        """
        pde=self.getPDE()

        pde.resetRightHandSideCoefficients()
        pde.setValue(Y=-self.__G*rho)
        phi=pde.getSolution()

        return phi

    def getDefect(self, rho, phi, gravity_force):
        """
        Returns the value of the defect

        :param rho: density distribution
        :type rho: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param gravity_force: gravity force
        :type gravity_force: ``Vector``
        :rtype: ``float``
        """
        return self._getDefect(gravity_force)

    def getGradient(self, rho, phi, gravity_force):
        """
        Returns the gradient of the defect with respect to density.

        :param rho: density distribution
        :type rho: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param gravity_force: gravity force
        :type gravity_force: ``Vector``
        :rtype: ``Scalar``
        """
        pde=self.getPDE()
        pde.resetRightHandSideCoefficients()
        pde.setValue(X=self.getDefectGradient(gravity_force))
        ZT=pde.getSolution()
        return ZT*(-self.__G)

class IsostaticPressure(object):
    """
    class to calculate isostatic pressure field correction due to gravity forces 
    """
    def __init__(self, domain,
                        level0=0, 
                        gravity0=-9.81 * U.m*U.sec**(-3),
                        background_density= 2670* U.kg*U.m**(-3),
                        gravity_constant=U.Gravitational_Constant, 
                        coordinates=None, 
                        tol=1e-8):
        """
        
        :param domain: domain of the model
        :type domain: `Domain`
        :param background_density: defines background_density in kg/m^3
        :type background_density: ``float``
        :param coordinates: defines coordinate system to be used
        :type coordinates: ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param level0: pressure for z>=`level0` is set to zero.
        :type level0: ``float``
        :param gravity0: vertical background gravity at `level0`
        :type gravity0: ``float``        
        """
        DIM=domain.getDim()
        self.__domain = domain
        self.__trafo=makeTranformation(domain, coordinates)
        
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()

        z = domain.getX()[DIM-1]
        self.__pde.setValue(q=whereNonNegative(z-level0))

        fw = self.__trafo.getScalingFactors()**2 * self.__trafo.getVolumeFactor()
        A=self.__pde.createCoefficient("A")
        for i in range(DIM): A[i,i]=fw[i]
        self.__pde.setValue(A=A)
        self.__g_b= 4*PI*gravity_constant/self.__trafo.getScalingFactors()[DIM-1]*background_density*(level0-z) + gravity0
        self.__rho_b=background_density
        
    def getPressure(self, g = None, rho=None):
        """
        return the pressure for gravity force anomaly `g` and 
        density anomaly `density`
        
        :param g: gravity anomaly data
        :type g: ``Vector`` 
        :param rho: gravity anomaly data
        :type rho: ``Scalar`` 
        :return: pressure distribution
        :rtype: ``Scalar`
        """
        if not g: g=Vector(0., Function(self.__domain))
        if not rho: rho=Scalar(0., Function(self.__domain))
        
        g2=(rho * self.__g_b) * [0,0,1] + self.__rho_b  * g + rho * g 
        d=self.__trafo.getScalingFactors()
        V= self.__trafo.getVolumeFactor()
        self.__pde.setValue(X= - g2*d*V )
        return self.__pde.getSolution()

    
class MagneticModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, w, B, background_magnetic_flux_density, coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
        """
        Creates a new magnetic model on the given domain with one or more
        surveys (w, B).

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Vector`` or list of ``Vector``
        :param B: magnetic field data
        :type B: ``Vector`` or list of ``Vector``
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param background_magnetic_flux_density: background magnetic flux
               density (in Tesla) with components (B_east, B_north, B_vertical)
        :type background_magnetic_flux_density: ``Vector`` or list of `float`
        :param coordinates: defines coordinate system to be used
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param fixPotentialAtBottom: if true potential is fixed to zero at the
                                     bottom of the domain in addition to the top
        :type fixPotentialAtBottom: ``bool``
        """
        super(MagneticModel, self).__init__(domain, w, B, coordinates, fixPotentialAtBottom, tol)
        background_magnetic_flux_density =interpolate(background_magnetic_flux_density, B[0].getFunctionSpace())
        if not self.getCoordinateTransformation().isCartesian():
            #self.__F = self.getCoordinateTransformation().getVolumeFactor()
            self.__B_r=background_magnetic_flux_density*self.getCoordinateTransformation().getScalingFactors()*self.getCoordinateTransformation().getVolumeFactor()
            self.__B_b=background_magnetic_flux_density/self.getCoordinateTransformation().getScalingFactors()

            A=self.getPDE().createCoefficient("A")
            fw=self.getCoordinateTransformation().getScalingFactors()**2*self.getCoordinateTransformation().getVolumeFactor()
            for i in range(self.getDomain().getDim()): A[i,i]=fw[i]
            self.getPDE().setValue(A=A)
        else: # cartesian
            self.getPDE().setValue(A=kronecker(self.getDomain()))
            self.__B_r=background_magnetic_flux_density
            self.__B_b=background_magnetic_flux_density

    def rescaleWeights(self, scale=1., k_scale=1.):
        """
        rescales the weights such that

        *sum_s integrate( ( w_i[s] *B_i[s]) (w_j[s]*1/L_j) * L**2 * (background_magnetic_flux_density_j[s]*1/L_j) * k_scale )=scale*

        :param scale: scale of data weighting factors
        :type scale: positive ``float``
        :param k_scale: scale of susceptibility.
        :type k_scale: ``Scalar``
        """
        self._rescaleWeights(scale, inner(self.__B_r,1/self.edge_lengths ) * k_scale)

    def getArguments(self, k):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param k: susceptibility
        :type k: ``Scalar``
        :return: scalar magnetic potential and corresponding magnetic field
        :rtype: ``Scalar``, ``Vector``
        """
        phi = self.getPotential(k)
        magnetic_flux_density = k * self.__B_b -grad(phi)
        return phi, magnetic_flux_density

    def getPotential(self, k):
        """
        Calculates the magnetic potential for a given susceptibility.

        :param k: susceptibility
        :type k: ``Scalar``
        :return: magnetic potential
        :rtype: ``Scalar``
        """
        pde=self.getPDE()

        pde.resetRightHandSideCoefficients()
        pde.setValue(X = k* self.__B_r)
        phi=pde.getSolution()

        return phi

    def getDefect(self, k, phi, magnetic_flux_density):
        """
        Returns the value of the defect.

        :param k: susceptibility
        :type k: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param magnetic_flux_density: magnetic field
        :type magnetic_flux_density: ``Vector``
        :rtype: ``float``
        """
        return self._getDefect(magnetic_flux_density)

    def getGradient(self, k, phi, magnetic_flux_density):
        """
        Returns the gradient of the defect with respect to susceptibility.

        :param k: susceptibility
        :type k: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param magnetic_flux_density: magnetic field
        :type magnetic_flux_density: ``Vector``
        :rtype: ``Scalar``
        """
        Y=self.getDefectGradient(magnetic_flux_density)
        pde=self.getPDE()
        pde.resetRightHandSideCoefficients()
        pde.setValue(X=Y)
        YT=pde.getSolution()
        return inner(grad(YT),self.__B_r) -inner(Y,self.__B_b)

class SelfDemagnetizationModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic inversion with self-demagnetization as
    described in the inversion cookbook.
    """
    def __init__(self, domain, w, B, background_magnetic_flux_density, coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
        """
        Creates a new magnetic model on the given domain with one or more
        surveys (w, B).

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Vector`` or list of ``Vector``
        :param B: magnetic field data
        :type B: ``Vector`` or list of ``Vector``
        :param background_magnetic_flux_density: background magnetic flux
               density (in Tesla) with components (B_east, B_north, B_vertical)
        :type background_magnetic_flux_density: ``Vector`` or list of `float`
        :param coordinates: defines coordinate system to be used
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param fixPotentialAtBottom: if true potential is fixed to zero at the
                                     bottom of the domain in addition to the top
        :type fixPotentialAtBottom: ``bool``
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        """
        super(SelfDemagnetizationModel, self).__init__(domain, w, B, coordinates, fixPotentialAtBottom, tol)
        #=========================================================
        background_magnetic_flux_density = interpolate(background_magnetic_flux_density, B[0].getFunctionSpace())
        if not self.getCoordinateTransformation().isCartesian():
            #self.__F = self.getCoordinateTransformation().getVolumeFactor()
            self.__B_r=background_magnetic_flux_density*self.getCoordinateTransformation().getScalingFactors()*self.getCoordinateTransformation().getVolumeFactor()
            self.__B_b=background_magnetic_flux_density/self.getCoordinateTransformation().getScalingFactors()

            self.__fw=self.getCoordinateTransformation().getScalingFactors()**2*self.getCoordinateTransformation().getVolumeFactor()
        else: # cartesian
            self.__fw=1
            self.__B_r=background_magnetic_flux_density
            self.__B_b=background_magnetic_flux_density

        # keep track of k used to build PDE:
        self.__last_k=None
        # this is just the initial set_up
        A=self.getPDE().createCoefficient("A")
        for i in range(self.getDomain().getDim()): A[i,i]=1.
        self.getPDE().setValue(A=A)

    def rescaleWeights(self, scale=1., k_scale=1.):
        """
        rescales the weights such that

        *sum_s integrate( ( w_i[s] *B_i[s]) (w_j[s]*1/L_j) * L**2 * (background_magnetic_flux_density_j[s]*1/L_j) * k_scale )=scale*

        :param scale: scale of data weighting factors
        :type scale: positive ``float``
        :param k_scale: scale of susceptibility.
        :type k_scale: ``Scalar``
        """
        self._rescaleWeights(scale, inner(self.__B_r,1/self.edge_lengths ) * k_scale)

    def getArguments(self, k):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param k: susceptibility
        :type k: ``Scalar``
        :return: scalar magnetic potential and corresponding magnetic field
        :rtype: ``Scalar``, ``Vector``
        """
        phi = self.getPotential(k)
        grad_phi=grad(phi)
        magnetic_flux_density = k * self.__B_b -(1+k)*grad_phi
        return phi, grad_phi, magnetic_flux_density

    def __updateA(self, k):
        """
        updates PDE coefficient if PDE is used with new k
        """
        pde=self.getPDE()
        if self.__last_k is not k:
           A=pde.getCoefficient('A')
           if self.getCoordinateTransformation().isCartesian():
               for i in range(self.getDomain().getDim()): A[i,i]=(1+k)
           else:
               for i in range(self.getDomain().getDim()): A[i,i]=(1+k)*self.__fw[i]
           self.__last_k = k
           pde.setValue(A=A)

    def getPotential(self, k):
        """
        Calculates the magnetic potential for a given susceptibility.

        :param k: susceptibility
        :type k: ``Scalar``
        :return: magnetic potential
        :rtype: ``Scalar``
        """
        self.__updateA(k)
        pde=self.getPDE()
        pde.resetRightHandSideCoefficients()
        pde.setValue(X = k* self.__B_r)
        phi=pde.getSolution()
        return phi

    def getDefect(self, k, phi, grad_phi, magnetic_flux_density):
        """
        Returns the value of the defect.

        :param k: susceptibility
        :type k: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param magnetic_flux_density: magnetic field
        :type magnetic_flux_density: ``Vector``
        :rtype: ``float``
        """
        return self._getDefect(magnetic_flux_density)

    def getGradient(self, k, phi, grad_phi, magnetic_flux_density):
        """
        Returns the gradient of the defect with respect to susceptibility.

        :param k: susceptibility
        :type k: ``Scalar``
        :param phi: corresponding potential
        :type phi: ``Scalar``
        :param magnetic_flux_density: magnetic field
        :type magnetic_flux_density: ``Vector``
        :rtype: ``Scalar``
        """
        self.__updateA(k)
        Y=self.getDefectGradient(magnetic_flux_density)
        pde=self.getPDE()
        pde.resetRightHandSideCoefficients()
        pde.setValue(X=(1+k)*Y)
        grad_YT=grad(pde.getSolution())

        if self.getCoordinateTransformation().isCartesian(): # then b_r=B_b
            return inner(grad_YT-Y, self.__B_r-grad_phi)
        else:
            return inner(grad_YT,self.__B_r-grad_phi)+inner(Y,grad_phi-self.__B_b)

class AcousticWaveForm(ForwardModel):
    """
    Forward Model for acoustic waveform inversion in the frequency domain.
    It defines a cost function:

        defect = 1/2 integrate( ( w * ( a * u - data ) ) ** 2 )

    where w are weighting factors, data are the measured data (as a 2-comp
    vector of real and imaginary part) for real frequency omega, and u is
    the corresponding result produced by the forward model.
    u (as a 2-comp vector) is the solution of the complex Helmholtz equation
    for frequency omega, source F and complex, inverse, squared p-velocity
    sigma:

        * -u_{ii} - omega**2 * sigma * u = F

    It is assumed that the exact scale of source F is unknown and the scaling
    factor a of F is calculated by minimizing the defect.
    """
    def __init__(self, domain, omega, w, data, F, coordinates=None, fixAtBottom=False, tol=1e-8, saveMemory=True, scaleF=True):
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
        self.__domain = domain
        self.__omega=omega
        self.scaleF=scaleF
        self.__trafo=makeTranformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian Coordinates are not supported yet.")
        if not isinstance(data, Data):
                raise ValueError("data must be an escript.Data object.")
        if not data.getFunctionSpace() == FunctionOnBoundary(domain):
                raise ValueError("data must be defined on boundary")
        if not data.getShape() == (2,):
                raise ValueError("data must have shape 2 (real and imaginary part).")
        if w is None:
            w = 1.
        if not isinstance(w, Data):
            w=Data(w, FunctionOnBoundary(domain))
        else:
            if not w.getFunctionSpace() == FunctionOnBoundary(domain):
                raise ValueError("Weights must be defined on boundary.")
            if not w.getShape() == ():
                raise ValueError("Weights must be scalar.")
        self.__weight = w
        self.__data = data
        if scaleF:
            A = integrate(self.__weight*length(self.__data)**2)
            if A>0:
                self.__data*=1./sqrt(A)

        self.__BX = boundingBox(domain)
        self.edge_lengths=np.asarray(boundingBoxEdgeLengths(domain))

        self.__F=Data()
        self.__f=Data()
        self.__f_dirac=Data()

        if not isinstance(F, Data):
            F=interpolate(F,  DiracDeltaFunctions(domain))
        if not F.getShape() == (2,):
            raise ValueError("Source must have shape (2,) (real and imaginary part).")

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

        *integrate( ( w omega**2 * sigma_scale * data * ((1/L_j)**2)**-1) +1 )/(data*omega**2 * ((1/L_j)**2)**-1) * sigma_scale )=scale*

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
            if self.scaleF: self.__data*=sqrt(A)
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

            if getEscriptParamInt("PASO_DIRECT")==0:
                raise ValueError("Either this build of escript or the current MPI configuration does not support direct solvers.")
            pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOff()
        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        return pde

    def getSourceScaling(self, u):
        """
        returns the scaling factor s required to rescale source F to minimize defect |s * u- data|^2
        :param u: value of pressure solution (real and imaginary part)
        :type u: ``Data`` of shape (2,)
        :rtype: `complex`
        """
        uTu=integrate(self.__weight * length(u)**2)
        uTar=integrate(self.__weight * ( u[0]*self.__data[0]+u[1]*self.__data[1]) )
        uTai=integrate(self.__weight * ( u[0]*self.__data[1]-u[1]*self.__data[0]) )
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
        :param uTar : = integrate( w  * (data[0]*u[0]+data[1]*u[1]))
        :type uTar: float
        :param uTai : = integrate( w  * (data[1]*u[0]-data[0]*u[1]))
        :type uTa: float
        :param uTu : = integrate( w  * (u,u))
        :type uTu: float

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
        :param uTar : = integrate( w  * (data[0]*u[0]+data[1]*u[1]))
        :type uTar: float
        :param uTai : = integrate( w  * (data[1]*u[0]-data[0]*u[1]))
        :type uTa: float
        :param uTu : = integrate( w  * (u,u))
        :type uTu: float
        """
        pde=self.setUpPDE()

        if self.scaleF and abs(uTu) >0:
            Z=((uTar**2+uTai**2)/uTu**2) *interpolate(u, self.__data.getFunctionSpace())
            Z[0]+= (-uTar/uTu) * self.__data[0]+ (-uTai/uTu) * self.__data[1]
            Z[1]+= (-uTar/uTu) * self.__data[1]+   uTai/uTu  * self.__data[0]

        else:
            Z= u - self.__data
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

class DcRes(ForwardModel):
    """
    Forward Model for dc resistivity, with a given source pair.
    The cost function is defined as:

        * defect = 1/2 (sum_s sum_pq w_pqs * ((phi_sp-phi_sq)-v_pqs)**2 *

    where p and q indate the 

    """
    def __init__(self, domain, locator, delphi_in, sampleTags, phiPrimary, sigmaPrimary, w=1., coordinates=None, tol=1e-8,saveMemory=True,b=None):
        
        """
        setup new ForwardModel
        :param domain: the domain of the model
        :type: escript domain
        :param locator: contains locator to the measurement pairs
        :type: `list` of ``Locator``
        :param: delphi_in: this is v_pq, the potential difference for the current source  and a set of measurement pairs. a list of measured potential differences is expected. Note this should be the secondary potential only. 
        :type delphi_in: tuple
        :param sampleTags:  tags of measurement points from which potential differences will be calculated.
        :type sampleTags: list of tuples
        :param phiPrimary: primary potential.
        :type phiPrimary: `Scalar`
        """
        super(DcRes, self).__init__()

        self.__domain=domain
        self.__tol = tol
        self.__locator=locator
        self.__trafo=makeTranformation(domain, coordinates) 
        self.__delphi_in=delphi_in
        self.__sampleTags = sampleTags
        self.__sigmaPrimary = sigmaPrimary

        if not isinstance(sampleTags, list):
            raise ValueError("sampleTags must be a list.")    
        if not len(sampleTags) == len(delphi_in):
            raise ValueError("sampleTags and delphi_in must have the same length.")  
        if not len(sampleTags)>0:
            raise ValueError("sampleTags list is empty.")    
        if not isinstance(sampleTags[0], tuple):
            raise ValueError("sampleTags must be a list of tuple.")    

        if isinstance(w, float) or isinstance(w, int):
               w =[ float(w) for z in delphi_in]
               self.__w=w
        if not len(w) == len(delphi_in):
               raise ValueError("Number of confidence factors and number of potential input values don't match.")

        self.__phiPrimary=phiPrimary
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian Coordinates are not supported yet.")
        if not isinstance(locator, Locator):
            raise ValueError("locator must be an escript locator object")    
        
        
        self.__pde=None        
        if not saveMemory:
            self.__pde=self.setUpPDE()

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

    def getPrimaryPotential(self):
        """
        returns the primary potential
        :rtype: `Data`
        """
        return self.__phiPrimary 
    def setUpPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        if self.__pde is None:
            dom=self.__domain
            x = dom.getX()
            DIM=dom.getDim()
            q=whereZero(x[DIM-1]-inf(x[DIM-1]))
            for i in xrange(DIM-1):
                xi=x[i]
                q+=whereZero(xi-inf(xi))+whereZero(xi-sup(xi))
            
            pde=LinearPDE(dom, numEquations=1)
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOn()
            A=pde.createCoefficient('A')
            X=pde.createCoefficient('X')
            pde.setValue(A=A,X=X,q=q)

            # else:
                # pde=LinearPDE(self.__domain, numEquations=1)   
                # A=pde.createCoefficient('A')
                # z = dom.getX()[DIM-1]
                # q=whereZero(z-inf(z))
                # r=0
                # pde.setValue(A=APrimary,y_dirac=y_dirac,d=alpha)

        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        return pde

    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: conductivity
        :type sigma: ``Data`` of shape (1,)
        :return: phi
        :rtype: ``Data`` of shape (1,)
        """
        # print "sigmaPrimary",self.__sigmaPrimary
        # print "sigma=",sigma
        dom=self.__domain
        pde=self.setUpPDE()
        X=(self.__sigmaPrimary - sigma) * grad(self.__phiPrimary)
        # print "+++++++++++++++++++"
        # print "sigma=",sigma
        # print "A=",A
        # print "X=",X
        # print "+++++++++++++++++++"
        pde.setValue(A=sigma * kronecker(dom),X=X)
        phi=pde.getSolution()        
        # print "got U Sol"
        #u.expand()
        # print "u=",u,"grad(u)=",grad(u)
        loc=self.__locator
        loc_phi=loc.getValue(phi)
        # print "read %s at %s"%(str(val),str(loc.getX()))
        return phi, loc_phi

   
    def getDefect(self, sigma, phi, loc_phi):
        """
        Returns the defect value.

        :param sigma: a suggestion for conductivity
        :type sigma: ``Data`` of shape (1,)
        :param phi: potential field
        :type phi: ``Data`` of shape (1,)
        
        :rtype: ``float``
        """
        # print ("getting Defect")
        # print "sigma=",sigma
        # print "placeholder=",placeholder
        val=loc_phi
        # print "val=",val
        length=len(val)
        # print self.__sampleTags[0] 
        if((self.__sampleTags[0][1]!="-" and (length%2) != 0) or (self.__sampleTags[0][1]!="-" and length/2 != len(self.__delphi_in))):
            raise ValueError("length of locator is wrong")

        delphi_calc=[]
        if self.__sampleTags[0][1]!="-":
            for i in range(0,length,2):
                delphi_calc.append(val[i+1]-val[i])
        else:
            for i in range(length):
                delphi_calc.append(val[i])
        A=0
        if (self.__sampleTags[0][1]!="-"):
            for i in range(length/2):
                A+=(self.__w[i]*(delphi_calc[i]-self.__delphi_in[i])**2)        
                # print "delphi_calc[i]=",delphi_calc[i],"self.__delphi_in[i]",self.__delphi_in[i] 
        else:
            for i in range(length):
                A+=(self.__w[i]*(delphi_calc[i]-self.__delphi_in[i])**2)        
                # A+=(self.__w[i]*(self.__delphi_in[i]-delphi_calc[i])**2)        
                # print "delphi_calc[i]=",delphi_calc[i],"self.__delphi_in[i]",self.__delphi_in[i] 
        print ("A/2=",A/2)


        return  A/2

    def getGradient(self, sigma, phi, loc_phi):
        """
        Returns the gradient of the defect with respect to density.

        :param sigma: a suggestison for conductivity
        :type sigma: ``Data`` of shape (1,)
        :param phi: potential field
        :type phi: ``Data`` of shape (1,)
        """
        # print ("getting gradient")
        # print "sigma=",sigma
        val=loc_phi
        sampleTags=self.__sampleTags

        jointSamples={}
        # print(sampleTags)
        for i in range(0,2*len(sampleTags),2): #2*len because sample tags is a list of tuples
            # print(i)
            if sampleTags[i/2][1]!="-":
                tmp=(val[i+1]-val[i]-self.__delphi_in[i/2])*self.__w[i]
            else:
                tmp=(val[i]-self.__delphi_in[i/2]) *self.__w[i]
                # tmp=self.__delphi_in[i/2]-val[i]
            # print ("in gradient","i=", i,"val[i]=",-val[i],"self.__delphi_in[i/2]=",self.__delphi_in[i/2])
            # print ("tmp=",tmp)
            if sampleTags[i/2][0] in jointSamples.keys():
                jointSamples[sampleTags[i/2][0]].append((sampleTags[i/2][1], -tmp))
            else:
                jointSamples[sampleTags[i/2][0]]=[(sampleTags[i/2][1],-tmp)]
            
            if sampleTags[i/2][1]!="-":
                if sampleTags[i/2][1] in jointSamples.keys():
                    jointSamples[sampleTags[i/2][1]].append((sampleTags[i/2][0], tmp))
                else:
                    jointSamples[sampleTags[i/2][1]]=[(sampleTags[i/2][0], tmp)]

        print ("jointSamples=",jointSamples)
        pde =self.setUpPDE()
        dom=self.__domain
        # conPrimary=self.__sigmaPrimary
        # APrimary = conPrimary * kronecker(dom)

        y_dirac = Scalar(0,DiracDeltaFunctions(dom))
        for i in jointSamples:
            total=0
            for j in jointSamples[i]:
                total+=j[1]
            print("setting y_dirac ", i, " to ", total)
            y_dirac.setTaggedValue(i,total)

        pde.setValue(A=sigma*kronecker(dom), y_dirac=y_dirac)
        u=pde.getSolution()
        retVal=-inner(grad(u),grad(phi+self.__phiPrimary))
        return retVal


class MT2DModelTEMode(ForwardModel):
    """
    Forward Model for two dimensional MT model in the TE mode for a given
    frequency omega.
    It defines a cost function:

      *  defect = 1/2 integrate( sum_s w^s * ( E_x - Z_XY^s * H_y ) ) ** 2  *

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
    def __init__(self, domain, omega, x, Z_XY, eta, w0=1., mu=4*PI*1e-7, E_x0=1,
        coordinates=None, fixAtBottom=False, tol=1e-8, saveMemory=True, directSolver=False):
        """
        initializes a new forward model.

        :param data: real and imaginary part of data
        :type data: ``Data`` of shape (2,)
        :param F: real and imaginary part of source given at Dirac points, on surface or at volume.
        :type F: ``Data`` of shape (2,)

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
        :param E_x0: value of E_x at top the domain and if `fixAtBottom` is set at the bottom of the domain.
        :type E_x0: ``float``, ``complex`` or ``Data`` of shape (2,)
        :param mu: permeability
        :type mu: ``float``
        :param coordinates: defines coordinate system to be used (not supported yet)
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param fixAtBottom: if true E_x is set E_x0 at the bottom of the domain
        :type fixAtBottom: ``bool``
        :param saveMemory: if true stiffness matrix is deleted after solution
                           of the PDE to minimize memory requests. This will
                           require more compute time as the matrix needs to be
                           reallocated.
        :type saveMemory: ``bool``
        """
        super(MT2DModelTEMode, self).__init__()
        self.__domain = domain
        DIM=domain.getDim()
        self.__trafo=makeTranformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian coordinates are not supported yet.")
        self.__omega=omega
        self.__mu=mu

        #====================================
        if not len(x) == len(Z_XY):
            raise ValueError("Number of data points and number of impedance values don't match.")
        self.__x=x
        f=-1./(complex(0,1)*omega*mu)
        self.__scaledZxy=[ z*f for z in Z_XY ]
        #====================================
        if isinstance(eta, float) or isinstance(eta, int) :
            self.__eta =[ float(eta) for z in Z_XY]
        else:
            self.__eta =eta
        if not len(self.__eta) == len(Z_XY):
                raise ValueError("Number of confidence radii and number of impedance values don't match.")
        #====================================
        if isinstance(w0, float) or isinstance(w0, int):
               w0 =[ float(w0) for z in Z_XY]
        if not len(w0) == len(Z_XY):
                raise ValueError("Number of confidence factors and number of impedance values don't match.")

        self.__wx0=[0.] * len(Z_XY)
        x=Function(domain).getX()
        totalS=0
        s = 0
        while s < len(self.__scaledZxy):
            totalS+=w0[s]
            f=integrate(self.getWeightingFactor(x, 1., self.__x[s], self.__eta[s]))
            if f < 2*PI*self.__eta[s]**2 * 0.1 :
                raise ValueError("Zero weight (almost) for data point %s. Change eta or refine mesh."%(s,))
            self.__wx0[s]=w0[s]/f
            s += 1
        if not totalS >0 : 
             raise ValueError("Scaling of weight factors failed as sum is zero.")
            
        self.__wx0=[ w/totalS for w in self.__wx0 ]
        #====================================
        if isinstance(E_x0, float) or isinstance(E_x0, int) :
            self.__E_x0 =Data((E_x0,0), Solution(domain))
        elif isinstance(E_x0, tuple):
            self.__E_x0 =Data((E_x0[0],E_x0[1]), Solution(domain))
        elif isinstance(E_x0, complex):
            self.__E_x0 =Data((E_x0.real,E_x0.imag), Solution(domain))
        else:
            if not E_x0.getShape() == (2,):
                raise ValueError("Expected shape of E_x0 is (2,)")
            self.__E_x0= E_x0
        #=======================
        self.__tol=tol
        self.__fixAtBottom=fixAtBottom
        self.__directSolver=directSolver
        self.__pde=None
        if not saveMemory:
            self.__pde=self.setUpPDE()

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
            pde=LinearPDE(self.__domain, numEquations=2)
            if(self.__directSolver == True):
                pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
            D=pde.createCoefficient('D')
            A=pde.createCoefficient('A')
            Y=pde.createCoefficient('Y')
            X=pde.createCoefficient('X')
            A[0,:,0,:]=kronecker(self.__domain.getDim())
            A[1,:,1,:]=kronecker(self.__domain.getDim())
            pde.setValue(A=A, D=D, X=X, Y=Y)
            DIM=self.__domain.getDim()
            z = self.__domain.getX()[DIM-1]
            q=whereZero(z-sup(z))
            if self.__fixAtBottom:
                q+=whereZero(z-inf(z))
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
        :return: E_x, E_x,z
        :rtype: ``Data`` of shape (2,)
        """
        pde=self.setUpPDE()
        D=pde.getCoefficient('D')
        f=self.__omega * self.__mu * sigma
        D[0,1]= - f
        D[1,0]=   f
        pde.setValue(D=D, r=self.__E_x0)
        u=pde.getSolution()        
        return u, grad(u)[:,1]

    def getWeightingFactor(self, x, wx0, x0, eta):
        """
        returns the weighting factor
        """
        return wx0 * exp(length(x-x0)**2*(-0.5/eta**2))

    def getDefect(self, sigma, E_x, dE_xdz):
        """
        Returns the defect value.

        :param sigma: a suggestion for complex 1/V**2
        :type sigma: ``Data`` of shape (2,)
        :param E_x: electric field
        :type E_x: ``Data`` of shape (2,)
        :param dE_xdz: vertical derivative of electric field
        :type dE_xdz: ``Data`` of shape (2,)

        :rtype: ``float``
        """
        x=Function(self.getDomain()).getX()
        A=0
        u0=E_x[0]
        u1=E_x[1]
        u01=dE_xdz[0]
        u11=dE_xdz[1]

        # this cane be done faster!
        s = 0
        while s < len(self.__scaledZxy):
            ws=self.getWeightingFactor(x, self.__wx0[s], self.__x[s], self.__eta[s])
            A+=integrate( ws * ( (u0-u01*self.__scaledZxy[s].real+u11*self.__scaledZxy[s].imag)**2 + (u1-u01*self.__scaledZxy[s].imag-u11*self.__scaledZxy[s].real)**2 ) )
            s += 1
        return  A/2

    def getGradient(self, sigma, E_x, dE_xdz):
        """
        Returns the gradient of the defect with respect to density.

        :param sigma: a suggestion for complex 1/V**2
        :type sigma: ``Data`` of shape (2,)
        :param E_x: electric field
        :type E_x: ``Data`` of shape (2,)
        :param dE_xdz: vertical derivative of electric field
        :type dE_xdz: ``Data`` of shape (2,)
        """
        x=Function(self.getDomain()).getX()
        pde=self.setUpPDE()

        u0=E_x[0]
        u1=E_x[1]
        u01=dE_xdz[0]
        u11=dE_xdz[1]

        D=pde.getCoefficient('D')
        Y=pde.getCoefficient('Y')
        X=pde.getCoefficient('X')

        f=self.__omega * self.__mu * sigma
        D[0,1]=  f
        D[1,0]= -f

        s = 0
        while s < len(self.__scaledZxy):
            ws=self.getWeightingFactor(x, self.__wx0[s], self.__x[s], self.__eta[s])
            Y[0]+=ws * (u0 - u01* self.__scaledZxy[s].real + u11* self.__scaledZxy[s].imag)
            Y[1]+=ws * (u1 - u01* self.__scaledZxy[s].imag - u11* self.__scaledZxy[s].real)
            X[0,1]+=ws * (u01* abs(self.__scaledZxy[s])**2 - u0* self.__scaledZxy[s].real - u1* self.__scaledZxy[s].imag )
            X[1,1]+=ws * (u11* abs(self.__scaledZxy[s])**2 + u0* self.__scaledZxy[s].imag - u1* self.__scaledZxy[s].real )
            s += 1
        pde.setValue(D=D, X=X, Y=Y)
        Zstar=pde.getSolution()
        return self.__omega * self.__mu * (Zstar[0]*u1-Zstar[1]*u0)
        
class Subsidence(ForwardModel):
    """
    Forward Model for subsidence inversion minimizing integrate( (inner(w,u)-d)**2) 
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
        rescales the weights such that
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
            

