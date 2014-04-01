
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

__all__ = ['ForwardModel','ForwardModelWithPotential','GravityModel','MagneticModel', 'SelfDemagnetizationModel', 'AcousticWaveForm']

from esys.escript import unitsSI as U
from esys.escript import Data, Vector, Scalar, Function, DiracDeltaFunctions, FunctionOnBoundary
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE, SolverOptions
from .coordinates import makeTranformation
from esys.escript.util import *
from math import pi as PI
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

        if not self.getCoordinateTransformation().isCartesian():
            self.__G = 4*PI*gravity_constant * self.getCoordinateTransformation().getVolumeFactor()

            fw=self.getCoordinateTransformation().getScalingFactors()**2*self.getCoordinateTransformation().getVolumeFactor()
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
        :param background_magnetic_flux_density: background magnetic flux density (in Teslar) with components (B_east, B_north, B_vertical)
        :type background_magnetic_flux_density: ``Vector`` or list of `float`
        :param coordinates: defines coordinate system to be used
        :type coordinates: ReferenceSystem` or `SpatialCoordinateTransformation`
        :param fixPotentialAtBottom: if true potential is fixed to zero at the bottom of the domain
                                     in addition to the top.
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
        
class  SelfDemagnetizationModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic inversion with self-demagnetization as described in the inversion
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
        :param background_magnetic_flux_density: background magnetic flux density (in Teslar) with components (B_east, B_north, B_vertical)
        :type background_magnetic_flux_density: ``Vector`` or list of `float`
        :param coordinates: defines coordinate system to be used
        :type coordinates: ReferenceSystem` or `SpatialCoordinateTransformation`
        :param fixPotentialAtBottom: if true potential is fixed to zero at the bottom of the domain
                                     in addition to the top.
        :type fixPotentialAtBottom: ``bool``
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

        # keep track on k used to build PDE:
        self.__last_k=None
        # this is just the initial set_up
        A=self.getPDE().createCoefficient("A")
        for i in range(self.getDomain().getDim()): A[i,i]=1
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
        if not self.__last_k == k:
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
    Forward Model for acoustic waveform inversion in the frequence domain 
    It defines a cost function:

        defect = 1/2 integrate( ( w * ( a * u - data ) ) ** 2 )

    where w are weighting factors, data are the measured data (as a 2-comp vector of real and imaginary part)  for real frequency omega, 
    and u is the coresponding result produced by the forward model. u (as a 2-comp vector) is the solution of the 
    complex Helmholtz equation for frequency omega, source F and complex, inverse, squared p-velocity sigma:
       * -u_{ii} - omega**2 * sigma * u = F
    It is assumed that the exact scale of source F is unknown and the scaling factor a of F is calculated by minimizing the
    defect 
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
        :param F: real and imaginary part of source given at Dirac points, on surface or at volume.
        :type F: ``Data`` of shape (2,)
        :param coordinates: defines coordinate system to be used (not supported yet)
        :type coordinates: ReferenceSystem` or `SpatialCoordinateTransformation`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param saveMemory: if true stiffness matrix is deleted after solution of PDE to 
                           minimize memory requests. This will require more compute time as
                           the matrix needs to be reallocated. 
        :type saveMemory: ``bool``
        :param scaleF: if true source F is scaled to minimize defect. 
        :type scaleF: ``bool``
        :param fixAtBottom: if true pressure is fixed to zero at the bottom of the domain 
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
                raise ValueError("data must be escript.Data object.")
        if not data.getFunctionSpace() == FunctionOnBoundary(domain):
                raise ValueError("data must be defined on boundary")                    
        if not data.getShape() == (2,):
                raise ValueError("data must have shape 2 (real and imaginary part).")
        if w == None: 
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
                raise ValueError("Sourcemust have shape 2 (real and imaginary part).")
                
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
        Return the underlying PDE.

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
        # assummuing  integrate(w * length(data)**2) =1 
        if self.scaleF and abs(uTu) >0:
           A=1.-(uTar**2 + uTai**2)/uTu  
        else:
           A =  integrate(self.__weight*length(self.__data)**2)- 2 * uTar + uTu
        return  A/2

    def getGradient(self,sigma, u,uTar, uTai, uTu ):
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





                
               
