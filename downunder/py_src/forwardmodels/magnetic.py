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
from __future__ import division, print_function

"""Forward models for magnetic fields"""

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['MagneticModel', 'SelfDemagnetizationModel', 'MagneticIntensityModel']

from .base import ForwardModelWithPotential
from esys.escript import Scalar
from esys.escript.util import *


class MagneticModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, w, B, background_magnetic_flux_density,
                 coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
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
        background_magnetic_flux_density=interpolate(background_magnetic_flux_density, self.getDataFunctionSpace() )
        if not self.getCoordinateTransformation().isCartesian():
            s = self.getCoordinateTransformation().getScalingFactors()
            v = self.getCoordinateTransformation().getVolumeFactor()
            self.__B_r = background_magnetic_flux_density * s * v
            self.__B_b = background_magnetic_flux_density / s

            A = self.getPDE().createCoefficient("A")
            fw = s**2 * v
            for i in range(self.getDomain().getDim()):
                A[i,i]=fw[i]
            self.getPDE().setValue(A=A)
        else: # cartesian
            self.getPDE().setValue(A=kronecker(self.getDomain()))
            self.__B_r = background_magnetic_flux_density
            self.__B_b = background_magnetic_flux_density

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
    def __init__(self, domain, w, B, background_magnetic_flux_density,
                 coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
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
        background_magnetic_flux_density = interpolate(background_magnetic_flux_density, self.getDataFunctionSpace())
        if not self.getCoordinateTransformation().isCartesian():
            s = self.getCoordinateTransformation().getScalingFactors()
            v = self.getCoordinateTransformation().getVolumeFactor()
            self.__B_r = background_magnetic_flux_density * s * v
            self.__B_b = background_magnetic_flux_density / s

            self.__fw = s**2 * v
        else: # cartesian
            self.__fw = 1
            self.__B_r = background_magnetic_flux_density
            self.__B_b = background_magnetic_flux_density

        # keep track of k used to build PDE:
        self.__last_k = None
        # this is just the initial set_up
        A=self.getPDE().createCoefficient("A")
        for i in range(self.getDomain().getDim()):
            A[i,i]=1.
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
               for i in range(self.getDomain().getDim()):
                   A[i,i] = 1+k
           else:
               for i in range(self.getDomain().getDim()):
                   A[i,i] = (1+k)*self.__fw[i]

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
        pde.setValue(X = k*self.__B_r)
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

class MagneticIntensityModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic intensity inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, w, b, background_magnetic_flux_density,
                 coordinates=None, fixPotentialAtBottom=False, tol=1e-8):
        """
        Creates a new magnetic intensity model on the given domain with one or more
        surveys (w, b).

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Scalar`` or list of ``Scalar``
        :param b: magnetic intensity field data
        :type b: ``Scalar`` or list of ``Scalar``
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        :param background_magnetic_flux_density: background magnetic flux
               density (in Tesla) with components (B_east, B_north, B_vertical)
        :type background_magnetic_flux_density: ``Vector`` or list of `float`
        :param coordinates: defines coordinate system to be used
        :type coordinates: None
        :param fixPotentialAtBottom: if true potential is fixed to zero at the
                                     bottom of the domain in addition to the top
        :type fixPotentialAtBottom: ``bool``
        """
        super(MagneticIntensityModel, self).__init__(domain, w, b, coordinates, fixPotentialAtBottom, tol)
        background_magnetic_flux_density=interpolate(background_magnetic_flux_density, self.getDataFunctionSpace())
        if not self.getCoordinateTransformation().isCartesian():   # need to be checked!
            s = self.getCoordinateTransformation().getScalingFactors()
            v = self.getCoordinateTransformation().getVolumeFactor()
            self.__B_r = background_magnetic_flux_density * s * v
            self.__B_b = background_magnetic_flux_density / s

            A = self.getPDE().createCoefficient("A")
            fw = s**2 * v
            for i in range(self.getDomain().getDim()):
                A[i,i]=fw[i]
            self.getPDE().setValue(A=A)
        else: # cartesian
            self.getPDE().setValue(A=kronecker(self.getDomain()))
            self.__B_r = background_magnetic_flux_density
            self.__B_b = background_magnetic_flux_density
        self.__normalized_B_b=normalize(self.__B_b)

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
        weights=self.getMisfitWeights()
        data=self.getData()
        A=0.
        for s in xrange(len(weights)):
            A += integrate( (weights[s]*(inner(self.__normalized_B_b, magnetic_flux_density)-data[s]) )**2 )
        return A/2

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
        weights=self.getMisfitWeights()
        data=self.getData()
        Y=Scalar(0.,magnetic_flux_density.getFunctionSpace())
        for s in xrange(len(weights)):
            Y+=weights[s]**2*(data[s]-inner(self.__normalized_B_b, magnetic_flux_density))
        pde=self.getPDE()
        pde.resetRightHandSideCoefficients()
        pde.setValue(X=Y*self.__normalized_B_b)
        YT=pde.getSolution()
        return inner(grad(YT),self.__B_r) -Y*inner(self.__normalized_B_b,self.__B_b)
