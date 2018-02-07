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

"""Forward model for gravity (Bouguer) anomaly"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['GravityModel']

from .base import ForwardModelWithPotential
from esys.escript import unitsSI as U
from esys.escript.util import *
from math import pi as PI


class GravityModel(ForwardModelWithPotential):
    """
    Forward Model for gravity inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, w, g, gravity_constant=U.Gravitational_Constant,
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
        :param fixPotentialAtBottom: if true potential is fixed to zero at the
                                     base of the domain in addition to the top
        :type fixPotentialAtBottom: ``bool``

        :note: It is advisable to call rescaleWeights() to rescale weights
               before starting the inversion.
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
        pde = self.getPDE()
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

