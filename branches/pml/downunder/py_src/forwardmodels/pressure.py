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

"""Isostatic Pressure calculation"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['IsostaticPressure']

from esys.downunder.coordinates import makeTransformation
from esys.escript import unitsSI as U
from esys.escript import Scalar, Vector, Function
from esys.escript.linearPDEs import LinearSinglePDE
from esys.escript.util import *
from math import pi as PI


class IsostaticPressure(object):
    """
    class to calculate isostatic pressure field correction due to gravity forces
    """
    def __init__(self, domain, p0=0., level0=0, gravity0=-9.81*U.m*U.sec**(-2),
                 background_density=2670* U.kg*U.m**(-3),
                 gravity_constant=U.Gravitational_Constant,
                 coordinates=None, tol=1e-8):
        """
        :param domain: domain of the model
        :type domain: `Domain`
        :param p0: pressure at level0
        :type p0: scalar `Data` or ``float``
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
        self.__trafo=makeTransformation(domain, coordinates)
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()

        z = domain.getX()[DIM-1]
        self.__pde.setValue(q=whereNonNegative(z-level0), r=p0)

        fw = self.__trafo.getScalingFactors()**2 * self.__trafo.getVolumeFactor()
        A=self.__pde.createCoefficient("A")
        for i in range(DIM): A[i,i]=fw[i]
        self.__pde.setValue(A=A)
        z = Function(domain).getX()[DIM-1]
        self.__g_b= 4*PI*gravity_constant/self.__trafo.getScalingFactors()[DIM-1]*background_density*(level0-z) + gravity0
        self.__rho_b=background_density

    def getPressure(self, g = None, rho=None):
        """
        return the pressure for gravity force anomaly `g` and
        density anomaly `rho`

        :param g: gravity anomaly data
        :type g: ``Vector``
        :param rho: gravity anomaly data
        :type rho: ``Scalar``
        :return: pressure distribution
        :rtype: ``Scalar``
        """
        if not g: g=Vector(0., Function(self.__domain))
        if not rho: rho=Scalar(0., Function(self.__domain))

        g2=(rho * self.__g_b)*[0,0,1] + self.__rho_b*g + rho*g
        # Tests need to be updated before the following is uncommented:
        #g2=((rho+self.__rho_b) * self.__g_b)*[0,0,1] + self.__rho_b*g + rho*g
        d=self.__trafo.getScalingFactors()
        V= self.__trafo.getVolumeFactor()
        self.__pde.setValue(X = -g2*d*V)
        #self.__pde.setValue(X = g2*d*V)
        return self.__pde.getSolution()

