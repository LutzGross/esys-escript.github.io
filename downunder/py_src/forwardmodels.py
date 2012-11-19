
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

"""Collection of forward models that define the inversion problem"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['ForwardModel','ForwardModelWithPotential','GravityModel','MagneticModel']

from esys.escript import unitsSI as U
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE
from esys.escript.util import *
from esys.escript import Data, Vector, Scalar, Function
from math import pi as PI


class ForwardModel(object):
    """
    An abstract forward model that can be plugged into a cost function.
    Subclasses need to implement `getValue()`, `getGradient()`, and possibly
    `getArguments()`.
    """
    def __init__(self):
        pass

    def getArguments(self, x):
        return ()

    def getValue(self, x, *args):
        raise NotImplementedError

    def getGradient(self, x, *args):
        raise NotImplementedError


class ForwardModelWithPotential(ForwardModel):
    """
    Base class for a forward model using a potential such as magnetic or
    gravity. It defines a cost function

        defect = 1/2 sum_s integrate( weight_i[s] * ( r_i - data_i[s] )**2 )

    where s runs over the survey, weight_i are weighting factors, data_i are
    the data, and r_i are the results produced by the forward model.
    It is assumed that the forward model is produced through postprocessing
    of the solution of a potential PDE.
    """
    def __init__(self, domain, weight, data,  tol=1e-8):
        """
        initialization.

        :param domain: domain of the model
        :type domain: `esys.escript.Domain`
        :param weight: data weighting factors
        :type weight: `Vector` or list of `Vector`
        :param data: data
        :type data: `Vector` or list of `Vector`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``

        """
        super(ForwardModelWithPotential, self).__init__()
        self.__domain = domain

        try:
            n=len(weight)
            m=len(data)
            if m != n:
                raise ValueError("Length of weight and data must be the same.")
            self.__weight = weight
            self.__data = data
        except TypeError:
            self.__weight = [weight]
            self.__data = [data]

        A=0
        for s in xrange(len(self.__weight)):
            A += integrate(inner(self.__weight[s], self.__data[s]**2))
        if A > 0:
	    A=vol(domain)/A 
            for s in xrange(len(self.__weight)):  self.__weight[s]*=A
            print "FORWARD SCALE = ",self.__weight[0]
        else:
            raise ValueError("No non-zero data contribution found.")


        BX = boundingBox(domain)
        DIM = domain.getDim()
        x = domain.getX()
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()
        #self.__pde.setValue(q=whereZero(x[DIM-1]-BX[DIM-1][1]))
        self.__pde.setValue(q=whereZero(x[DIM-1]-BX[DIM-1][1])+whereZero(x[DIM-1]-BX[DIM-1][0]))
    def getDomain(self):
        """
        Returns the domain of the forward model.
        """
        return self.__domain

    def getPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        return self.__pde

    def getDefect(self, result):
        """
        Returns the defect value.

        :param result: a result vector
        :type result: `Vector`
        :rtype: ``float``
        """
        A=0.
        for s in xrange(len(self.__weight)):
            A += integrate(inner(self.__weight[s], (result-self.__data[s])**2))
        return A/2

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
    def __init__(self, domain, chi, g,  gravity_constant=U.Gravitational_Constant, tol=1e-8):
        """
        Creates a new gravity model on the given domain with one or more
        surveys (chi, g).

        :param domain: domain of the model
        :type domain: `esys.escript.Domain`
        :param chi: data weighting factors
        :type chi: `Vector` or list of `Vector`
        :param g: gravity anomaly data
        :type g: `Vector` or list of `Vector`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        """
        super(GravityModel, self).__init__(domain, chi, g, tol)

        self.__G = gravity_constant
        self.getPDE().setValue(A=kronecker(self.getDomain()))

    def getArguments(self, rho):
        """
        Returns precomputed values shared by getValue() and getGradient().

        :param rho: a suggestion for the density distribution
        :type rho: `Scalar`
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
        :type rho: `Scalar`
        :return: gravity potential
        :rtype: ``Scalar``
        """
        pde=self.getPDE()

        pde.resetRightHandSideCoefficients()
        pde.setValue(Y=-4.*PI*self.__G*rho)
        phi=pde.getSolution()

        return phi

    def getValue(self, rho, phi, gravity_force):
        """
        Returns the value of the defect

        :param rho: density distribution
        :type rho: `Scalar`
        :param phi: corresponding potential
        :type phi: `Scalar`
        :param gravity_force: gravity force
        :type gravity_force: `Vector`
        :rtype: ``float``
        """
        return self.getDefect(gravity_force)

    def getGradient(self, rho, phi, gravity_force):
        """
        Returns the gradient of the defect with respect to density.

        :param rho: density distribution
        :type rho: `Scalar`
        :param phi: corresponding potential
        :type phi: `Scalar`
        :param gravity_force: gravity force
        :type gravity_force: `Vector`
        :rtype: `Scalar`
        """
        pde=self.getPDE()
        g, chi = self.getSurvey()

        Z=0.
        for s in xrange(len(chi)):
            Z = chi[s] * (g[s]-gravity_force) + Z

        pde.resetRightHandSideCoefficients()
        pde.setValue(X=Z)
        ZT=pde.getSolution()
        return ZT*(-4*PI*self.__G)


class MagneticModel(ForwardModelWithPotential):
    """
    Forward Model for magnetic inversion as described in the inversion
    cookbook.
    """
    def __init__(self, domain, chi, B, background_field,  tol=1e-8):
        """
        Creates a new magnetic model on the given domain with one or more
        surveys (chi, B).

        :param domain: domain of the model
        :type domain: `Domain`
        :param chi: data weighting factors
        :type chi: `Vector` or list of `Vector`
        :param B: magnetic field data
        :type B: `Vector` or list of `Vector`
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        """
        super(MagneticModel, self).__init__(domain, chi, B, tol)
        print "B = ",B[0]
        print "background_field = ",background_field
        
        self.__background_field=interpolate(background_field, B[0].getFunctionSpace())
        self.getPDE().setValue(A=kronecker(self.getDomain()))

    def getArguments(self, k):
        """
        Returns precomputed values shared by getValue() and getGradient().

        :param k: susceptibility
        :type k: `Scalar`
        :return: scalar magnetic potential and corresponding magnetic field
        :rtype: ``Scalar``, ``Vector``
        """
        phi = self.getPotential(k)
        magnetic_field = k * self.__background_field -grad(phi)
        return phi, magnetic_field

    def getPotential(self, k):
        """
        Calculates the magnetic potential for a given susceptibility.

        :param k: susceptibility
        :type k: `Scalar`
        :return: magnetic potential
        :rtype: ``Scalar``
        """
        pde=self.getPDE()

        pde.resetRightHandSideCoefficients()
        pde.setValue(X = k* self.__background_field)
        phi=pde.getSolution()

        return phi

    def getValue(self, k, phi, magnetic_field):
        """
        Returns the value of the defect.

        :param k: susceptibility
        :type k: `Scalar`
        :param phi: corresponding potential
        :type phi: `Scalar`
        :param magnetic_field: magnetic field
        :type magnetic_field: `Vector`
        :rtype: ``float``
        """
        return self.getDefect(magnetic_field)

    def getGradient(self, k, phi, magnetic_field):
        """
        Returns the gradient of the defect with respect to susceptibility.

        :param k: susceptibility
        :type k: `Scalar`
        :param phi: corresponding potential
        :type phi: `Scalar`
        :param magnetic_field: magnetic field
        :type magnetic_field: `Vector`
        :rtype: `Scalar`
        """
        pde=self.getPDE()
        B, chi = self.getSurvey()

        Y=0.
        for s in xrange(len(chi)):
            Y = chi[s] * (magnetic_field-B[s]) + Y

        pde.resetRightHandSideCoefficients()
        pde.setValue(X=Y)

        #print "test VTK file generated"
        #from esys.weipa import saveVTK
        #saveVTK("Y.vtu",RHS=pde.getRightHandSide(), Y=Y)

        YT=pde.getSolution()
        print "YT = ",YT
        return inner(Y-grad(YT),self.__background_field)

