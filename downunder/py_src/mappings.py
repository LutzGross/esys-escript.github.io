
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

"""Collection of parametrizations that map physical values to model parameters
   and back"""

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['Mapping', 'DensityMapping', 'SusceptibilityMapping', 'BoundedRangeMapping', 'LinearMapping']

from esys.escript import inf, sup, log, tanh, boundingBoxEdgeLengths, clip
import esys.escript.unitsSI as U

class Mapping(object):
    """
    An abstract mapping class to map level set functions *m* to physical
    parameters *p*.
    """

    def __init__(self, *args):
        pass

    def __call__(self, m):
        """
        short for getValue(m).
        """
        return self.getValue(m)

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        raise NotImplementedError

    def getDerivative(self, m):
        """
        returns the value for the derivative of the mapping for m
        """
        raise NotImplementedError

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for physical parameter p
        """
        raise NotImplementedError

    def getTypicalDerivative(self):
        """
        returns a typical value for the derivative
        """
        raise NotImplementedError


class LinearMapping(Mapping):
    """
    Maps a parameter by a linear transformation p = a * m + p0
    """

    def __init__(self, a=1, p0=0):
        self.__a=a
        self.__p0=p0
        self.__a_inv = 1/a

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        return self.__a*m + self.__p0

    def getDerivative(self, m):
        """
        returns the value for the derivative of the mapping for m
        """
        return self.__a

    def getInverse(self, p):
        """
        returns the value of the inverse of the mapping for s
        """
        return self.__a_inv * ( p - self.__p0)

    def getTypicalDerivative(self):
        """
        returns a typical value for the derivative
        """
        return self.__a

class DensityMapping(LinearMapping):
    """
    Density mapping with depth weighting

    *rho =  rho0 + drho  * ( (x_2 - z0)/l_z)^(beta/2) ) * m*

    """
    def __init__(self, domain, z0=None, rho0=None, drho=None,  beta=None):
        """
        initializes the mapping

        :param domain: domain of the mapping
        :type domain: ``Domain``
        :param z0: depth weighting offset. If not present no depth scaling is applied.
        :type z0: scalar
        :param rho0: reference density, defaults to 0
        :type rho0: scalar
        :param drho: density scale. By default density of granite = 2750kg/m**3 is used.
        :type drho: scalar
        :param beta: depth weighting exponent, defaults to 2
        :type beta: ``float``
        """
        if rho0 is None: rho0=0.
        if drho is None: drho=2750*U.kg/U.m**3
        if beta is None: beta=2.

        self.domain=domain
        if z0 is not None:
            DIM=self.domain.getDim()
            l_z=boundingBoxEdgeLengths(domain)[DIM-1]
            a = drho * ( clip(z0-domain.getX()[DIM-1], minval=0)/l_z )**(beta/2)
        else:
            a = drho
        super(DensityMapping,self).__init__(a=a, p0=rho0)

class SusceptibilityMapping(LinearMapping):
    """
    Susceptibility mapping with depth weighting

    *k =  k0 + dk  * ( (x_2 - z0)/l_z)^(beta/2) ) * m*

    """
    def __init__(self, domain, z0=None, k0=None, dk=None,  beta=None):
        """
        set up mapping

        :param domain: domain of the mapping
        :type domain: ``Domain``
        :param z0: depth weighting offset. If not present no depth scaling is applied.
        :type z0: scalar
        :param k0: reference density, defaults to 0
        :type k0: scalar
        :param dk: susceptibility scale, defaults to 1
        :type dk: scalar
        :param beta: depth weighting exponent, defaults to 2
        :type beta: ``float``
        """
        if k0 is None: k0=0.
        if beta is None: beta=2.
        if dk is None: dk=1.
        self.domain=domain

        if z0 is not None:
            DIM=self.domain.getDim()
            l_z=boundingBoxEdgeLengths(domain)[DIM-1]
            a = dk * ( clip(z0-domain.getX()[DIM-1] , minval=0)/l_z )**(beta/2)
        else:
            a = dk
        super(SusceptibilityMapping,self).__init__(a=a, p0=k0)

# needs REVISION
class BoundedRangeMapping(Mapping):
    """
    Maps an unbounded parameter to a bounded range. The mapping is smooth and
    continuous.
    """

    def __init__(self, s_min=0, s_max=1):
        if not s_min < s_max:
            raise ValueError("value for s_min must be less than the value for s_max.")
        self.s_min=s_min
        self.s_max=s_max

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        return (self.s_max+self.s_min)/2. + (self.s_max-self.s_min)/2. * tanh(m)

    def getDerivative(self, m):
        """
        returns the value for the derivative of the mapping for m
        """
        return ((self.s_max-self.s_min)/2.) * (1.-tanh(m)**2.)

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for s
        """
        if not (inf(s) > self.s_min and sup(s) < self.s_max):
            raise ValueError("s is out of range [%f,%f]"%(inf(s),sup(s)))

        return 1./2. * log( (s-self.s_min)/(self.s_max-s) )

