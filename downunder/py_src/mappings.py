##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

"""Collection of parametrizations that map physical values to model parameters
   and back"""

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['Mapping', 'DensityMapping', 'SusceptibilityMapping',\
           'BoundedRangeMapping', 'LinearMapping', 'AcousticVelocityMapping',\
           'MTMapping']

from esys.escript import inf, sup, log, tanh, boundingBoxEdgeLengths, clip, atan2, sin, cos, sqrt, exp, whereZero
import esys.escript.unitsSI as U
from numpy import pi
import logging

class Mapping(object):
    """
    An abstract mapping class to map level set functions *m* to physical
    parameters *p*.
    """

    def __init__(self, *args):
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)

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

    def __init__(self, a=1., p0=0.):
        self.__a=a
        self.__p0=p0
        self.__a_inv = 1./a

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

class AcousticVelocityMapping(Mapping):
    """
    Maps a p-velocity and Q-index to slowness square sigma=(V*(1-i*1/(2*Q))^{-2}
    in the form sigma=e^{Mr+m[0])}*( cos(Mi+m[1])) + i * sin(Mi+m[1])
    """
    def __init__(self, V_prior, Q_prior):
        """
        initializes the mapping

        :param V_prior: a-priori p-wave velocity
        :param Q_prior: a-priori Q-index (must be positive)
        """
        over2Q=1./(2*Q_prior)
        # sigma_prior=1/(V_prior*(1-I*over2Q))**2 = 1/( V_prior * (1+over2Q**2)) **2 * ( (1-over2Q**2) + I * 2* over2Q )
        self.Mr=log( sqrt((1-over2Q**2)**2+(2*over2Q)**2)/(V_prior*(1+over2Q**2))**2 )
        self.Mi=atan2(2*over2Q, 1-over2Q**2)

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        return exp(m[0]+self.Mr)*(cos(m[1]+self.Mi)*[1,0]+ sin(m[1]+self.Mi)*[0,1])

    def getDerivative(self, m):
        """
        returns the value for the derivative of the mapping for m
        """
        e=exp(m[0]+self.Mr)
        return (e*cos(m[1]+self.Mi))*[[1,0],[0,1]]+(e*sin(m[1]+self.Mi))*[[0,-1],[1,0]]

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for s
        """
        # self.logger.info("m0:"+str(log(s[0]**2+s[1]**2)/2))
        return (log(s[0]**2+s[1]**2)/2-self.Mr)*[1., 0 ] + (atan2(s[1],s[0])-self.Mi)*[0, 1. ]

class DcResMapping(Mapping):
    """DcResMapping
        sigmoid mapping
        s=a/(1+e^(-k*m))    
    """
    def __init__(self, sigma_prior, k=1., a=0.01, minVal=1/1000.):
        self.__sigma0=sigma_prior
        self.__k=k
        self.a=a
        self.minVal=minVal
    def getValue(self, m):
        print ("in get value inf(m)=",inf(m)," sup(m)=", sup(m))
        # s=self.__sigma0 + (self.__sigma0 * self.__k*m)
        # s=self.__sigma0*exp(self.__k*m)
        #### use sigmoid mapping
        s=(self.__sigma0*self.a / (1+exp(-self.__k * m))) + self.minVal
        print ("in get value inf(s)=",inf(s)," sup(s)=", sup(s))
        if sup(s)!=0 and inf(s)!=0:
            print ("in get value 1/inf(s)=",1./inf(s)," 1/sup(s)=", 1./sup(s))
        return s
       
    def getDerivative(self, m):
        """
        returns the derivative of the mapping for m
        """
        # return self.__sigma0 * self.__k
        # return self.__sigma0*self.__k*exp(self.__k*m)
        return (self.__sigma0*self.__k*self.a*exp(-self.__k*m))/ (1+exp(-self.__k * m))**2

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for s
        """
        # return (s-self.__sigma0) / (self.__sigma0 * self.__k)
        if inf(((self.__sigma0*self.a)/s)) <= 1.:
            raise ValueError("sigma 0*a/s < 1 this is not valid as log cannot be 0 or negative")
        m= - 1./self.__k * log(((self.__sigma0*self.a)/(s-self.minVal))-1)
        print ("inv(s)=",m)
        return m 
        

class MTMapping(Mapping):
    """
    mt mapping

    sigma=sigma0*exp(a*m)
    """

    def __init__(self,sigma_prior, a=1.):
        """
        initializes the mapping

        :param sigma_prior: a a-priori conductivity
        """
        self.__sigma0 = sigma_prior
        self.__a = a

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        return self.__sigma0 * exp(self.__a*m)

    def getDerivative(self, m):
        """
        returns the derivative of the mapping for m
        """
        return self.__sigma0*self.__a*exp(self.__a*m)

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for s
        """
        ms=whereZero(s)
        ms0=whereZero(self.__sigma0)
        m=1/self.__a* log(((1-ms)*s+ms*1)/((1-ms0)*self.__sigma0+ms0*1)) * (1-ms)*(1-ms0) 
        return m

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

