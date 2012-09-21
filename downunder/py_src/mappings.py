
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

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['Mapping', 'BoundedRangeMapping', 'ScalingMapping']

from esys.escript import inf, sup, log, tanh

class Mapping(object):
    """
    An abstract mapping class.
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
        returns the value of the inverse of the mapping for s
        """
        raise NotImplementedError


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


class ScalingMapping(Mapping):
    """
    Maps a parameter by scaling it with a constant.
    """

    def __init__(self, alpha=1):
        if alpha==0:
            raise ValueError("Scaling parameter must not be 0.")
        self.__alpha=alpha

    def getValue(self, m):
        """
        returns the value of the mapping for m
        """
        return self.__alpha*m

    def getDerivative(self, m):
        """
        returns the value for the derivative of the mapping for m
        """
        return self.__alpha

    def getInverse(self, s):
        """
        returns the value of the inverse of the mapping for s
        """
        return s/self.__alpha

