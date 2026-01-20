
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


from esys.ripley import MultiRectangle, MultiBrick

class MultiResolutionDomain(object):
    """ Constructs domains of varying resolutions that are guaranteed to be
    compatible for cross-domain interpolation. The parameters supplied will be
    used to construct the coarsest resolution. No coarser domain can be
    constructed.
    
    Each domain of finer resolution will have the number of elements in every
    axis of the coarsest domain multiplied by ``2**n``, where ``n`` is the level of
    subdivision.
    """
    def __init__(self, dim, **kwargs):
        """
        :param dim: the spatial dimension of the domain to create
        :type dim: `int`
        :param kwargs: the arguments normally passed to a constructor of
                Rectangle or Brick, including as the number of elements ``n0=...``, ``n1=...``, etc.
        """
        self.__kwargs = kwargs
        self.__levels = {}
        self.__dim = dim
        self.__generateDomain(0)

    def __generateDomain(self, level):
        if self.__dim == 2:
            self.__levels[level] = self.__newRectangle(2**level)
        elif self.__dim == 3:
            self.__levels[level] = self.__newBrick(2**level)
        else:
            raise ValueError("Invalid spatial dimension of domain: %d"%self.__dim)
        return self.__levels[level]
    
    def __newRectangle(self, subdivisions):
        n0 = self.__kwargs['n0']
        n1 = self.__kwargs['n1']
        l0 = self.__kwargs.get('l0', 1.)
        l1 = self.__kwargs.get('l1', 1.)
        d0 = self.__kwargs.get('d0', -1)
        d1 = self.__kwargs.get('d1', -1)
        diracPoints = self.__kwargs.get('diracPoints', [])
        tags = self.__kwargs.get('diracTags', [])
        return MultiRectangle(n0, n1, l0, l1, d0, d1, diracPoints, tags, subdivisions)

    def __newBrick(self, subdivisions):
        n0 = self.__kwargs['n0']
        n1 = self.__kwargs['n1']
        n2 = self.__kwargs['n2']
        l0 = self.__kwargs.get('l0', 1.)
        l1 = self.__kwargs.get('l1', 1.)
        l2 = self.__kwargs.get('l2', 1.)
        d0 = self.__kwargs.get('d0', -1)
        d1 = self.__kwargs.get('d1', -1)
        d2 = self.__kwargs.get('d2', -1)
        diracPoints = self.__kwargs.get('diracPoints', [])
        tags = self.__kwargs.get('diracTags', [])
        return MultiBrick(n0, n1, n2, l0, l1, l2, d0, d1, d2, diracPoints, tags, subdivisions)

    def getMaxDepth(self):
        """ Returns the level of the finest domain created so far """
        return len(self.__levels) - 1
    
    def getLevel(self, level):
        """ Returns a domain with each element subdivided ``level`` times 
        
        :param level: the number of times to subdivide each element
        :type level: `int`
        """
        if int(level) != level or level < 0:
            raise ValueError("level must be a non-negative integer")
        dom = self.__levels.get(level)
        if not dom:
            dom = self.__generateDomain(level)
        return dom
