
##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
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

from __future__ import division, print_function

from esys.ripley import MultiRectangle, MultiBrick

class MultiResolutionDomain(object):
    def __init__(self, dim, **kwargs):
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
        escriptworld = self.__kwargs.get('escriptWorld', None)
        return MultiRectangle(n0, n1, l0, l1, d0, d1, diracPoints, tags,
                escriptworld, subdivisions)

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
        escriptworld = self.__kwargs.get('escriptWorld', None)
        return MultiBrick(n0, n1, n2, l0, l1, l2, d0, d1, d2, diracPoints, tags,
                escriptworld, subdivisions)

    def getMaxDepth(self):
        return len(self.__levels)
    
    def getLevel(self, level):
        if int(level) != level or level < 0:
            raise ValueError("level must be a non-negative integer")
        dom = self.__levels.get(level)
        if not dom:
            dom = self.__generateDomain(level)
        return dom
