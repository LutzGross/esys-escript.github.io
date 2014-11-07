
##############################################################################
#
# Copyright (c) 2014 by University of Queensland
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

"""
Domain couplers allow coupling of two domains, each from a different domain
family. This coupling enables the interpolation from one to another.

"""

__copyright__="""Copyright (c) 2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escriptcore.escriptcpp import getMPISizeWorld, MPIBarrierWorld
from esys.ripley import Rectangle as rRectangle
from esys.ripley import Brick as rBrick
from esys.speckley import Rectangle as sRectangle
from esys.speckley import Brick as sBrick

class SpeckleyToRipley(object):
    DEFAULT_lengths = (1.,1.,1.)
    DEFAULT_order = 2

    def __init__(self, dimensions, pointsPerDim, lengths=None, diracPoints=[],
                    diracTags=[], order=None):
        self.ranks = getMPISizeWorld()

        if dimensions not in [2,3]:
            raise ValueError("SpeckleyToRipley: requires dimension of 2 or 3")
        self.dim = dimensions

        if len(pointsPerDim) != dimensions:
            raise ValueError("SpeckleyToRipley: requires point estimate for each dimension")
        for point in pointsPerDim:
            if point < 2:
                raise ValueError("SpeckleyToRipley: requires at least 2 points per dimension")

        if lengths is None:
            lengths = self.DEFAULT_lengths[:self.dim]
        else:
            if len(lengths) != dimensions:
                raise ValueError("SpeckleyToRipley: requires length for each dimension")
            for length in lengths:
                if length <= 0:
                    raise ValueError("SpeckleyToRipley: requires positive lengths")
        self.lengths = lengths

        if order is None:
            order = self.DEFAULT_order
        self.order = order
        self.diracPoints = diracPoints
        self.diracTags = diracTags

        self.createDomains(pointsPerDim)


    def createDomains(self, pointsPerDim):
        splitDim = pointsPerDim.index(max(pointsPerDim))
        divs = [1]*self.dim
        divs[splitDim] = self.ranks
        speckleyElements = []
        ripleyElements = []
        for i in range(self.dim):
            points = pointsPerDim[i]
            if i == splitDim:
                ripleyElements.append(points + self.ranks - (points % self.ranks) - 1)
            else:
                ripleyElements.append(points)
            speck = ripleyElements[i]/self.order % self.ranks
            if i == splitDim:
                if speck != 0:
                    speck = self.ranks - speck
            speck += ripleyElements[i]/self.order
            if speck < 2:
                speck = 2
            speckleyElements.append(speck)

        self.speckleyShape = tuple(speckleyElements)
        self.ripleyShape = tuple(ripleyElements)
        if self.dim == 2:
            l0,l1 = self.lengths
            ex,ey = speckleyElements
            self.speckley = sRectangle(self.order, ex, ey,
                    d0=divs[0], d1=divs[1], l0=l0, l1=l1,
                    diracPoints=self.diracPoints, diracTags=self.diracTags)

            ex,ey = ripleyElements
            self.ripley = rRectangle(ex, ey,
                    d0=divs[0], d1=divs[1], l0=l0, l1=l1,
                    diracPoints=self.diracPoints, diracTags=self.diracTags)
        else:
            l0,l1,l2 = self.lengths
            ex,ey,ez = speckleyElements
            self.speckley = sBrick(self.order, ex, ey, ez,
                    d0=divs[0], d1=divs[1], d2=divs[2], l0=l0, l1=l1, l2=l2,
                    diracPoints=self.diracPoints, diracTags=self.diracTags)

            ex,ey,ez = ripleyElements
            self.ripley = rBrick(ex, ey, ez,
                    d0=divs[0], d1=divs[1], d2=divs[2], l0=l0, l1=l1, l2=l2,
                    diracPoints=self.diracPoints, diracTags=self.diracTags)

    def getShapes(self):
        return (self.speckleyShape, self.ripleyShape)

    def getSpeckleyShape(self):
        return self.speckleyShape

    def getRipleyShape(self):
        return self.ripleyShape

    def getDomains(self):
        return (self.speckley, self.ripley)

    def getSpeckley(self):
        return self.speckley

    def getRipley(self):
        return self.ripley


