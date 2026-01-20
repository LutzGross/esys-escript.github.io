
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""
These domain couplers create a matching domain from each of two different domain
families, allowing more simple interpolation across the two. These domains
must already support interpolation in at least one direction.

"""


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escriptcore.escriptcpp import getMPISizeWorld, MPIBarrierWorld
try:
    from esys.ripley import Rectangle as rRectangle
    from esys.ripley import Brick as rBrick
except ImportError:
    raise ImportError("Missing required ripley module")
try:
    from esys.speckley import Rectangle as sRectangle
    from esys.speckley import Brick as sBrick
except ImportError:
    raise ImportError("Missing required speckley module")
class SpeckleyToRipley(object):
    """
    A class for creating and storing a matching pair of domains from speckley
    and ripley families.
    """

    DEFAULT_lengths = (1.,1.,1.)
    DEFAULT_order = 2

    def __init__(self, dimensions, pointsPerDim, lengths=None, diracPoints=[],
                    diracTags=[], order=None):
        """
        Initialises the coupler, creating domains.

       :param dimensions: whether 2-dimensional or 3-dimensional
       :type dimensions: ``int``
       :param pointsPerDim: the number of data points (not elements) in each
                            dimension
       :type pointsPerDim: ``tuple`` of ``int``
       :param lengths: the length of the domain, defaults to 1 in each dimension
       :type lengths: ``tuple`` of ``int``
       :param diracPoints: set of dirac point locations
       :type diracPoints: ``tuple`` of ``tuple`` of ``float``
       :param diracTags: tag name for each point in `diracPoints`
       :type diracTags: ``tuple`` of ``string``
       :param order: element order of the speckley domain, defaults to 2
       :type order: ``int``
        """
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
        """
        Creates a pair of domains with the previously supplied information and
        number of points.

       :param pointsPerDim: the number of data points (not elements) in each
                            dimension
       :type pointsPerDim: ``tuple`` of ``int``
        """
        
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
        """
        Returns the shape of the domains

       :return: A ``tuple`` of the two shape ``tuples`` in the form
                (speckley, ripley)
       :rtype: ``tuple`` of ``tuple`` of ``int``
        """
        return (self.speckleyShape, self.ripleyShape)

    def getSpeckleyShape(self):
        """
        Returns the shape of the speckley domain

       :return: A ``tuple`` containing the number of elements in each dimension
       :rtype: ``tuple`` of ``int``
        """
        return self.speckleyShape

    def getRipleyShape(self):
        """
        Returns the shape of the ripley domain.

       :return: A ``tuple`` containing the number of elements in each dimension
       :rtype: ``tuple`` of ``int``
        """
        return self.ripleyShape

    def getDomains(self):
        """
        Returns both domains as a tuple.

        :return: A tuple containing (speckley_domain, ripley_domain)
        :rtype: ``tuple``
        """
        return (self.speckley, self.ripley)

    def getSpeckley(self):
        """
        Returns the speckley domain

       :return: The speckley domain
       :rtype: ``SpeckleyDomain``
        """
        return self.speckley

    def getRipley(self):
        """
        Returns the `ripley` domain

       :return: The `ripley` domain
       :rtype: ``RipleyDomain``
        """
        return self.ripley


