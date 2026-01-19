
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os, sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.escript.linearPDEs import LameEquation
from esys.ripley import MultiResolutionDomain
from run_diracOnRipley import Test_RipleyDiracPoints

mpiSize = getMPISizeWorld()

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='.'

brickLevel = 2
rectLevel = 2

def Rectangle(**kwargs):
    m = MultiResolutionDomain(2, **kwargs)
    return m.getLevel(rectLevel - 1)

def Brick(**kwargs):
    m = MultiResolutionDomain(3, **kwargs)
    return m.getLevel(brickLevel - 1)

@unittest.skipIf(mpiSize > 1, "Multiresolution domains require single process")
class Test_DiracPointsOnMultiResolutionDomains(Test_RipleyDiracPoints):

    def setup(self):
        # constants
        self.numRanks = getMPISizeWorld()
        self.rank = getMPIRankWorld()
        self.shortEdge = 3
        self.longFactor = 5
        self.longEdge = self.longFactor*self.numRanks-1
        self.empty = "(data contains no samples)\n"

    def getRectRefs(self, xLong):
        Ex = self.longEdge+1
        Ey = self.shortEdge+1
        if not xLong:
            Ex, Ey = Ey, Ex
        result = [[-1 for j in range(Ex)] for i in range(Ey)]
        ref = 0
        if xLong:
            for rankx in range(self.numRanks):
                for y in range(Ey):
                    for x in range(self.longFactor):
                        old = [ref%Ex, ref//Ex]
                        new = [i*rectLevel for i in old]
                        node = new[0] + new[1]*(rectLevel*Ex-1)
                        result[y][x+self.longFactor*rankx] = node
                        ref += 1
        else:
            for y in range(Ey):
                for x in range(Ex):
                    old = [ref%Ex, ref//Ex]
                    new = [i*rectLevel for i in old]
                    node = new[0] + new[1]*(rectLevel*Ex-1)
                    result[y][x] = node
                    ref += 1
        return result
            
    def getBrickRefs(self, splitAxis, dims):
        dims = [i+1 for i in dims]
        results = [[[-1 for z in range(dims[2])] for y in range(dims[1])] for x in range(dims[0])]
        ref = 0
        rankDim = [i for i in dims]
        rankDim[splitAxis] = dims[splitAxis]//self.numRanks
        rc = [0, 0, 0] #rank counters
        for rank in range(self.numRanks):
            for z in range(rankDim[2]):
                for y in range(rankDim[1]):
                    for x in range(rankDim[0]):
                        old = [ref%dims[0], (ref%(dims[0]*dims[1]))//dims[0], ref//(dims[0]*dims[1])]
                        new = [i*brickLevel for i in old]
                        node = new[0] + new[1]*(brickLevel*(dims[0]-1)+1) + new[2]*(brickLevel*(dims[0]-1) + 1)*(brickLevel*(dims[1]-1)+1)
                        results[x+rc[0]][y+rc[1]][z+rc[2]] = node
                        ref += 1
            rc[splitAxis] += rankDim[splitAxis]
        return results

    def generateRects(self, a, b):
        diracLoc = [a,b]
        edges = [self.longEdge, self.shortEdge]
        rects = []
        for i in range(2):
            rects.append(Rectangle(n0=edges[0], n1=edges[1],
                      l0=edges[0], l1=edges[1],
                      d0=self.numRanks, diracPoints=[tuple(diracLoc)],
                      diracTags=["test"]))
            diracLoc = diracLoc[::-1]
            edges = edges[::-1]
        return rects

    def generateBricks(self, a, b, c):
        diracLoc = [a,b,c]
        bricks = []
        edges = [self.longEdge, self.shortEdge, self.shortEdge]
        for i in range(3):
            bricks.append(Brick(n0=edges[0], n1=edges[1], n2=edges[2],
                l0=edges[0], l1=edges[1], l2=edges[2],
                d0=self.numRanks,
                diracPoints=[tuple(diracLoc)], diracTags=["test"]))
            diracLoc = diracLoc[2:] + diracLoc[:2]
            edges = edges[2:] + edges[:2]
        tmp = [self.shortEdge]*3
        dims = [tmp[:], tmp[:], tmp[:]]
        for i in range(3):
            dims[i][i] = self.longEdge
        return bricks, dims

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)


