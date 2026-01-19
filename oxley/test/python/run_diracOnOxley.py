
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
from esys.oxley import Rectangle, Brick

try:
     OXLEY_WORKDIR=os.environ['OXLEY_WORKDIR']
except KeyError:
     OXLEY_WORKDIR='.'

class Test_OxleyDiracPoints(unittest.TestCase):

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
                    for x in range(3):
                        result[y][x+3*rankx] = ref
                        ref += 1
        else:
            for y in range(Ey):
                for x in range(Ex):
                    result[y][x] = ref
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
                        results[x+rc[0]][y+rc[1]][z+rc[2]] = ref
                        ref += 1
            rc[splitAxis] += rankDim[splitAxis]
        return results

    def generateRects(self, a, b):
        rectX = Rectangle(self.longEdge, self.shortEdge, l0=self.longEdge,
                l1=self.shortEdge, d0=self.numRanks, diracPoints=[(a,b)], diracTags=["test"])
        rectY = Rectangle(self.shortEdge, self.longEdge, l0=self.shortEdge,
                l1=self.longEdge, d1=self.numRanks, diracPoints=[(b,a)], diracTags=["test"])
        return [rectX, rectY]

    def generateBricks(self, a, b, c):
        brickX = Brick(self.longEdge, 5, 5, 
                l0=self.longEdge, l1=5, l2=5,
                d0=self.numRanks,
                diracPoints=[(a,b,c)], diracTags=["test"])
        brickY = Brick(5, self.longEdge, 5, 
                l0=5, l1=self.longEdge, l2=self.shortEdge,
                d1=self.numRanks,
                diracPoints=[(c,a,b)], diracTags=["test"])
        brickZ = Brick(5, 5, self.longEdge,
                l0=5, l1=5, l2=self.longEdge, 
                d2=self.numRanks,
                diracPoints=[(b,c,a)], diracTags=["test"])
        dims = [ [self.longEdge, 5, 5], [5, self.longEdge, 5], 
                      [5, 5, self.longEdge]]
        return [brickX, brickY, brickZ], dims

    def owns(self, refs, y, x, xLong):
        if xLong and self.rank*3 - 1 <= x <= self.rank*3 + 3:
            return True
        elif not xLong and self.rank*3 - 1 <= y <= self.rank*3 + 3:
            return True
        return False

    def rectMessage(self, loc, xLong):
        if xLong:
            if not (0 <= loc[0] <= self.longEdge and 0 <= loc[1] <= self.shortEdge):
                return self.empty
        else:
            if not (0 <= loc[1] <= self.longEdge and 0 <= loc[0] <= self.shortEdge):
                return self.empty
        x = int(round(loc[0] - 0.01))
        y = int(round(loc[1] - 0.01))
        refs = self.getRectRefs(xLong)
        if self.owns(refs, y, x, xLong):
            return "( id: 0, ref: {0}, pnt: 0) (0) {1}\n( id: 0, ref: {0}, pnt: 0) (1) {2}".format(refs[y][x], x, y)
        return self.empty

    def brickMessage(self, loc, splitAxis, refs, dims):
        dims = [i + 1 for i in dims[splitAxis]]
        for n in range(3):
            if not 0 <= loc[n] <= dims[n]:
                return self.empty
        loc = [int(round(i - 1e-8)) for i in loc]
        if self.rank*3 - 1 <= loc[splitAxis] <= self.rank*3 + 3:
            x, y, z = loc
            return "( id: 0, ref: {0}, pnt: 0) (0) {1}\n( id: 0, ref: {0}, pnt: 0) (1) {2}\n( id: 0, ref: {0}, pnt: 0) (2) {3}".format(refs[x][y][z], x, y, z)
        return self.empty

    def setUp(self):
        # constants
        self.numRanks = getMPISizeWorld()
        self.rank = getMPIRankWorld()
        self.shortEdge = 5
        self.longFactor = 3
        self.longEdge = self.longFactor*self.numRanks-1
        self.empty = "(data contains no samples)\n"

    def tearDown(self):
        pass

    def test_Creation(self):
        r = self.numRanks
        el = self.numRanks*3-1

        #test bad types
        with self.assertRaises(TypeError):
            Rectangle(5, el, d1=r,  diracPoints=(.0,0.), diracTags=["test"])
        with self.assertRaises(TypeError):
            Rectangle(5, el, d1=r,  diracPoints=[(.0,0.)], diracTags=("test"))
        with self.assertRaises(TypeError):
            Rectangle(5, el, d1=r,  diracPoints=[.0], diracTags=["test"])
        with self.assertRaises(TypeError):
            Rectangle(5, el, d1=r,  diracPoints=[.0,.0], diracTags=["test"])

        with self.assertRaises(TypeError):
            Brick(5, el, 5, d1=r, diracPoints=(.0,0.,0.), diracTags=["test"])
        with self.assertRaises(TypeError):
            Brick(5, el, 5, d1=r, diracPoints=[(.0,0.,0.)], diracTags=("test"))
        with self.assertRaises(TypeError):
            Brick(5, el, 5, d1=r, diracPoints=[.0,0.], diracTags=["test"])
            
        #test bad arg lengths
        with self.assertRaises(RuntimeError):
            Rectangle(5, el, d1=r,  diracPoints=[(.0,)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Rectangle(5, el, d1=r,  diracPoints=[(.0,1.)], diracTags=[])
        with self.assertRaises(RuntimeError):
            Rectangle(5, el, d1=r,  diracPoints=[(.0,0.)], diracTags=["test", "break"])
        with self.assertRaises(RuntimeError):
            Rectangle(5, el, d1=r,  diracPoints=[(.0,0.,0.)], diracTags=["test"])

        with self.assertRaises(RuntimeError):
            Brick(5, el, 5, d1=r, diracPoints=[(.0,0.,0.,0.)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Brick(5, el, 5, d1=r, diracPoints=[(.0,0.,)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Brick(5, el, 5, d1=r, diracPoints=[(.0,)], diracTags=["test"])

    @unittest.skip("Oxley Brick with Dirac points causes segfault in renumberNodes() - see issue #118")
    def test_BrickInterpolation(self):
        for a in range(-1, (self.longEdge+self.numRanks)*2, self.numRanks*2):
            a = a//2.
            c = 2.5
            for b in range(self.shortEdge):
                doms, dims = self.generateBricks(a,b,c)
                for n, dom in enumerate(doms):
                    points = [a,b,c]
                    for i in range(n): #rotate points to match
                        points = [points.pop()] + points
                    refs = self.getBrickRefs(n, dims[n])
                    i = interpolate(dom.getX(), DiracDeltaFunctions(dom))
                    expected = self.brickMessage(tuple(points), n, refs, dims)
                    got = str(i)
                    global_result = getMPIWorldSum(1 if got != expected else 0)
                    self.assertEqual(got, expected, 
                            "{0} not mapped correctly on rank {1} for d{2} splits".format(tuple(points), self.rank, n) + \
                            "\nexpected === \n{0}\ngot === \n{1}".format(expected,got)+ \
                            "\nlongedge=%d, shortedge=%d\n"%(self.longEdge, self.shortEdge))
                    #remaining ranks must also exit, otherwise we'll lock up
                    self.assertEqual(global_result, 0, "One or more ranks failed")

    @unittest.skip("Oxley accepts out-of-bounds Dirac points - see issue #118")
    def test_RectangleInterpolation(self):
        for a in range(-1, (self.longEdge+self.numRanks)*2, self.numRanks*2):
            a = a//2.
            for b in range(self.shortEdge):
                rectX, rectY = self.generateRects(a,b)
                i = interpolate(rectX.getX(), DiracDeltaFunctions(rectX))
                got = str(i)
                refs = self.getRectRefs(False)
                expected = self.rectMessage((a,b),True)
                global_result = getMPIWorldSum(1 if got != expected else 0)
                self.assertEqual(got, expected,
                        "({0},{1}) not mapped correctly on rank {2} for d0 splits".format(a,b, self.rank) + \
                        "\nexpected === \n{0}\ngot === \n{1}".format(expected, got)+\
                        "\nlongedge=%d, shortedge=%d\n"%(self.longEdge, self.shortEdge))
                #remaining ranks must also exit, otherwise we'll lock up
                self.assertEqual(global_result, 0, "One or more ranks failed")
                i = interpolate(rectY.getX(), DiracDeltaFunctions(rectY))
                expected = self.rectMessage((b,a),False)
                got = str(i)
                global_result = getMPIWorldSum(1 if got != expected else 0)
                self.assertEqual(got, expected, 
                        "({0},{1}) not mapped correctly on rank {2} for d1 splits".format(b,a, self.rank) + \
                        "\nexpected === \n{0}\ngot === \n{1}".format(expected, got)+\
                        "\nlongedge={0}, shortedge={1}\n".format(self.longEdge, self.shortEdge))
                #remaining ranks must also exit, otherwise we'll lock up
                self.assertEqual(global_result, 0, "One or more ranks failed")

    def test_DDF_to_Continuous_2D(self):
        expected_value = self.numRanks*11
        rectX, rectY = self.generateRects(self.longEdge,self.shortEdge)
        for dom, point in [(rectX, (self.longEdge, self.shortEdge)),
                           (rectY, (self.shortEdge, self.longEdge))]:
            xDDF = Data(0, DiracDeltaFunctions(dom))
            xDDF.setTaggedValue("test", expected_value)

            X = dom.getX()
            expected_data = whereZero(X[0]-point[0]) * whereZero(X[1]-point[1]) * expected_value
            cont = interpolate(xDDF, ContinuousFunction(dom))
            result = Lsup(expected_data - cont)
            self.assertLess(result, 1e-15,
                    "Interpolation failure, expected zero, got %g"%result)

    @unittest.skip("Oxley Brick with Dirac points causes segfault in renumberNodes() - see issue #118")
    def test_DDF_to_Continuous_3D(self):
        expected_value = self.numRanks*11
        doms, dims = self.generateBricks(self.longEdge,
                self.shortEdge, self.shortEdge)
        for dom, point in zip(doms, dims):
            xDDF = Data(0, DiracDeltaFunctions(dom))
            xDDF.setTaggedValue("test", expected_value)

            X = dom.getX()
            expected_data = whereZero(X[0]-point[0]) \
                    * whereZero(X[1]-point[1]) * whereZero(X[2]-point[2]) \
                    * expected_value
            cont = interpolate(xDDF, ContinuousFunction(dom))
            result = Lsup(expected_data - cont)
            self.assertLess(result, 1e-15,
                    "Interpolation failure, expected zero, got %g"%result)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)


