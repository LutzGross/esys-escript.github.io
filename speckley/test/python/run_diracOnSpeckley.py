
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
from esys.speckley import Rectangle, Brick

try:
     SPECKLEY_WORKDIR=os.environ['SPECKLEY_WORKDIR']
except KeyError:
     SPECKLEY_WORKDIR='.'

class Test_SpeckleyDiracPoints(unittest.TestCase):
    # constants
    numRanks = getMPISizeWorld()
    rank = getMPIRankWorld()
    shortEdge = 5
    longEdge = 3*numRanks

    def generateRects(self, order, a, b):
        rectX = Rectangle(order, self.longEdge, self.shortEdge, l0=self.longEdge,
                l1=self.shortEdge, d0=self.numRanks, diracPoints=[(a,b)], diracTags=["test"])
        rectY = Rectangle(order, self.shortEdge, self.longEdge, l0=self.shortEdge,
                l1=self.longEdge, d1=self.numRanks, diracPoints=[(b,a)], diracTags=["test"])
        return [rectX, rectY]

    def generateBricks(self, order, a, b, c):
        brickX = Brick(order, self.longEdge, 5, 5, 
                l0=self.longEdge, l1=5, l2=5,
                d0=self.numRanks,
                diracPoints=[(a,b,c)], diracTags=["test"])
        brickY = Brick(order, 5, self.longEdge, 5, 
                l0=5, l1=self.longEdge, l2=self.shortEdge,
                d1=self.numRanks,
                diracPoints=[(c,a,b)], diracTags=["test"])
        brickZ = Brick(order, 5, 5, self.longEdge,
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

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Creation(self):
        r = self.numRanks
        el = self.numRanks*3

        #test bad types
        with self.assertRaises(TypeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=(.0,0.), diracTags=["test"])
        with self.assertRaises(TypeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[(.0,0.)], diracTags=("test"))
        with self.assertRaises(TypeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[.0], diracTags=["test"])
        with self.assertRaises(TypeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[.0,.0], diracTags=["test"])

        with self.assertRaises(TypeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=(.0,0.,0.), diracTags=["test"])
        with self.assertRaises(TypeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=[(.0,0.,0.)], diracTags=("test"))
        with self.assertRaises(TypeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=[.0,0.], diracTags=["test"])
            
        #test bad arg lengths
        with self.assertRaises(RuntimeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[(.0,)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[(.0,1.)], diracTags=[])
        with self.assertRaises(RuntimeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[(.0,0.)], diracTags=["test", "break"])
        with self.assertRaises(RuntimeError):
            Rectangle(2, 5, el, d1=r,  diracPoints=[(.0,0.,0.)], diracTags=["test"])

        with self.assertRaises(RuntimeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=[(.0,0.,0.,0.)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=[(.0,0.,)], diracTags=["test"])
        with self.assertRaises(RuntimeError):
            Brick(2, 5, el, 5, d1=r, diracPoints=[(.0,)], diracTags=["test"])

    def test_DDF_to_Continuous_2D(self):
        expected_value = self.numRanks*11
        for order in range(2, 11):
            rectX, rectY = self.generateRects(order, self.longEdge, self.shortEdge)
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

    def test_DDF_to_Continuous_3D(self):
        expected_value = self.numRanks*11
        for order in range(2, 11):
            doms, dims = self.generateBricks(order, self.longEdge,
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


