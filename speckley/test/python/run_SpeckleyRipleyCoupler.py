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

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.speckley import Rectangle as sRectangle, Brick as sBrick

#domain families may not be present
HAS_RIPLEY = True
HAS_FINLEY = True
try:
    from esys.ripley import Rectangle as rRectangle, Brick as rBrick
except ImportError as e:
    HAS_RIPLEY = False
    
try:
    from esys.finley import Rectangle as fRectangle
except ImportError as e:
    HAS_FINLEY = False

if HAS_RIPLEY:
    from esys.escriptcore.domainCouplers import SpeckleyToRipley


@unittest.skipIf(not HAS_RIPLEY, "Ripley domains not present")
class Test_ripleyCoupler(unittest.TestCase):

    def calculateVariance(self, s, r):
        actual = r.getX()
        sinput = s.getX()
        sX = interpolate(sinput, Function(s))
        return actual - interpolate(sX, Function(r)) #actual - interpo...

    def badInterpolations(self, speckley, ripley):
        FS = Function(speckley)
        FR = Function(ripley)
        #bad speck -> good rip
        with self.assertRaises(RuntimeError):
            interpolate(speckley.getX(), FR)
        with self.assertRaises(RuntimeError):
            interpolate(Data(5, ReducedFunction(speckley)), FR)
        #good speck -> bad rip
        with self.assertRaises(RuntimeError):
            interpolate(Data(5, FS), ReducedFunction(ripley))
        with self.assertRaises(RuntimeError):
            interpolate(Data(5, FS), ContinuousFunction(ripley))

    def test_Rectangle_non_Function(self):
        for order in range(2, 11):
            coupler = SpeckleyToRipley(2, (2*getMPISizeWorld()*order,order),
                    order=order, lengths=[3.*getMPISizeWorld(),2.])
            self.badInterpolations(coupler.getSpeckley(), coupler.getRipley())

    def test_Brick_non_Function(self):
        for order in range(2, 11):
            #values here are arbitrary, just has to be Bricks
            coupler = SpeckleyToRipley(3, (2*getMPISizeWorld()*order,order,order),
                    order=order, lengths=[3.*getMPISizeWorld(),2.,2.])
            self.badInterpolations(coupler.getSpeckley(), coupler.getRipley())


    def test_Rectangle(self):
        for order in range(2,11):
            coupler = SpeckleyToRipley(2, (2*getMPISizeWorld()*order,order),
                    order=order, lengths=[3.*getMPISizeWorld(),2.])
            d = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            res = Lsup(d)
            self.assertLess(res, 1e-10,
                    "Rectangles of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%res)

            coupler = SpeckleyToRipley(2, (order, 2*getMPISizeWorld()*order),
                    order=order, lengths=[2.,3.*getMPISizeWorld()])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Rectangles of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%Lsup(a))

    @unittest.skipIf(getMPISizeWorld() < 4, "requires at least 4 ranks")
    def test_MultiDimSplit_Rectangle(self):
        ranks = getMPISizeWorld()
        squares = [4,9,16,25,36] #could do more, unlikely to test though
        subdivs = [0,0]
        loop = 0
        if ranks not in squares:
            if ranks % 3 == 0:
                subdivs[0] = 3
                subdivs[1] = ranks//3
                loop = 2
            elif ranks % 2 == 0:
                subdivs[0] = 2
                subdivs[1] = ranks//2
                loop = 1
            else:
                raise unittest.SkipTest("inappropriate number of ranks")
        else:
            subdivs = [squares.index(ranks)+2]*2
            loop = 2

        for order in range(2,11):
            divs = subdivs[:]
            for i in range(loop):
                r = rRectangle(2*divs[0]*order-1, 2*divs[1]*order-1,
                        d0=divs[0], d1=divs[1])
                s = sRectangle(order, divs[0]*2, divs[1]*2,
                        d0=divs[0], d1=divs[1])
                d = self.calculateVariance(s, r)
                res = Lsup(d)
                self.assertLess(res, 1e-10,
                        "".join(["Rectangles of order %d with "%order,
                                "subdivisons %s "%divs,
                                "failed to interpolate correctly",
                                ", variance of %g"%res]))
                divs.append(divs.pop(0))

    @unittest.skipIf(getMPISizeWorld() < 4, "requires at least 4 ranks")
    def test_BiaxialSplits_Brick(self):
        ranks = getMPISizeWorld()
        squares = [4,9,16,25,36] #could do more, unlikely to test though
        error_msg = "Brick of order {0} and subdivisions {1} failed to"+\
                " interpolate correctly, variance of {2}"
        subdivs = [0,0,0]
        loop = 3
        if ranks not in squares:
            if ranks % 3 == 0:
                subdivs[0] = 3
                subdivs[1] = ranks//3
                subdivs[2] = 1
            elif ranks % 2 == 0:
                subdivs[0] = 2
                subdivs[1] = ranks//2
                subdivs[2] = 1
            else:
                raise unittest.SkipTest("inappropriate number of ranks")
        else:
            subdivs = [squares.index(ranks)+2]*2 + [1]
            loop = 1

        for order in range(2,11):
            for i in range(loop):
                divs = subdivs[i:] + subdivs[:i]
                r = rBrick(2*divs[0]*order-1, 2*divs[1]*order-1, 2*divs[2]*order-1,
                        d0=divs[0], d1=divs[1], d2=divs[2])
                s = sBrick(order, divs[0]*2, divs[1]*2, divs[2]*2,
                        d0=divs[0], d1=divs[1], d2=divs[2])
                d = self.calculateVariance(s, r)
                res = Lsup(d)
                self.assertLess(res, 1e-10, error_msg.format(order, divs, res))

    @unittest.skipIf(getMPISizeWorld() < 8, "requires at least 8 ranks")
    def test_TriaxialSplits_Brick(self):
        ranks = getMPISizeWorld()
        cubes = [8,27,64] #could do more, unlikely to test though
        error_msg = "Brick of order {0} and subdivisions {1} failed to"+\
                " interpolate correctly, variance of {2}"
        subdivs = [0,0,0]
        loop = 3
        if ranks not in cubes:
            if ranks % 2 == 0:
                subdivs[0] = 2
                if (ranks // 2) % 2 == 0:
                    subdivs[1] = 2
                    subdivs[2] = ranks // 4
                elif (ranks // 2) % 3 == 0:
                    subdivs[1] = 3
                    subdivs[2] = ranks // 6
                else:
                    raise unittest.SkipTest("inappropriate number of ranks")
            elif ranks % 3 == 0:
                subdivs[0] = 3
                if (ranks // 2) % 2 == 0:
                    subdivs[1] = 2
                    subdivs[2] = ranks // 4
                elif (ranks // 2) % 3 == 0:
                    subdivs[1] = 3
                    subdivs[2] = ranks // 6
                else:
                    raise unittest.SkipTest("inappropriate number of ranks")
            else:
                raise unittest.SkipTest("inappropriate number of ranks")
        else:
            subdivs = [cubes.index(ranks)+2]*3
            loop = 1

        for order in range(2,11):
            for i in range(loop):
                divs = subdivs[i:] + subdivs[:i]
                r = rBrick(2*divs[0]*order-1, 2*divs[1]*order-1, 2*divs[2]*order-1,
                        d0=divs[0], d1=divs[1], d2=divs[2])
                s = sBrick(order, divs[0]*2, divs[1]*2, divs[2]*2,
                        d0=divs[0], d1=divs[1], d2=divs[2])
                res = Lsup(self.calculateVariance(s, r))
                self.assertLess(res, 1e-10, error_msg.format(order, divs, res))

    def test_Brick(self):
        elements = [3*getMPISizeWorld(), 2, 1]
        lengths = [3.*getMPISizeWorld(), 4., 5.]
        labels = ["x", "y", "z"]
        error_message = "Brick {0}-split of order {1} failed to interpolate " +\
                        "correctly, variance of {2}"
        for i in range(len(labels)):
            e = tuple(elements[i:] + elements[:i])
            l = tuple(lengths[i:] + lengths[:i])
            for order in range(2,11):
                ele = tuple([j*order for j in e])
                coupler = SpeckleyToRipley(3, ele, order=order, lengths=l)
                a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
                self.assertLess(Lsup(a), 1e-10,
                        error_message.format(labels[i], order, Lsup(a)))

    def test_mismatch_errors(self):
        ranks = getMPISizeWorld()
        r = rRectangle(2*ranks - 1, 2, l0=1., d0=ranks)
        s = sRectangle(2, 2*ranks, 2, l0=2., d0=ranks)
        with self.assertRaises(RuntimeError): #length mismatch
            self.calculateVariance(r, s)
        r = rBrick(2*ranks - 1, 2, 2, l1=2., d0=ranks)
        s = sBrick(2, 2*ranks, 2, 2, l1=1., d0=ranks)
        with self.assertRaises(RuntimeError): #length mismatch
            self.calculateVariance(r, s)
        r = rRectangle(2*ranks - 1, 2, l0=1., l1=2., d0=ranks)
        s = sBrick(2, 2*ranks, 2, 2, l0=2., l1=1., d0=ranks)
        with self.assertRaises(RuntimeError): #dimensionality mismatch
            self.calculateVariance(r, s)
        if getMPISizeWorld() > 1:
            r = rRectangle(2*ranks - 1, 2, l0=1., d0=ranks)
            s = sRectangle(2, 2, 2*ranks, l0=2., d1=ranks)
            with self.assertRaises(RuntimeError): #subdivision mismatch
                self.calculateVariance(r, s)

    @unittest.skipIf(not HAS_FINLEY, "Finley domains not present")
    def test_bad_domains(self):
        ranks = getMPISizeWorld()
        f = fRectangle(2*ranks - 1, 2, l0=1., l1=2., d0=ranks)
        s = sRectangle(2, 2*ranks, 2, l0=2., d0=ranks)
        with self.assertRaises(RuntimeError): #finley is not ripley
            self.calculateVariance(f, s)

    def test_creation(self):
        with self.assertRaises(RuntimeError): #invalid order
            coupler = SpeckleyToRipley(2, pointsPerDim=(4,2*getMPISizeWorld()), order=1)
        with self.assertRaises(RuntimeError): #invalid order
            coupler = SpeckleyToRipley(2, pointsPerDim=(40,20*getMPISizeWorld()), order=11)
        with self.assertRaises(ValueError): #invalid dimension
            coupler = SpeckleyToRipley(1, pointsPerDim=(4,2*getMPISizeWorld()))
        with self.assertRaises(ValueError): #incalid dimension
            coupler = SpeckleyToRipley(4, pointsPerDim=(4,2*getMPISizeWorld()))


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

