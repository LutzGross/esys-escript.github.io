
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.escriptcore.domainCouplers import SpeckleyToRipley
from esys.ripley import Rectangle as rRectangle, Brick as rBrick
from esys.speckley import Rectangle as sRectangle, Brick as sBrick
from esys.finley import Rectangle as fRectangle

class Test_ripleyCoupler(unittest.TestCase):

    def calculateVariance(self, s, r):
        actual = r.getX()
        sinput = s.getX()
        sX = interpolate(sinput, Function(s))
        return actual - interpolate(sX, Function(r)) #actual - interpo...

    def test_Rectangle(self):
        for order in range(2,11):
            coupler = SpeckleyToRipley(2, (2*getMPISizeWorld()*order,order),
                    order=order, lengths=[3.*getMPISizeWorld(),2.])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Rectangles of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%Lsup(a))
            
            coupler = SpeckleyToRipley(2, (order, 2*getMPISizeWorld()*order),
                    order=order, lengths=[2.,3.*getMPISizeWorld()])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Rectangles of order %d failed to interpolate correctly"%order)

    @unittest.skipIf(getMPISizeWorld() < 4, "requires at least 4 ranks")
    def test_MultiDimSplit_Rectangle(self):
        ranks = getMPISizeWorld()
        squares = [4,9,16,25,36] #could do more, unlikely to test though
        subdivs = [0,0]
        loop = 0
        if ranks not in squares:
            if ranks % 2 == 0:
                subdivs[0] = 2
                subdivs[1] = ranks//2
                loop = 1
            elif ranks % 3 == 0:
                subdivs[0] = 3
                subdivs[1] = ranks//3
                loop = 2
            else:
                raise unittest.SkipTest("inappropriate number of ranks")
        else:
            subdivs = [squares.index(ranks)+1]*2

        for order in range(2,11):
            divs = subdivs[:]        
            for i in range(loop):
                r = rRectangle(2*divs[0]*order-1, 2*divs[1]*order-1,
                        d0=divs[0], d1=divs[1])
                s = sRectangle(order, divs[0]*2, divs[1]*2,
                        d0=divs[0], d1=divs[1])
                self.assertLess(Lsup(self.calculateVariance(s, r)), 1e-10,
                        "".join(["Rectangles of order %d with "%order,
                                "subdivisons %s "%subdivs,
                                "failed to interpolate correctly"]))
                divs.append(divs.pop(0))
                
    def test_Brick(self):
        for order in range(2,11):
            coupler = SpeckleyToRipley(3, (3*getMPISizeWorld()*order, 2*order, order),
                    order=order, lengths=[3.*getMPISizeWorld(),4.,5.])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Brick x-split of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%Lsup(a))
            
            coupler = SpeckleyToRipley(3, (order, 3*getMPISizeWorld()*order, 2*order),
                    order=order, lengths=[5.,3.*getMPISizeWorld(),4.])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Brick y-split of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%Lsup(a))

            coupler = SpeckleyToRipley(3, (2*order, order, 3*getMPISizeWorld()*order),
                    order=order, lengths=[4.,5.,3.*getMPISizeWorld()])
            a = self.calculateVariance(coupler.getSpeckley(), coupler.getRipley())
            self.assertLess(Lsup(a), 1e-10,
                    "Brick z-split of order %d failed to interpolate correctly"%order +\
                    ", variance of %g"%Lsup(a))


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

