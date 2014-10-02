
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

class Test_ripleyCoupler(unittest.TestCase):

    def test_Rectangle(self):
        for order in range(2,11):
            coupler = SpeckleyToRipley(2, (2*order,order), order=order,
                    lengths=[2.*order,4.*order])
            s = coupler.getSpeckley()
            r = coupler.getRipley()
            sinput = s.getX()
            rinput = r.getX()
            sX = interpolate(sinput, Function(s))
            a = interpolate(sX, Function(r))
            self.assertLess(Lsup(a - interpolate(rinput, Function(r))), 1e-10,
                    "Rectangles of order %d failed to interpolate correctly"%order)

    def test_Brick(self):
        for order in range(2,11):
            coupler = SpeckleyToRipley(3, (4*order, 2*order,order),
                    order=order, lengths=[2.*order,4.*order,8.*order])
            s = coupler.getSpeckley()
            r = coupler.getRipley()
            sinput = s.getX()
            rinput = r.getX()
            sX = interpolate(sinput, Function(s))
            a = interpolate(sX, Function(r))
            self.assertLess(Lsup(a - interpolate(rinput, Function(r))), 1e-10,
                    "Bricks of order %d failed to interpolate correctly"%order)

    def test_creation(self):
        with self.assertRaises(RuntimeError):
            coupler = SpeckleyToRipley(2, pointsPerDim=(4,2), order=1)
        with self.assertRaises(RuntimeError):
            coupler = SpeckleyToRipley(2, pointsPerDim=(40,20), order=11)
        with self.assertRaises(ValueError):
            coupler = SpeckleyToRipley(1, pointsPerDim=(4,2))
        with self.assertRaises(ValueError):
            coupler = SpeckleyToRipley(4, pointsPerDim=(4,2), order=1)
        

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

