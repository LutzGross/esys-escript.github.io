
##############################################################################
#
# Copyright (c) 2012-2018 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2012-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import logging
import numpy as np
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
from esys.downunder.minimizers import *
from esys.downunder.costfunctions import CostFunction

# number of dimensions for the test function
N=3

# this is mainly to avoid warning messages
lslogger=logging.getLogger('inv').setLevel(logging.CRITICAL)

# Rosenbrock test function to be minimized. The minimum is 0 and lies at
# [1,1,...,1].
class RosenFunc(CostFunction):
    def __init__(self):
        super(RosenFunc, self).__init__()
    def getDualProduct(self, f0, f1):
        return np.dot(f0, f1)
    def getNorm(self,x):
        return (abs(x.max()))
    def getGradient(self, x, *args):
        xm = x[1:-1]
        xm_m1 = x[:-2]
        xm_p1 = x[2:]
        der = np.zeros_like(x)
        der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
        der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
        der[-1] = 200*(x[-1]-x[-2]**2)
        return der
    def getValue(self, x, *args):
        return np.sum(100.0*(x[1:]-x[:-1]**2.)**2. + (1-x[:-1])**2.)

class TestMinimizerLBFGS(unittest.TestCase):
    def setUp(self):
        self.f=RosenFunc()
        self.minimizer=MinimizerLBFGS(self.f)
        self.x0=np.array([2.]*N)
        self.xstar=np.array([1.]*N)

    def test_max_iterations(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(1)
        self.assertRaises(MinimizerMaxIterReached, self.minimizer.run,self.x0)

    def test_solution(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(120)
        x=self.minimizer.run(self.x0)
        xx=self.minimizer.getResult()
        self.assertEqual(np.amax(abs(x-xx)), 0.)
        # We should be able to get a solution in under 100 iterations
        self.assertLess(np.amax(abs(x-self.xstar)), 1e-7)

    def test_callback(self):
        n=[0]
        def callback(**args):
            n[0]=n[0]+1
        self.minimizer.setCallback(callback)
        self.minimizer.setTolerance(1e-8)
        self.minimizer.setMaxIterations(10)
        try:
            x=self.minimizer.run(self.x0)
        except MinimizerMaxIterReached:
            pass
        # callback should be called once for each iteration (including 0th)
        self.assertEqual(n[0], 11)

    def test_solution_with_interpolation_order1(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(120)
        self.minimizer.setOptions(interpolationOrder=1)
        x=self.minimizer.run(self.x0)

        self.assertLess(np.amax(abs(x-self.xstar)), 1e-7)
        
    def test_solution_with_interpolation_order2(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(200)
        self.minimizer.setOptions(interpolationOrder=2)
        
        xx=self.minimizer.run(self.x0)

        self.assertLess(np.amax(abs(xx-self.xstar)), 1e-7)
        
    def test_solution_with_interpolation_order3(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(200)
        self.minimizer.setOptions(interpolationOrder=3)
        xxx=self.minimizer.run(self.x0)
        self.assertLess(np.amax(abs(xxx-self.xstar)), 1e-7)
    
    def test_solution_with_interpolation_order4(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(200)
        self.minimizer.setOptions(interpolationOrder=4)
        xxxx=self.minimizer.run(self.x0)
        self.assertLess(np.amax(abs(xxxx-self.xstar)), 1e-7)


class TestMinimizerNLCG(unittest.TestCase):
    def setUp(self):
        self.f=RosenFunc()
        self.minimizer=MinimizerNLCG(self.f)
        self.x0=np.array([2.]*N)
        self.xstar=np.array([1.]*N)

    def test_max_iterations(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(1)
        self.assertRaises(MinimizerMaxIterReached, self.minimizer.run, self.x0)

    def test_solution(self):
        self.minimizer.setTolerance(1e-8)
        self.minimizer.setMaxIterations(2000)
        x=self.minimizer.run(self.x0)
        xx=self.minimizer.getResult()
        self.assertEqual(np.amax(abs(x-xx)), 0.)
        # We should be able to get a solution to set tolerance in #iterations
        self.assertAlmostEqual(np.amax(abs(x-self.xstar)), 0., places=3)

    def test_callback(self):
        n=[0]
        def callback(**args):
            n[0]=n[0]+1
        self.minimizer.setCallback(callback)
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(10)
        try:
            x=self.minimizer.run(self.x0)
        except MinimizerMaxIterReached:
            pass
        # callback should be called once for each iteration (including 0th)
        self.assertEqual(n[0], 11)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

