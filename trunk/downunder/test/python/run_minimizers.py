
##############################################################################
#
# Copyright (c) 2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import logging
import numpy as np
import unittest
import sys
from esys.downunder.minimizers import *
from esys.downunder.costfunctions import CostFunction

# number of dimensions for the test function
N=10

# this is mainly to avoid warning messages
logger=logging.getLogger('inv')
logger.setLevel(logging.INFO)
handler=logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

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
        self.minimizer.setTolerance(1e-8)
        self.minimizer.setMaxIterations(100)
        reason=self.minimizer.run(self.x0)
        x=self.minimizer.getResult()
        # We should be able to get a solution in under 100 iterations
        self.assertEqual(reason, MinimizerLBFGS.TOLERANCE_REACHED)
        self.assertAlmostEqual(np.amax(abs(x-self.xstar)), 0.)

    def test_callback(self):
        n=[0]
        def callback(k, x, fg, gf):
            n[0]=n[0]+1
        self.minimizer.setCallback(callback)
        self.minimizer.setTolerance(1e-8)
        self.minimizer.setMaxIterations(10)
        try:
            reason=self.minimizer.run(self.x0)
        except MinimizerMaxIterReached:
            pass
        # callback should be called once for each iteration (including 0th)
        self.assertEqual(n[0], 11)

class TestMinimizerBFGS(unittest.TestCase):
    def setUp(self):
        self.f=RosenFunc()
        self.minimizer=MinimizerBFGS(self.f)
        self.x0=np.array([2.]*N)
        self.xstar=np.array([1.]*N)

    def test_max_iterations(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(1)
        reason=self.minimizer.run(self.x0)
        self.assertEqual(reason, MinimizerBFGS.MAX_ITERATIONS_REACHED)

    def test_solution(self):
        self.minimizer.setTolerance(1e-6)
        self.minimizer.setMaxIterations(100)
        self.minimizer.setOptions(initialHessian=1e-3)
        reason=self.minimizer.run(self.x0)
        x=self.minimizer.getResult()
        # We should be able to get a solution in under 100 iterations
        self.assertEqual(reason, MinimizerBFGS.TOLERANCE_REACHED)
        self.assertAlmostEqual(np.amax(abs(x-self.xstar)), 0.)

    def test_callback(self):
        n=[0]
        def callback(k, x, fg, gf):
            n[0]=n[0]+1
        self.minimizer.setCallback(callback)
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(10)
        reason=self.minimizer.run(self.x0)
        # callback should be called once for each iteration (including 0th)
        self.assertEqual(n[0], 11)

class TestMinimizerNLCG(unittest.TestCase):
    def setUp(self):
        self.f=RosenFunc()
        self.minimizer=MinimizerNLCG(self.f)
        self.x0=np.array([2.]*N)
        self.xstar=np.array([1.]*N)

    def test_max_iterations(self):
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(1)
        reason=self.minimizer.run(self.x0)
        self.assertEqual(reason, MinimizerNLCG.MAX_ITERATIONS_REACHED)

    def test_solution(self):
        self.minimizer.setTolerance(1e-4)
        self.minimizer.setMaxIterations(400)
        reason=self.minimizer.run(self.x0)
        x=self.minimizer.getResult()
        # We should be able to get a solution to set tolerance in #iterations
        self.assertEqual(reason, MinimizerNLCG.TOLERANCE_REACHED)
        self.assertAlmostEqual(np.amax(abs(x-self.xstar)), 0., places=3)

    def test_callback(self):
        n=[0]
        def callback(k, x, fg, gf):
            n[0]=n[0]+1
        self.minimizer.setCallback(callback)
        self.minimizer.setTolerance(1e-10)
        self.minimizer.setMaxIterations(10)
        reason=self.minimizer.run(self.x0)
        # callback should be called once for each iteration (including 0th)
        self.assertEqual(n[0], 11)


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMinimizerLBFGS))
    suite.addTest(unittest.makeSuite(TestMinimizerBFGS))
    suite.addTest(unittest.makeSuite(TestMinimizerNLCG))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

