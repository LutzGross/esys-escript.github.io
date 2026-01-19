
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

from esys.escript.minimizer import MinimizerLBFGS, CostFunction, MinimizerException
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
import numpy as np

def createTest(p=2, dimension=3, numPoints=30):
    """
    this creates a simple test to find the point with the smallest distance from set of data points
    in the l^{p}-norm
    """
    import numpy as np
    from scipy.optimize import minimize_scalar

    theData = np.random.randint(-10, 10 + 1, size=(dimension, numPoints)).astype(float)
    def F0(x, *args):
        return np.sum(abs(x - theData[args[0], :]) ** p)

    x_true = np.zeros((dimension,))
    for i in range(dimension):
        R = minimize_scalar(F0, args=(i,), tol=1e-10)
        x_true[i] = R.x
    print("p = %g" % p)
    print("data = %s" % np.array_str(theData, precision=15, max_line_width=300))
    print("m_true = %s" % np.array_str(x_true, precision=15))
    return p, theData, x_true

class MinDist(CostFunction):
        def __init__(self, p=2, data=[], m_true=[]):
            super().__init__()

            self.p=p
            self.data=np.array(data)
            self.m_true = np.array(m_true)
            assert self.m_true.shape[0] == self.data.shape[0]
            self.dimension = self.m_true.shape[0]

        def getDualProduct(self, m, g):
            return np.dot(m, g)

        def getValue(self, m, *args):
            d,n = self.data.shape
            out=0
            for i in range(n):
                out+=np.sum(abs(m-self.data[:,i])**self.p)
            return out

        def getGradient(self, m, *args):
            d, n = self.data.shape
            out=np.zeros((d, ))
            for i in range(n):
                out+=self.p*abs(m-self.data[:,i])**(self.p-1)*np.sign(m-self.data[:,i])
            return out

        def getNorm(self, m):
            return np.linalg.norm(m)

        def getInverseHessianApproximation(self, r, m, *args, initializeHessian=False):
            if initializeHessian :
                d, n = self.data.shape
                self.h=np.zeros((d, ))
                for i in range(n):
                    self.h+=abs(m-self.data[:,i])**(self.p-2) #*self.p*(self.p-1)
            s = np.linalg.solve(np.diag(self.h), r)
            return s

class Test_MinimizerBFGS(unittest.TestCase):
    TEST1 = { 'p' : 2.5,
              'data': [[  5.,    3.,   -2.,   -8.,    0.,  -10.,   -3.,    8.,   -7.,    1.,   -7.,    4.,    1.,   10.,  -10.,   -1.,    3.,    6.,   -7.,    3.,    9.,    0.,    2.,   -8.,    4.,    7.,    1.,    4.,    5.,   -7., ],
                     [  2.,   -2.,    3.,   10.,    2.,   -3.,   -2.,   -9.,  -10.,   -8.,    4.,   10.,   -6.,    3.,   -3.,    9.,   -2.,   -6.,   -6.,    8.,    5.,    9.,   -9.,    7.,    6.,    3.,    6.,    5.,    8.,   -1., ],
                     [  9.,    8.,    5.,    9.,   -8.,    5.,    7.,    7.,   -2.,    4.,   -5.,   -5.,  -10.,    3.,   -7.,    9.,    2.,   -1.,   -3.,   -4.,    3.,   -3.,    9.,   10.,  -10.,   -5.,    0.,   -9.,    6.,    4., ]] ,
             'm_true' : [-0.101614684980357,  0.884081483782368,  0.759171570197073] }

    TEST2 = {'p' : 5,
            'data': [[-3., -10.,   7.,   3., -10., -5., 6., -1.,  10., -7.,  6., -1.,  6., 8., -8., 3., 3., -3., -2., -1., 10., -7., 6.,  5.,  1.,  0.,  8.,  5., -1., 8.],
                        [9., -7., -1., 3., 7., 7., -10., -2., -7., 5., -7., -3., -6., -2., -8.,  10., 4., 2., -10., -1., 10., 2., -7., -7., -10., -5., 0., -4., 3.,  10.],
                        [-4.,  5., -3., -1., 7., -4., 1., -4., -3., 8., -7., -3., -5., 0., -7., -1., 1., -5., 4., 10., 0., -4., 7., 3., -5., 4., -6., -1., -8., 5.]],
            'm_true':  [0.349084852003848, -0.144257203323146,  0.388851745831033] }

    TEST3 = { 'p':  9,
            'data': [[7., 2., -8., 9.,  2., 4., -4., 10., 3., 1., -4., 8., -5., -9., -8., 6., 10., -7., -7., 4., -3., -5., 1.,   10.,    1.,  -7.,    6.,    7.,    5.,  -5.],
                    [5.,  - 10.,  - 5.,    3.,  - 10.,    4.,  - 3.,  - 6.,  - 8.,  - 7.,    8.,    4.,  - 3.,    7.,  - 8.,  - 3.,  - 7.,    5.,   10.,    1.,  - 1., 9.,    5.,  - 6.,    8.,  - 3.,    6.,  - 4.,  - 9.,    8.],
                    [-3., 5., 10., 9., 9., 0., 5.,  - 4.,  - 3., 5.,  - 5.,  - 1., 0.,  - 6.,  - 7., 3.,  - 7.,  - 3., 10., 5., 7.,  - 8., 0., 7.,  - 5., 2.,  - 6., 10., 0., 7.]],
            'm_true': [0.756394642236015, -0.210685240100182,  1.315069006925737] }
    def testTEST1(self):
        F = MinDist(**self.TEST1)
        solve = MinimizerLBFGS(F, m_tol=1e-7, grad_tol=None, iterMax=100)
        m0 = np.zeros((F.dimension,))
        solve.run(m0)
        m = solve.getResult()
        self.assertTrue(np.allclose(m, F.m_true))

    def testTEST2(self):
        F = MinDist(**self.TEST2)
        solve = MinimizerLBFGS(F, m_tol=1e-7, grad_tol=None,  iterMax=100)
        m0 = np.zeros((F.dimension,))
        solve.run(m0)
        m = solve.getResult()
        self.assertTrue(np.allclose(m, F.m_true))

    def testTEST3(self):
        F = MinDist(**self.TEST3)
        solve = MinimizerLBFGS(F, m_tol=1e-7, grad_tol=None,  iterMax=100)
        m0 = np.zeros((F.dimension,))
        solve.run(m0)
        m = solve.getResult()
        self.assertTrue(np.allclose(m, F.m_true))

    def testTEST4(self):
        F = MinDist(**self.TEST1)
        solve = MinimizerLBFGS(F, m_tol=1e-10, grad_tol=1e-8, iterMax=100)
        m0 = np.zeros((F.dimension,))
        solve.run(m0)
        m = solve.getResult()
        #print((m-F.m_true)/F.m_true)
        self.assertTrue(np.allclose(m, F.m_true, rtol=1e-3))

    def testOptions(self):
        """
        test for setting options
        """
        F = MinDist(**self.TEST1)
        solve = MinimizerLBFGS(F, m_tol=1e-2, iterMax=23)
        self.assertEqual(solve.getOptions()['m_tol'], 1e-2)
        self.assertEqual(solve.getOptions()['iterMax'], 23)

        solve.setTolerance(m_tol=2., grad_tol=1.)
        solve.setIterMax(iterMax=1234)
        self.assertEqual(solve.getOptions()['m_tol'], 2.)
        self.assertEqual(solve.getOptions()['grad_tol'], 1.)
        self.assertEqual(solve.getOptions()['iterMax'], 1234)
        # test options:
        self.assertRaises(KeyError, solve.setOptions, BOB=-1)

        OPTS={'m_tol' :1,
        'grad_tol' : 2,
        'iterMax' : 3,
        'truncation' : 4,
        'restart' : 5,
        'relAlphaMin' : 6,
        'scaleSearchDirection' : 7,
        'initialAlpha' : 8}

        solve.setOptions(**OPTS)
        newOPTS=solve.getOptions()
        self.assertTrue('line_search' in newOPTS.keys())
        del newOPTS['line_search']
        self.assertDictEqual(newOPTS, OPTS)

        # linesearch:
        self.assertRaises(KeyError, solve.getLineSearch().setOptions, BOB=-1)
        lsOPTS={
        'phiEpsilon': 1,
        'alphaMin': 3,
        'alphaMax': 4,
        'overStepFactor': 5,
        'iterMax': 6,
        'c1': 7,
        'c2': 8,
        'inter_order': 2,
        'inter_tol': 9,
        'inter_iterMax': 10,
        'alphaOffset': 11,
        'alphaWidthMin': 12,
        'zoom_iterMax': 13,
        'zoom_reductionMin': 0.5
        }
        solve.getLineSearch().setOptions(**lsOPTS)
        newOPTS=solve.getOptions()
        self.assertDictEqual(newOPTS['line_search'], lsOPTS)


        # test of call back is called:
        def CB(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
            raise ValueError
        self.assertRaises(TypeError, solve.setCallback, 5)
        solve.setCallback(CB)
        self.assertRaises(ValueError, solve.run, np.zeros((F.dimension,)))

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

