
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for the escript.symbolic module

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Cihan Altinay"

from esys.escript import Data, RandomData, FunctionSpace
from esys.escript.symbolic import *
import unittest

class SymbolicTestCase(unittest.TestCase):

    def test_Evaluator(self):
        e=Evaluator()
        self.assertEqual(len(e), 0, "empty evaluator returns wrong length")
        self.assertEqual(e.evaluate(), (), "result of evaluate() not empty")

        x=Symbol('x')
        e=Evaluator(x*x, x**3)
        self.assertEqual(len(e), 2, "evaluator returns wrong length")
        self.assertEqual(e[0], x*x, "first expression wrong")
        self.assertEqual(e[1], x**3, "second expression wrong")

        f=e.addExpression(x**4)
        self.assertEqual(len(e), 3, "wrong length after addExpression()")
        self.assertEqual(e, f, "addExpression() did not return self")
        self.assertEqual(e[2], x**4, "third expression wrong")

        e+=x**5
        self.assertEqual(len(e), 4, "wrong length after += operator")
        self.assertEqual(e[3], x**5, "fourth expression wrong")

        self.assertRaises(RuntimeError, e.evaluate)
        f=e.subs(x=2)
        self.assertEqual(e, f, "subs() did not return self")
        self.assertEqual(e.evaluate(), (4,8,16,32), "wrong result after subs()")
        self.assertEqual(e(x=3), (9,27,81,243), "wrong result after __call__")

        xx=RandomData((4,), FunctionSpace())
        ref=[d.toListOfTuples() for d in (xx**2, xx**3, xx**4, xx**5)]
        res=e(x=xx)
        for d in res:
            self.assertTrue(isinstance(d, Data), "substituted expression not a Data object")
        res=[x.toListOfTuples() for x in res]
        self.assertEqual(res, ref, "wrong result after substitution with Data object")

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def test_Functions(self):
        import esys.escript as es
        x=Symbol('x')
        xx=RandomData((4,4), FunctionSpace())
        TOL_DIGITS=8

        y=wherePositive(x)
        ref=es.wherePositive(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for wherePositive()")

        y=whereNonPositive(x)
        ref=es.whereNonPositive(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for whereNonPositive()")

        y=whereNegative(x)
        ref=es.whereNegative(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for whereNegative()")

        y=whereNonNegative(x)
        ref=es.whereNonNegative(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for whereNonNegative()")

        y=whereZero(x)
        ref=es.whereZero(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for whereZero()")

        y=whereNonZero(x)
        ref=es.whereNonZero(xx).toListOfTuples()
        res=Evaluator(y)(x=xx)[0].toListOfTuples()
        self.assertEqual(res, ref, "wrong result for whereNonZero()")

        y=log10(x)
        ref=es.log10(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for log10()")

        y=inverse(x)
        ref=es.inverse(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for inverse()")

        y=minval(x)
        ref=es.minval(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for minval()")

        y=maxval(x)
        ref=es.maxval(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for maxval()")

        y=trace(x)
        ref=es.trace(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for trace()")

        y=transpose(x)
        ref=es.transpose(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for transpose()")

        y=symmetric(x)
        ref=es.symmetric(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for symmetric()")

        y=nonsymmetric(x)
        ref=es.nonsymmetric(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for nonsymmetric()")

        y=swap_axes(x)
        ref=es.swap_axes(xx)
        res=Evaluator(y)(x=xx)[0]
        self.assertAlmostEqual(es.Lsup(res-ref), 0.0, TOL_DIGITS, "wrong result for swap_axes()")

        #y=grad(x)
        #ref=es.grad(xx).toListOfTuples()
        #res=Evaluator(y)(x=xx)[0].toListOfTuples()
        #self.assertEqual(res, ref, "wrong result for grad()")


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SymbolicTestCase))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

