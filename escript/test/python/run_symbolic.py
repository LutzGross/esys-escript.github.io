
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

from esys.escript import *
import unittest

class SymbolicTestCase(unittest.TestCase):

    # number of digits that have to match for results to be considered equal
    TOL_DIGITS=8

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
    def test_wherePositive_Symbol(self):
        x=Symbol('x')
        y=wherePositive(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=wherePositive(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_whereNonPositive_Symbol(self):
        x=Symbol('x')
        y=whereNonPositive(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=whereNonPositive(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_whereNegative_Symbol(self):
        x=Symbol('x')
        y=whereNegative(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=whereNegative(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_whereNonNegative_Symbol(self):
        x=Symbol('x')
        y=whereNonNegative(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=whereNonNegative(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_whereZero_Symbol(self):
        x=Symbol('x')
        y=whereZero(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=whereZero(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_whereNonZero_Symbol(self):
        x=Symbol('x')
        y=whereNonZero(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=whereNonZero(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_log10_Symbol(self):
        x=Symbol('x')
        y=log10(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=log10(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_inverse_Symbol(self):
        x=Symbol('x')
        y=inverse(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=inverse(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_minval_Symbol(self):
        x=Symbol('x')
        y=minval(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=minval(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_maxval_Symbol(self):
        x=Symbol('x')
        y=maxval(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=maxval(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_trace_Symbol(self):
        x=Symbol('x')
        y=trace(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=trace(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_transpose_Symbol(self):
        x=Symbol('x')
        y=transpose(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=transpose(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_symmetric_Symbol(self):
        x=Symbol('x')
        y=symmetric(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=symmetric(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_nonsymmetric_Symbol(self):
        x=Symbol('x')
        y=nonsymmetric(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=nonsymmetric(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_swapaxes_Symbol(self):
        x=Symbol('x')
        y=swap_axes(x, 1, 0)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=swap_axes(xx, 1, 0)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_sin_Symbol(self):
        x=Symbol('x')
        y=sin(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=sin(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_cos_Symbol(self):
        x=Symbol('x')
        y=cos(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=cos(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_tan_Symbol(self):
        x=Symbol('x')
        y=tan(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=tan(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_asin_Symbol(self):
        x=Symbol('x')
        y=asin(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=asin(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_acos_Symbol(self):
        x=Symbol('x')
        y=acos(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=acos(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_atan_Symbol(self):
        x=Symbol('x')
        y=atan(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=atan(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_sinh_Symbol(self):
        x=Symbol('x')
        y=sinh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=sinh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_cosh_Symbol(self):
        x=Symbol('x')
        y=cosh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=cosh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_tanh_Symbol(self):
        x=Symbol('x')
        y=tanh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=tanh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_asinh_Symbol(self):
        x=Symbol('x')
        y=asinh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=asinh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_acosh_Symbol(self):
        x=Symbol('x')
        y=acosh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=acosh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_atanh_Symbol(self):
        x=Symbol('x')
        y=atanh(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=atanh(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_exp_Symbol(self):
        x=Symbol('x')
        y=exp(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=exp(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_sqrt_Symbol(self):
        x=Symbol('x')
        y=sqrt(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=sqrt(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_log_Symbol(self):
        x=Symbol('x')
        y=log(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=log(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_sign_Symbol(self):
        x=Symbol('x')
        y=sign(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=sign(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_abs_Symbol(self):
        x=Symbol('x')
        y=abs(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=abs(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_inner_Symbol(self):
        x,y=symbols('xy')
        z=inner(x,y)
        self.assertTrue(isinstance(z, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        yy=RandomData((4,4), FunctionSpace())
        ref=inner(xx,yy)
        res=Evaluator(z)(x=xx, y=yy)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_outer_Symbol(self):
        x,y=symbols('xy')
        z=outer(x,y)
        self.assertTrue(isinstance(z, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        yy=RandomData((4,4), FunctionSpace())
        ref=outer(xx,yy)
        res=Evaluator(z)(x=xx, y=yy)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_clip_Symbol(self):
        x=Symbol('x')
        y=clip(x, 0.4, 0.6)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=clip(xx, 0.4, 0.6)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_length_Symbol(self):
        x=Symbol('x')
        y=length(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        ref=length(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_maximum_Symbol(self):
        x,y=symbols('xy')
        z=maximum(x,y)
        self.assertTrue(isinstance(z, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        yy=RandomData((4,4), FunctionSpace())
        ref=maximum(xx,yy)
        res=Evaluator(z)(x=xx, y=yy)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_minimum_Symbol(self):
        x,y=symbols('xy')
        z=minimum(x,y)
        self.assertTrue(isinstance(z, Basic), "wrong type of result")
        xx=RandomData((4,4), FunctionSpace())
        yy=RandomData((4,4), FunctionSpace())
        ref=minimum(xx,yy)
        res=Evaluator(z)(x=xx, y=yy)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, self.TOL_DIGITS, "wrong result")


if __name__ == "__main__":
    import sys
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SymbolicTestCase))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

