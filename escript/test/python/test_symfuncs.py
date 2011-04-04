
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
Test suite for the symbols.py module

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

class Test_symfuncs(unittest.TestCase):
    def test_Symbolic_grad(self):
        x=Symbol('x')
        y=grad(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=self.domain.getX()
        ref=grad(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, int(-log10(self.RES_TOL)), "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_Symbolic_integrate(self):
        x=Symbol('x')
        fs=Symbol('fs')
        y=integrate(x,fs)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=self.domain.getX()
        _fs=Function(self.domain)
        ref=integrate(xx, _fs)
        res=Evaluator(y)(x=xx,fs=_fs)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, int(-log10(self.RES_TOL)), "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_Symbolic_interpolate(self):
        x=Symbol('x')
        fs=Symbol('fs')
        y=interpolate(x, fs)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=self.domain.getX()
        _fs=Function(self.domain)
        ref=interpolate(xx, _fs)
        res=Evaluator(y)(x=xx, fs=_fs)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, int(-log10(self.RES_TOL)), "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_Symbolic_div(self):
        x=Symbol('x')
        fs=Symbol('fs')
        y=div(x, fs)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=self.domain.getX()
        _fs=Function(self.domain)
        ref=div(xx, _fs)
        res=Evaluator(y)(x=xx, fs=_fs)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, int(-log10(self.RES_TOL)), "wrong result")
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def test_Symbolic_L2(self):
        x=Symbol('x')
        y=L2(x)
        self.assertTrue(isinstance(y, Basic), "wrong type of result")
        xx=self.domain.getX()
        ref=L2(xx)
        res=Evaluator(y)(x=xx)
        self.assertAlmostEqual(Lsup(res-ref), 0.0, int(-log10(self.RES_TOL)), "wrong result")

