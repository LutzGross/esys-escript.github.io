
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
extra functions that can be used in symbolic expressions
"""

from sympy.core.function import Function
from sympy import S

class wherePositive(Function):
    """Returns:
         1 where expr > 0
         0 else
    """
    nargs = 1
    is_bounded = True
    is_negative = False
    is_real = True

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg.is_positive: return S.One
        if arg.is_negative or arg is S.Zero: return S.Zero
        if arg.is_Function:
            if arg.func is wherePositive: return arg
            if arg.func is whereNegative: return arg

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

class whereNonPositive(Function):
    """Returns:
         0 where expr > 0
         1 else
    """
    nargs = 1
    is_bounded = True
    is_negative = False
    is_real = True

    @classmethod
    def eval(cls, arg):
        return 1-wherePositive(arg)

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

class whereNegative(Function):
    """Returns:
         1 where expr < 0
         0 else
    """
    nargs = 1
    is_bounded = True
    is_negative = False
    is_real = True

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg.is_nonnegative: return S.Zero
        if arg.is_negative: return S.One
        if arg.is_Function:
            if arg.func is wherePositive: return S.Zero

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

class whereNonNegative(Function):
    """Returns:
         0 where expr < 0
         1 else
    """
    nargs = 1
    is_bounded = True
    is_negative = False
    is_real = True

    @classmethod
    def eval(cls, arg):
        return 1-whereNegative(arg)

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

#
# vim: expandtab shiftwidth=4:
