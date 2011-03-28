
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

class whereZero(Function):
    """Returns:
         1 where expr == 0
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
        if arg.is_zero: return S.One
        if arg.is_nonzero: return S.Zero
        if arg.is_Function:
            if arg.func is whereZero: return 1-arg

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

class whereNonZero(Function):
    """Returns:
         0 where expr == 0
         1 else
    """
    nargs = 1
    is_bounded = True
    is_negative = False
    is_real = True

    @classmethod
    def eval(cls, arg):
        return 1-whereZero(arg)

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, x):
        return S.Zero

class log10(Function):
    """Returns the base-10 logarithm of the argument (same as log(x,10))
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.exponential import log
        return log(arg,10)

class inverse(Function):
    """Returns the inverse of the argument
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            return 1./arg

class minval(Function):
    """Returns the minimum value over all components of the argument
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            return arg

class maxval(Function):
    """Returns the maximum value over all components of the argument
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            return arg

class trace(Function):
    """Returns the trace of the argument with optional axis_offset
    """
    nargs = (1,2)

class transpose(Function):
    """Returns the transpose of the argument
    """
    nargs = (1,2)

class symmetric(Function):
    """Returns the symmetric part of the argument
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Function:
            if arg.func is symmetric: return arg
            if arg.func is nonsymmetric: return S.Zero
        elif arg.is_Number:
            return arg

class nonsymmetric(Function):
    """Returns the non-symmetric part of the argument
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Function:
            if arg.func is nonsymmetric: return arg
            if arg.func is symmetric: return S.Zero
        elif arg.is_Number:
            return arg

class swap_axes(Function):
    """Returns the 'swap' of the argument
    """
    nargs = (1,2,3)

class grad(Function):
    """Returns the spatial gradient of the argument
    """
    nargs = (1,2)

#
# vim: expandtab shiftwidth=4:
