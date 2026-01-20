
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
extra functions that can be used in symbolic expressions
"""

import sympy
from sympy.functions import *
from sympy import S

__author__="Cihan Altinay"

class wherePositive(sympy.Function):
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

class whereNonPositive(sympy.Function):
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

class whereNegative(sympy.Function):
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

class whereNonNegative(sympy.Function):
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

class whereZero(sympy.Function):
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

class whereNonZero(sympy.Function):
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

class log10(sympy.Function):
    """Returns the base-10 logarithm of the argument (same as log(x,10))
    """
    nargs = 1

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.exponential import log
        return log(arg,10)

class clip(sympy.Function):
    """Returns the argument clipped to a minimum and maximum value
    """
    nargs = (1,2,3)

class grad_n(sympy.Function):
    """Returns the spatial gradient of the argument
    """
    nargs = (2,3)

    def __str__(self):
        return "("+str(self.args[0])+"),"+str(self.args[1])

    @classmethod
    def eval(cls, *args):
        if args[0].is_zero: return S.Zero

class eigenvalues(sympy.Function):
    """Returns the Eigenvalues of the argument
    """
    pass

class eigenvalues_and_eigenvectors(sympy.Function):
    """Returns the Eigenvalues and Eigenvectors of the argument
    """
    pass

class minval(sympy.Function):
    """Returns the minimum value over all components of the argument
    """
    pass

class maxval(sympy.Function):
    """Returns the maximum value over all components of the argument
    """
    pass

class maximum(sympy.Function):
    """Returns the maximum over the arguments
    """
    pass

class minimum(sympy.Function):
    """Returns the minimum over the arguments
    """
    pass

class integrate(sympy.Function):
    """Returns the integral of the argument
    """
    pass

class interpolate(sympy.Function):
    """Returns the argument interpolated on the function space provided
    """
    pass

class L2(sympy.Function):
    """Returns the L2 norm of the argument
    """
    pass

class abs(sympy.Function):
    """Returns the absolute value of the argument
    """
    pass

#
# vim: expandtab shiftwidth=4:
