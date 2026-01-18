# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""
Elementary mathematical functions and conditional masks for escript.

This module provides mathematical functions (sin, cos, exp, sqrt, etc.) and
conditional mask functions (wherePositive, whereNegative, etc.) that work
with escript Data objects, numpy arrays, and symbolic expressions.
"""

from __future__ import print_function, division

__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://launchpad.net/escript-finley"

import math
import cmath
import numpy

from . import escriptcpp as escore
from .start import HAVE_SYMBOLS

if HAVE_SYMBOLS:
    from . import symboliccore as sym
else:
    class sym:
        Symbol = type(memoryview(b'1'))  # This type should never be passed

# Machine epsilon for tolerance calculations
EPSILON = escore.getMachinePrecision()

__all__ = [
    'log10', 'wherePositive', 'whereNegative', 'whereNonNegative',
    'whereNonPositive', 'whereZero', 'whereNonZero', 'Abs', 'erf',
    'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'atan2',
    'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh',
    'exp', 'sqrt', 'log', 'sign', 'minval', 'maxval'
]


def log10(arg):
    """
    Returns base-10 logarithm of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.log10(arg)
    elif isinstance(arg, escore.Data):
        return arg._log10()
    elif isinstance(arg, complex):
        return cmath.log10(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.log10(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.log10)
    else:
        raise TypeError("log10: Unknown argument type.")


def wherePositive(arg):
    """
    Returns mask of positive values of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("wherePositive: operation not supported for complex")
        out = numpy.greater(arg, numpy.zeros(arg.shape, numpy.double)) * 1.
        return out
    elif isinstance(arg, escore.Data):
        return arg._wherePositive()
    elif isinstance(arg, float):
        if arg > 0:
            return 1.
        else:
            return 0.
    elif isinstance(arg, int):
        if arg > 0:
            return 1.
        else:
            return 0.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.wherePositive)
    elif isinstance(arg, complex):
        raise TypeError("wherePositive: operation not supported for complex")
    else:
        raise TypeError("wherePositive: Unknown argument type.")


def whereNegative(arg):
    """
    Returns mask of negative values of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("whereNegative: operation not supported for complex")
        out = numpy.less(arg, numpy.zeros(arg.shape, numpy.double)) * 1.
        return out
    elif isinstance(arg, escore.Data):
        return arg._whereNegative()
    elif isinstance(arg, float):
        if arg < 0:
            return 1.
        else:
            return 0.
    elif isinstance(arg, int):
        if arg < 0:
            return 1.
        else:
            return 0.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.whereNegative)
    elif isinstance(arg, complex):
        raise TypeError("whereNegative: operation not supported for complex")
    else:
        raise TypeError("whereNegative: Unknown argument type.")


def whereNonNegative(arg):
    """
    Returns mask of non-negative values of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("whereNonNegative: operation not supported for complex")
        out = numpy.greater_equal(arg, numpy.zeros(arg.shape, numpy.double)) * 1.
        return out
    elif isinstance(arg, escore.Data):
        return arg._whereNonNegative()
    elif isinstance(arg, float):
        if arg < 0:
            return 0.
        else:
            return 1.
    elif isinstance(arg, int):
        if arg < 0:
            return 0.
        else:
            return 1.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.whereNonNegative)
    elif isinstance(arg, complex):
        raise TypeError("whereNonNegative: operation not supported for complex")
    else:
        raise TypeError("whereNonNegative: Unknown argument type.")


def whereNonPositive(arg):
    """
    Returns mask of non-positive values of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("whereNonPositive: operation not supported for complex")
        out = numpy.less_equal(arg, numpy.zeros(arg.shape, numpy.double)) * 1.
        return out
    elif isinstance(arg, escore.Data):
        return arg._whereNonPositive()
    elif isinstance(arg, float):
        if arg > 0:
            return 0.
        else:
            return 1.
    elif isinstance(arg, int):
        if arg > 0:
            return 0.
        else:
            return 1.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.whereNonPositive)
    elif isinstance(arg, complex):
        raise TypeError("whereNonPositive: operation not supported for complex")
    else:
        raise TypeError("whereNonPositive: Unknown argument type.")


def whereZero(arg, tol=None, rtol=math.sqrt(EPSILON)):
    """
    Returns mask of zero entries of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :param tol: absolute tolerance. Values with absolute value less than tol are accepted
                as zero. If ``tol`` is not present ``rtol * Lsup(arg)`` is used.
    :type tol: ``float``
    :param rtol: relative tolerance used to define the absolute tolerance if ``tol`` is not present.
    :type rtol: non-negative ``float``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise ValueError: if ``rtol`` is negative.
    :raise TypeError: if the type of the argument is not expected
    """
    if tol is None and not isinstance(arg, sym.Symbol):
        if rtol < 0:
            raise ValueError("rtol must be non-negative.")
        # Import Lsup locally to avoid circular import
        from . import util
        tol = util.Lsup(arg) * rtol
    if isinstance(arg, numpy.ndarray):
        out = numpy.less_equal(abs(arg) - tol, numpy.zeros(arg.shape, numpy.double)) * 1.
        if isinstance(out, float):
            out = numpy.array(out, dtype=numpy.double)
        return out
    elif isinstance(arg, escore.Data):
        return arg._whereZero(tol)
    elif isinstance(arg, float) or isinstance(arg, complex) or isinstance(arg, int):
        if abs(arg) <= tol:
            return 1.
        else:
            return 0.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.whereZero)
    else:
        raise TypeError("whereZero: Unknown argument type.")


def whereNonZero(arg, tol=0.):
    """
    Returns mask of values different from zero of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :param tol: absolute tolerance. Values with absolute value less than or equal
                to tol are considered zero.
    :type tol: ``float``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        out = numpy.greater(abs(arg) - tol, numpy.zeros(arg.shape, numpy.double)) * 1.
        if isinstance(out, float):
            out = numpy.array(out, dtype=numpy.double)
        return out
    elif isinstance(arg, escore.Data):
        return arg._whereNonZero(tol)
    elif isinstance(arg, float) or isinstance(arg, complex) or isinstance(arg, int):
        if abs(arg) > tol:
            return 1.
        else:
            return 0.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.whereNonZero)
    else:
        raise TypeError("whereNonZero: Unknown argument type.")


def Abs(arg):
    """
    Returns the absolute value of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.abs)
    else:
        return abs(arg)


def erf(arg):
    """
    Returns the error function *erf* of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, escore.Data):
        return arg._erf()
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.erf)
    else:
        raise TypeError("erf: Unknown argument type.")


def sin(arg):
    """
    Returns sine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.sin(arg)
    elif isinstance(arg, escore.Data):
        return arg._sin()
    elif isinstance(arg, complex):
        return cmath.sin(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.sin(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.sin)
    else:
        raise TypeError("sin: Unknown argument type.")


def cos(arg):
    """
    Returns cosine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.cos(arg)
    elif isinstance(arg, escore.Data):
        return arg._cos()
    elif isinstance(arg, complex):
        return cmath.cos(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.cos(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.cos)
    else:
        raise TypeError("cos: Unknown argument type.")


def tan(arg):
    """
    Returns tangent of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.tan(arg)
    elif isinstance(arg, escore.Data):
        return arg._tan()
    elif isinstance(arg, complex):
        return cmath.tan(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.tan(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.tan)
    else:
        raise TypeError("tan: Unknown argument type.")


def asin(arg):
    """
    Returns the inverse sine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arcsin(arg)
    elif isinstance(arg, escore.Data):
        return arg._asin()
    elif isinstance(arg, complex):
        return cmath.asin(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.asin(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.asin)
    else:
        raise TypeError("asin: Unknown argument type.")


def acos(arg):
    """
    Returns the inverse cosine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arccos(arg)
    elif isinstance(arg, escore.Data):
        return arg._acos()
    elif isinstance(arg, complex):
        return cmath.acos(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.acos(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.acos)
    else:
        raise TypeError("acos: Unknown argument type.")


def atan(arg):
    """
    Returns inverse tangent of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arctan(arg)
    elif isinstance(arg, escore.Data):
        return arg._atan()
    elif isinstance(arg, complex):
        return cmath.atan(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.atan(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.atan)
    else:
        raise TypeError("atan: Unknown argument type.")


def atan2(arg0, arg1):
    """
    Returns inverse tangent of ``arg0/arg1``, handling quadrant correctly.

    :param arg0: numerator argument (y coordinate)
    :type arg0: ``float``, `escript.Data`, ``numpy.ndarray``
    :param arg1: denominator argument (x coordinate)
    :type arg1: ``float``, `escript.Data`, ``numpy.ndarray``
    :return: angle in radians between -pi and pi
    :rtype: ``float``, `escript.Data`, ``numpy.ndarray`` depending on input types
    """
    m = whereZero(arg1, rtol=EPSILON)
    m2 = whereNegative(arg1 * arg0)
    s = atan(arg0 / (arg1 + m)) * (1 - m) + (numpy.pi / 2) * (1 - 2 * m2) * m
    s += (wherePositive(arg0) * whereNegative(arg1) - whereNegative(arg0) * whereNegative(arg1)) * numpy.pi
    return s


def sinh(arg):
    """
    Returns the hyperbolic sine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.sinh(arg)
    elif isinstance(arg, escore.Data):
        return arg._sinh()
    elif isinstance(arg, complex):
        return cmath.sinh(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.sinh(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.sinh)
    else:
        raise TypeError("sinh: Unknown argument type.")


def cosh(arg):
    """
    Returns the hyperbolic cosine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.cosh(arg)
    elif isinstance(arg, escore.Data):
        return arg._cosh()
    elif isinstance(arg, complex):
        return cmath.cosh(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.cosh(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.cosh)
    else:
        raise TypeError("cosh: Unknown argument type.")


def tanh(arg):
    """
    Returns the hyperbolic tangent of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.tanh(arg)
    elif isinstance(arg, escore.Data):
        return arg._tanh()
    elif isinstance(arg, complex):
        return cmath.tanh(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.tanh(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.tanh)
    else:
        raise TypeError("tanh: Unknown argument type.")


def asinh(arg):
    """
    Returns the inverse hyperbolic sine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arcsinh(arg)
    elif isinstance(arg, escore.Data):
        return arg._asinh()
    elif isinstance(arg, complex):
        return numpy.arcsinh(complex(arg))
    elif isinstance(arg, float) or isinstance(arg, int):
        return numpy.arcsinh(float(arg))
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.asinh)
    else:
        raise TypeError("asinh: Unknown argument type.")


def acosh(arg):
    """
    Returns the inverse hyperbolic cosine of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arccosh(arg)
    elif isinstance(arg, escore.Data):
        return arg._acosh()
    elif isinstance(arg, complex):
        return numpy.arccosh(complex(arg))
    elif isinstance(arg, float) or isinstance(arg, int):
        return numpy.arccosh(float(arg))
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.acosh)
    else:
        raise TypeError("acosh: Unknown argument type.")


def atanh(arg):
    """
    Returns the inverse hyperbolic tangent of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.arctanh(arg)
    elif isinstance(arg, escore.Data):
        return arg._atanh()
    elif isinstance(arg, complex):
        return numpy.arctanh(complex(arg))
    elif isinstance(arg, float) or isinstance(arg, int):
        return numpy.arctanh(float(arg))
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.atanh)
    else:
        raise TypeError("atanh: Unknown argument type.")


def exp(arg):
    """
    Returns *e* to the power of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of arg
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.exp(arg)
    elif isinstance(arg, escore.Data):
        return arg._exp()
    elif isinstance(arg, complex):
        return cmath.exp(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.exp(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.exp)
    else:
        raise TypeError("exp: Unknown argument type.")


def sqrt(arg):
    """
    Returns the square root of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
            depending on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.sqrt(arg)
    elif isinstance(arg, escore.Data):
        return arg._sqrt()
    elif isinstance(arg, complex):
        return cmath.sqrt(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.sqrt(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.sqrt)
    else:
        raise TypeError("sqrt: Unknown argument type.")


def log(arg):
    """
    Returns the natural logarithm of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        return numpy.log(arg)
    elif isinstance(arg, escore.Data):
        return arg._log()
    elif isinstance(arg, complex):
        return cmath.log(arg)
    elif isinstance(arg, float) or isinstance(arg, int):
        return math.log(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.log)
    else:
        raise TypeError("log: Unknown argument type.")


def sign(arg):
    """
    Returns the sign of argument ``arg``.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
            on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("sign: operation not supported for complex")
        return wherePositive(arg) - whereNegative(arg)
    elif isinstance(arg, escore.Data):
        return arg._sign()
    elif isinstance(arg, complex):
        raise TypeError("sign: operation not supported for complex")
    elif isinstance(arg, float) or isinstance(arg, int):
        if arg > 0:
            return 1.
        elif arg < 0:
            return -1.
        else:
            return 0.
    elif isinstance(arg, sym.Symbol):
        return arg.applyfunc(sym.symfn.sign)
    else:
        raise TypeError("sign: Unknown argument type.")


def minval(arg):
    """
    Returns the minimum value over all components of ``arg`` at each data point.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("minval: operation not supported for complex")
        if arg.ndim == 0:
            return float(arg)
        else:
            return arg.min()
    elif isinstance(arg, escore.Data):
        return arg._minval()
    elif isinstance(arg, float):
        return arg
    elif isinstance(arg, int):
        return float(arg)
    elif isinstance(arg, sym.Symbol):
        return sym.symfn.minval(arg)
    else:
        raise TypeError("minval: Unknown argument type.")


def maxval(arg):
    """
    Returns the maximum value over all components of ``arg`` at each data point.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
    :raise TypeError: if the type of the argument is not expected
    """
    if isinstance(arg, numpy.ndarray):
        if arg.dtype.kind == 'c':
            raise TypeError("maxval: operation not supported for complex")
        if arg.ndim == 0:
            return float(arg)
        else:
            return arg.max()
    elif isinstance(arg, escore.Data):
        return arg._maxval()
    elif isinstance(arg, float):
        return arg
    elif isinstance(arg, int):
        return float(arg)
    elif isinstance(arg, sym.Symbol):
        return sym.symfn.maxval(arg)
    else:
        raise TypeError("maxval: Unknown argument type.")
