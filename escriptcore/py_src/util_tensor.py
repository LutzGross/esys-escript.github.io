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
Tensor operations for escript.

This module provides tensor construction, algebra, arithmetic, and product
operations for escript Data objects, numpy arrays, and symbolic expressions.
"""


__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://launchpad.net/escript-finley"

import numpy

from . import escriptcpp as escore
from .escriptcpp import C_GeneralTensorProduct
from .start import HAVE_SYMBOLS
from .util_math import sqrt, wherePositive, whereNegative, whereNonPositive, whereNonNegative

if HAVE_SYMBOLS:
    from . import symboliccore as sym
else:
    class sym:
        Symbol = type(memoryview(b'1'))

__all__ = [
    # Construction
    'kronecker', 'identity', 'zeros', 'identityTensor', 'identityTensor4', 'unitVector',
    # Algebra
    'length', 'trace', 'transpose', 'swap_axes', 'symmetric', 'nonsymmetric',
    'antisymmetric', 'hermitian', 'antihermitian', 'inverse', 'escript_inverse',
    'eigenvalues', 'eigenvalues_and_eigenvectors',
    # Arithmetic
    'mult', 'maximum', 'minimum', 'clip', 'cross',
    # Products
    'inner', 'outer', 'matrixmult', 'matrix_mult', 'tensormult', 'tensor_mult',
    'generalTensorProduct', 'escript_generalTensorProduct',
    'transposed_matrix_mult', 'transposed_tensor_mult',
    'generalTransposedTensorProduct', 'escript_generalTransposedTensorProduct',
    'matrix_transposed_mult', 'tensor_transposed_mult',
    'generalTensorTransposedProduct', 'escript_generalTensorTransposedProduct'
]


# =========================================================================
#   Tensor Construction
# =========================================================================

def kronecker(d=3):
    """
    Returns the kronecker delta-symbol.

    :param d: dimension or an object that has the ``getDim`` method defining the
              dimension
    :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
    :return: the object u of rank 2 with *u[i,j]=1* for *i=j* and *u[i,j]=0*
             otherwise
    :rtype: ``numpy.ndarray`` or `escript.Data` of rank 2
    """
    return identityTensor(d)


def identity(shape=()):
    """
    Returns the ``shape`` x ``shape`` identity tensor.

    :param shape: input shape for the identity tensor
    :type shape: ``tuple`` of ``int``
    :return: array whose shape is shape x shape where *u[i,k]=1* for *i=k* and
             *u[i,k]=0* otherwise for len(shape)=1. If len(shape)=2:
             *u[i,j,k,l]=1* for *i=k and j=l* and *u[i,j,k,l]=0* otherwise.
    :rtype: ``numpy.ndarray`` of rank 1, rank 2 or rank 4
    :raise ValueError: if len(shape)>2
    """
    if len(shape) > 0:
        out = numpy.zeros(shape + shape, numpy.float64)
        if len(shape) == 1:
            for i0 in range(shape[0]):
                out[i0, i0] = 1.
        elif len(shape) == 2:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    out[i0, i1, i0, i1] = 1.
        else:
            raise ValueError("identity: length of shape is restricted to 2.")
    else:
        out = 1.
    return out


def zeros(shape=()):
    """
    Returns the ``shape`` zero tensor.

    :param shape: input shape for the zero tensor
    :type shape: ``tuple`` of ``int``
    :return: array of shape filled with zeros
    :rtype: ``numpy.ndarray``
    """
    if len(shape) > 0:
        out = numpy.zeros(shape, numpy.float64)
    else:
        out = 0.
    return out


def identityTensor(d=3):
    """
    Returns the ``d`` x ``d`` identity matrix.

    :param d: dimension or an object that has the ``getDim`` method defining the
              dimension
    :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
    :return: the object u of rank 2 with *u[i,j]=1* for *i=j* and *u[i,j]=0*
             otherwise
    :rtype: ``numpy.ndarray`` or `escript.Data` of rank 2
    """
    if isinstance(d, escore.FunctionSpace):
        return escore.Data(identity((d.getDim(),)), d)
    elif isinstance(d, escore.Domain):
        return identity((d.getDim(),))
    else:
        return identity((d,))


def identityTensor4(d=3):
    """
    Returns the ``d`` x ``d`` x ``d`` x ``d`` identity tensor.

    :param d: dimension or an object that has the ``getDim`` method defining the
              dimension
    :type d: ``int`` or any object with a ``getDim`` method
    :return: the object u of rank 4 with *u[i,j,k,l]=1* for *i=k and j=l* and
             *u[i,j,k,l]=0* otherwise
    :rtype: ``numpy.ndarray`` or `escript.Data` of rank 4
    """
    if isinstance(d, escore.FunctionSpace):
        return escore.Data(identity((d.getDim(), d.getDim())), d)
    elif isinstance(d, escore.Domain):
        return identity((d.getDim(), d.getDim()))
    else:
        return identity((d, d))


def unitVector(i=0, d=3):
    """
    Returns a unit vector u of dimension d whose non-zero element is at index i.

    :param i: index for non-zero element
    :type i: ``int``
    :param d: dimension or an object that has the ``getDim`` method defining the
              dimension
    :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
    :return: the object u of rank 1 with *u[j]=1* for *j=index* and *u[j]=0*
             otherwise
    :rtype: ``numpy.ndarray`` or `escript.Data` of rank 1
    """
    return kronecker(d)[i]


# =========================================================================
#   Tensor Algebra
# =========================================================================

def length(arg):
    """
    Returns the length (Euclidean norm) of argument ``arg`` at each data point.

    :param arg: argument
    :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
    :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
    """
    a = abs(arg)
    return sqrt(inner(a, a))


def trace(arg, axis_offset=0):
    """
    Returns the trace of ``arg`` which is the sum of ``arg[k,k]`` over k.

    :param arg: argument
    :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``
    :param axis_offset: ``axis_offset`` to components to sum over. ``axis_offset``
                        must be non-negative and less than the rank of ``arg`` +1.
                        The dimensions of component ``axis_offset`` and
                        axis_offset+1 must be equal.
    :type axis_offset: ``int``
    :return: trace of arg. The rank of the returned object is rank of ``arg``
             minus 2.
    :rtype: `escript.Data`, `Symbol` or ``numpy.ndarray`` depending on the type
            of ``arg``
    """
    if isinstance(arg, numpy.ndarray):
        sh = arg.shape
        if len(sh) < 2:
            raise ValueError("rank of argument must be greater than 1")
        if axis_offset < 0 or axis_offset > len(sh) - 2:
            raise ValueError("axis_offset must be between 0 and %d" % (len(sh) - 2))
        s1 = 1
        for i in range(axis_offset):
            s1 *= sh[i]
        s2 = 1
        for i in range(axis_offset + 2, len(sh)):
            s2 *= sh[i]
        if not sh[axis_offset] == sh[axis_offset + 1]:
            raise ValueError("dimensions of component %d and %d must match." % (axis_offset, axis_offset + 1))
        arg_reshaped = numpy.reshape(arg, (s1, sh[axis_offset], sh[axis_offset], s2))
        out = numpy.zeros([s1, s2], numpy.double)
        for i1 in range(s1):
            for i2 in range(s2):
                for j in range(sh[axis_offset]):
                    out[i1, i2] += arg_reshaped[i1, j, j, i2]
        out.resize(sh[:axis_offset] + sh[axis_offset + 2:])
        return out
    elif isinstance(arg, escore.Data):
        if arg.getRank() < 2:
            raise ValueError("rank of argument must be greater than 1")
        if axis_offset < 0 or axis_offset > arg.getRank() - 2:
            raise ValueError("axis_offset must be between 0 and %d" % (arg.getRank() - 2))
        s = list(arg.getShape())
        if not s[axis_offset] == s[axis_offset + 1]:
            raise ValueError("dimensions of component %d and %d must match." % (axis_offset, axis_offset + 1))
        return arg._trace(axis_offset)
    elif isinstance(arg, sym.Symbol):
        if arg.getRank() < 2:
            raise ValueError("rank of argument must be greater than 1")
        if axis_offset < 0 or axis_offset > arg.getRank() - 2:
            raise ValueError("axis_offset must be between 0 and %d" % (arg.getRank() - 2))
        s = list(arg.getShape())
        if not s[axis_offset] == s[axis_offset + 1]:
            raise ValueError("dimensions of component %d and %d must match." % (axis_offset, axis_offset + 1))
        return arg.trace(axis_offset)
    elif isinstance(arg, complex):
        raise TypeError("illegal argument type complex.")
    elif isinstance(arg, float):
        raise TypeError("illegal argument type float.")
    elif isinstance(arg, int):
        raise TypeError("illegal argument type int.")
    else:
        raise TypeError("Unknown argument type.")


def transpose(arg, axis_offset=None):
    """
    Returns the transpose of ``arg`` by swapping the first ``axis_offset`` and
    the last ``rank-axis_offset`` components.

    :param arg: argument
    :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``, ``float``, ``int``
    :param axis_offset: the first ``axis_offset`` components are swapped with the
                        rest. ``axis_offset`` must be non-negative and less or
                        equal to the rank of ``arg``. If ``axis_offset`` is not
                        present ``int(r/2)`` where r is the rank of ``arg`` is
                        used.
    :type axis_offset: ``int``
    :return: transpose of ``arg``
    :rtype: `escript.Data`, `Symbol`, ``numpy.ndarray``, ``float``, ``int``
            depending on the type of ``arg``
    """
    if isinstance(arg, numpy.ndarray):
        if axis_offset is None:
            axis_offset = int(arg.ndim / 2)
        return numpy.transpose(arg, axes=list(range(axis_offset, arg.ndim)) + list(range(0, axis_offset)))
    elif isinstance(arg, escore.Data):
        r = arg.getRank()
        if axis_offset is None:
            axis_offset = int(r / 2)
        if axis_offset < 0 or axis_offset > r:
            raise ValueError("axis_offset must be between 0 and %s" % r)
        return arg._transpose(axis_offset)
    elif isinstance(arg, complex):
        if not (axis_offset == 0 or axis_offset is None):
            raise ValueError("axis_offset must be 0 for complex argument")
        return arg
    elif isinstance(arg, float):
        if not (axis_offset == 0 or axis_offset is None):
            raise ValueError("axis_offset must be 0 for float argument")
        return arg
    elif isinstance(arg, int):
        if not (axis_offset == 0 or axis_offset is None):
            raise ValueError("axis_offset must be 0 for int argument")
        return float(arg)
    elif isinstance(arg, sym.Symbol):
        r = arg.getRank()
        if axis_offset is None:
            axis_offset = int(r / 2)
        if axis_offset < 0 or axis_offset > r:
            raise ValueError("axis_offset must be between 0 and %s" % r)
        return arg.transpose(axis_offset)
    else:
        raise TypeError("Unknown argument type.")


def swap_axes(arg, axis0=0, axis1=1):
    """
    Returns the swap of ``arg`` by swapping the components ``axis0`` and ``axis1``.

    :param arg: argument
    :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``
    :param axis0: first axis. ``axis0`` must be non-negative and less than the
                  rank of ``arg``.
    :type axis0: ``int``
    :param axis1: second axis. ``axis1`` must be non-negative and less than the
                  rank of ``arg``.
    :type axis1: ``int``
    :return: ``arg`` with swapped components
    :rtype: `escript.Data`, `Symbol` or ``numpy.ndarray`` depending on the type
            of ``arg``
    """
    if axis0 > axis1:
        axis0, axis1 = axis1, axis0
    if isinstance(arg, numpy.ndarray):
        return numpy.swapaxes(arg, axis0, axis1)
    elif isinstance(arg, escore.Data):
        return arg._swap_axes(axis0, axis1)
    elif isinstance(arg, sym.Symbol):
        return arg.swap_axes(axis0, axis1)
    elif isinstance(arg, complex):
        raise TypeError("complex argument is not supported.")
    elif isinstance(arg, float):
        raise TypeError("float argument is not supported.")
    elif isinstance(arg, int):
        raise TypeError("int argument is not supported.")
    else:
        raise TypeError("Unknown argument type.")


def symmetric(arg):
    """
    Returns the symmetric part of the square matrix ``arg``. That is,
    *(arg+transpose(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: symmetric part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg, numpy.ndarray):
        if arg.ndim == 2:
            if not (arg.shape[0] == arg.shape[1]):
                raise ValueError("argument must be square.")
        elif arg.ndim == 4:
            if not (arg.shape[0] == arg.shape[2] and arg.shape[1] == arg.shape[3]):
                raise ValueError("argument must be square.")
        else:
            raise ValueError("rank 2 or 4 is required.")
        return (arg + transpose(arg)) / 2
    elif isinstance(arg, escore.Data):
        if arg.getRank() == 2:
            if not (arg.getShape()[0] == arg.getShape()[1]):
                raise ValueError("argument must be square.")
            return arg._symmetric()
        elif arg.getRank() == 4:
            if not (arg.getShape()[0] == arg.getShape()[2] and arg.getShape()[1] == arg.getShape()[3]):
                raise ValueError("argument must be square.")
            return arg._symmetric()
        else:
            raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank() == 2:
            if arg.getShape()[0] != arg.getShape()[1]:
                raise ValueError("symmetric: argument must be square.")
        elif arg.getRank() == 4:
            if arg.getShape()[0] != arg.getShape()[2] or arg.getShape()[1] != arg.getShape()[3]:
                raise ValueError("symmetric: argument must be square.")
        else:
            raise ValueError("symmetric: rank 2 or 4 is required.")
        return (arg + transpose(arg)) / 2
    elif isinstance(arg, complex):
        return arg
    elif isinstance(arg, float):
        return arg
    elif isinstance(arg, int):
        return float(arg)
    else:
        raise TypeError("symmetric: Unknown argument type.")


def nonsymmetric(arg):
    """
    Deprecated alias for antisymmetric.
    """
    return antisymmetric(arg)


def antisymmetric(arg):
    """
    Returns the anti-symmetric part of the square matrix ``arg``. That is,
    *(arg-transpose(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: anti-symmetric part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg, numpy.ndarray):
        if arg.ndim == 2:
            if not (arg.shape[0] == arg.shape[1]):
                raise ValueError("antisymmetric: argument must be square.")
        elif arg.ndim == 4:
            if not (arg.shape[0] == arg.shape[2] and arg.shape[1] == arg.shape[3]):
                raise ValueError("antisymmetric: argument must be square.")
        else:
            raise ValueError("antisymmetric: rank 2 or 4 is required.")
        return (arg - transpose(arg)) / 2
    elif isinstance(arg, escore.Data):
        if arg.getRank() == 2:
            if not (arg.getShape()[0] == arg.getShape()[1]):
                raise ValueError("argument must be square.")
            return arg._antisymmetric()
        elif arg.getRank() == 4:
            if not (arg.getShape()[0] == arg.getShape()[2] and arg.getShape()[1] == arg.getShape()[3]):
                raise ValueError("argument must be square.")
            return arg._antisymmetric()
        else:
            raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank() == 2:
            if arg.getShape()[0] != arg.getShape()[1]:
                raise ValueError("antisymmetric: argument must be square.")
        elif arg.getRank() == 4:
            if arg.getShape()[0] != arg.getShape()[2] or arg.getShape()[1] != arg.getShape()[3]:
                raise ValueError("antisymmetric: argument must be square.")
        else:
            raise ValueError("antisymmetric: rank 2 or 4 is required.")
        return (arg - transpose(arg)) / 2
    elif isinstance(arg, complex):
        return complex(0)
    elif isinstance(arg, float):
        return float(0)
    elif isinstance(arg, int):
        return float(0)
    else:
        raise TypeError("antisymmetric: Unknown argument type.")


def hermitian(arg):
    """
    Returns the hermitian part of the square matrix ``arg``. That is,
    *(arg+adjoint(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: hermitian part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg, numpy.ndarray):
        if arg.ndim == 2:
            if not (arg.shape[0] == arg.shape[1]):
                raise ValueError("argument must be square.")
        elif arg.ndim == 4:
            if not (arg.shape[0] == arg.shape[2] and arg.shape[1] == arg.shape[3]):
                raise ValueError("argument must be square.")
        else:
            raise ValueError("rank 2 or 4 is required.")
        return (arg + transpose(arg).conj()) / 2
    elif isinstance(arg, escore.Data):
        if arg.getRank() == 2:
            if not (arg.getShape()[0] == arg.getShape()[1]):
                raise ValueError("argument must be square.")
            return arg._hermitian()
        elif arg.getRank() == 4:
            if not (arg.getShape()[0] == arg.getShape()[2] and arg.getShape()[1] == arg.getShape()[3]):
                raise ValueError("argument must be square.")
            return arg._hermitian()
        else:
            raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank() == 2:
            if arg.getShape()[0] != arg.getShape()[1]:
                raise ValueError("hermitian: argument must be square.")
        elif arg.getRank() == 4:
            if arg.getShape()[0] != arg.getShape()[2] or arg.getShape()[1] != arg.getShape()[3]:
                raise ValueError("hermitian: argument must be square.")
        else:
            raise ValueError("hermitian: rank 2 or 4 is required.")
        # Note: adjoint not available, using transpose as approximation
        return (arg + transpose(arg)) / 2
    elif isinstance(arg, complex):
        return complex(arg.real)
    elif isinstance(arg, float):
        return arg
    elif isinstance(arg, int):
        return float(arg)
    else:
        raise TypeError("hermitian: Unknown argument type.")


def antihermitian(arg):
    """
    Returns the anti-hermitian part of the square matrix ``arg``. That is,
    *(arg-adjoint(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: anti-hermitian part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg, numpy.ndarray):
        if arg.ndim == 2:
            if not (arg.shape[0] == arg.shape[1]):
                raise ValueError("antihermitian: argument must be square.")
        elif arg.ndim == 4:
            if not (arg.shape[0] == arg.shape[2] and arg.shape[1] == arg.shape[3]):
                raise ValueError("antihermitian: argument must be square.")
        else:
            raise ValueError("antihermitian: rank 2 or 4 is required.")
        return (arg - transpose(arg).conj()) / 2
    elif isinstance(arg, escore.Data):
        if arg.getRank() == 2:
            if not (arg.getShape()[0] == arg.getShape()[1]):
                raise ValueError("argument must be square.")
            return arg._antihermitian()
        elif arg.getRank() == 4:
            if not (arg.getShape()[0] == arg.getShape()[2] and arg.getShape()[1] == arg.getShape()[3]):
                raise ValueError("argument must be square.")
            return arg._antihermitian()
        else:
            raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank() == 2:
            if arg.getShape()[0] != arg.getShape()[1]:
                raise ValueError("antihermitian: argument must be square.")
        elif arg.getRank() == 4:
            if arg.getShape()[0] != arg.getShape()[2] or arg.getShape()[1] != arg.getShape()[3]:
                raise ValueError("antihermitian: argument must be square.")
        else:
            raise ValueError("antihermitian: rank 2 or 4 is required.")
        return (arg - hermitian(arg)) / 2
    elif isinstance(arg, complex):
        return complex(arg.imag * 1j)
    elif isinstance(arg, float):
        return float(0)
    elif isinstance(arg, int):
        return float(0)
    else:
        raise TypeError("antihermitian: Unknown argument type.")


def inverse(arg):
    """
    Returns the inverse of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: inverse of the argument. ``matrix_mult(inverse(arg),arg)`` will be
             almost equal to ``kronecker(arg.getShape()[0])``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :note: for `escript.Data` objects the dimension is restricted to 3.
    """
    import numpy.linalg
    if isinstance(arg, numpy.ndarray):
        return numpy.linalg.tensorinv(arg, ind=1)
    elif isinstance(arg, escore.Data):
        return escript_inverse(arg)
    elif isinstance(arg, complex):
        return 1. / arg
    elif isinstance(arg, float):
        return 1. / arg
    elif isinstance(arg, int):
        return 1. / float(arg)
    elif isinstance(arg, sym.Symbol):
        return arg.inverse()
    else:
        raise TypeError("inverse: Unknown argument type.")


def escript_inverse(arg):
    """
    Returns the inverse of an escript Data matrix.

    :param arg: square matrix Data object (dimension restricted to 3)
    :type arg: `escript.Data`
    :return: inverse of the matrix
    :rtype: `escript.Data`
    """
    return arg._inverse()


def eigenvalues(arg):
    """
    Returns the eigenvalues of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                ``transpose(arg)==arg`` (this is not checked).
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the eigenvalues in increasing order
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :note: for `escript.Data` and `Symbol` objects the dimension is
           restricted to 3.
    """
    if isinstance(arg, numpy.ndarray):
        out = numpy.linalg.eigvals((arg + numpy.transpose(arg)) / 2.)
        out.sort()
        return out
    elif isinstance(arg, escore.Data):
        return arg._eigenvalues()
    elif isinstance(arg, complex):
        return arg
    elif isinstance(arg, float):
        return arg
    elif isinstance(arg, int):
        return float(arg)
    elif isinstance(arg, sym.Symbol):
        return sym.symfn.eigenvalues(arg)
    else:
        raise TypeError("eigenvalues: Unknown argument type.")


def eigenvalues_and_eigenvectors(arg):
    """
    Returns the eigenvalues and eigenvectors of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                ``transpose(arg)==arg`` (this is not checked).
    :type arg: `escript.Data`
    :return: the eigenvalues and eigenvectors. The eigenvalues are ordered by
             increasing value. The eigenvectors are orthogonal and normalized.
             If V are the eigenvectors then V[:,i] is the eigenvector
             corresponding to the i-th eigenvalue.
    :rtype: `tuple` of `escript.Data`
    :note: The dimension is restricted to 3.
    """
    if isinstance(arg, numpy.ndarray):
        raise TypeError("eigenvalues_and_eigenvectors does not support numpy.ndarray arguments")
    elif isinstance(arg, escore.Data):
        return arg._eigenvalues_and_eigenvectors()
    elif isinstance(arg, complex):
        return (numpy.array([[arg]], numpy.cdouble_), numpy.ones((1, 1), numpy.cdouble_))
    elif isinstance(arg, float):
        return (numpy.array([[arg]], numpy.float_), numpy.ones((1, 1), numpy.float_))
    elif isinstance(arg, int):
        return (numpy.array([[arg]], numpy.float_), numpy.ones((1, 1), numpy.float_))
    elif isinstance(arg, sym.Symbol):
        return sym.symfn.eigenvalues_and_eigenvectors(arg)
    else:
        raise TypeError("eigenvalues: Unknown argument type.")


# =========================================================================
#   Tensor Arithmetic
# =========================================================================

def mult(arg0, arg1):
    """
    Product of ``arg0`` and ``arg1``.

    :param arg0: first term
    :type arg0: `Symbol`, ``float``, ``int``, `escript.Data` or ``numpy.ndarray``
    :param arg1: second term
    :type arg1: `Symbol`, ``float``, ``int``, `escript.Data` or ``numpy.ndarray``
    :return: the product of ``arg0`` and ``arg1``
    :rtype: `Symbol`, ``float``, ``int``, `escript.Data` or ``numpy.ndarray``
    :note: The shape of both arguments is matched according to the rules
           used in `matchShape`.
    """
    # Import locally to avoid circular import
    from . import util
    args = util.matchShape(arg0, arg1)
    if util.testForZero(args[0]) or util.testForZero(args[1]):
        return numpy.zeros(util.getShape(args[0]), numpy.double)
    else:
        if isinstance(args[0], numpy.ndarray):
            return args[1] * args[0]
        else:
            return args[0] * args[1]


def maximum(*args):
    """
    The maximum over arguments ``args``.

    :param args: arguments
    :type args: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
    :return: an object which in each entry gives the maximum of the
             corresponding values in ``args``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
            depending on the input
    """
    if max([isinstance(v, sym.Symbol) for v in args]):
        return sym.symfn.maximum(*args)
    out = None
    for a in args:
        if out is None:
            out = a * 1.
        else:
            if isinstance(out, escore.Data) and isinstance(a, escore.Data):
                if out.getRank() == 0 and a.getRank() > 0:
                    res = a.copy()
                    res.copyWithMask(out, wherePositive(out - a))
                    out = res
                else:
                    out.copyWithMask(a, wherePositive(a - out))
            else:
                if isinstance(a, numpy.ndarray):
                    diff = -out + a
                else:
                    diff = a - out
                temp = mult(whereNonPositive(diff), out) + mult(wherePositive(diff), a)
                if isinstance(out, numpy.ndarray) and isinstance(a, numpy.ndarray):
                    temp = numpy.array(temp)
                out = temp
    return out


def minimum(*args):
    """
    The minimum over arguments ``args``.

    :param args: arguments
    :type args: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
    :return: an object which gives in each entry the minimum of the
             corresponding values in ``args``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
            depending on the input
    """
    if max([isinstance(v, sym.Symbol) for v in args]):
        return sym.symfn.minimum(*args)
    out = None
    for a in args:
        if out is None:
            if isinstance(a, numpy.ndarray):
                out = a.copy()
            else:
                out = a * 1.
        else:
            if isinstance(out, escore.Data) and isinstance(a, escore.Data):
                if out.getRank() == 0 and a.getRank() > 0:
                    res = a.copy()
                    res.copyWithMask(out, whereNegative(out - a))
                    out = res
                else:
                    out.copyWithMask(a, whereNegative(a - out))
            else:
                if isinstance(a, numpy.ndarray):
                    diff = -out + a
                else:
                    diff = a - out
                temp = mult(whereNonNegative(diff), out) + mult(whereNegative(diff), a)
                if isinstance(out, numpy.ndarray) and isinstance(a, numpy.ndarray):
                    temp = numpy.array(temp)
                out = temp
    return out


def clip(arg, minval=None, maxval=None):
    """
    Cuts the values of ``arg`` between ``minval`` and ``maxval``.

    :param arg: argument
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
    :param minval: lower range. If None no lower range is applied
    :type minval: ``float`` or ``None``
    :param maxval: upper range. If None no upper range is applied
    :type maxval: ``float`` or ``None``
    :return: an object that contains all values from ``arg`` between ``minval``
             and ``maxval``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or ``float``
            depending on the input
    :raise ValueError: if ``minval>maxval``
    """
    if isinstance(arg, sym.Symbol):
        clip_item = lambda item: sym.symfn.clip(item, minval, maxval)
        return arg.applyfunc(clip_item)
    if minval is not None and maxval is not None:
        if minval > maxval:
            raise ValueError("minval = %s must be less than maxval %s" % (minval, maxval))
    if minval is None:
        tmp = arg
    else:
        tmp = maximum(minval, arg)
    if maxval is None:
        return tmp
    else:
        return minimum(tmp, maxval)


def cross(arg0, arg1):
    """
    Cross product of the two arguments ``arg0`` and ``arg1`` which need to be shape (3,).

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the cross product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on ``arg0``
    :raise ValueError: if the shapes of the arguments are not (3,)
    """
    from . import util
    if isinstance(arg0, numpy.ndarray) and isinstance(arg1, numpy.ndarray):
        out = numpy.cross(arg0, arg1)
    else:
        sh0 = util.getShape(arg0)
        sh1 = util.getShape(arg1)
        if not sh0 == (3,):
            raise ValueError("cross: arg0 needs to be of shape (3,)")
        if not sh1 == (3,):
            raise ValueError("cross: arg1 needs to be of shape (3,)")

        if isinstance(arg0, sym.Symbol):
            if isinstance(arg1, sym.Symbol):
                out = sym.Symbol(arg0.name + "x" + arg1.name, (3,))
            else:
                out = sym.Symbol(arg0.name + "x" + str(type(arg1)), (3,))
        elif isinstance(arg0, escore.Data):
            if isinstance(arg1, sym.Symbol):
                out = sym.Symbol(str(type(arg0)) + "x" + arg1.name, (3,))
            else:
                out = escore.Data(0., (3,), arg0.getFunctionSpace())
        elif isinstance(arg1, escore.Data):
            if isinstance(arg0, sym.Symbol):
                out = sym.Symbol(str(type(arg0)) + "x" + arg1.name, (3,))
            else:
                out = escore.Data(0., (3,), arg1.getFunctionSpace())
        else:
            raise TypeError("cross: argument type not supported")

        out[0] = arg0[1] * arg1[2] - arg0[2] * arg1[1]
        out[1] = arg0[2] * arg1[0] - arg0[0] * arg1[2]
        out[2] = arg0[0] * arg1[1] - arg0[1] * arg1[0]
    return out


# =========================================================================
#   Tensor Products
# =========================================================================

def inner(arg0, arg1):
    """
    Inner product of the two arguments. The inner product is defined as:

    `out=Sigma_s arg0[s]*arg1[s]`

    where s runs through ``arg0.Shape``.

    ``arg0`` and ``arg1`` must have the same shape.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the inner product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``
            depending on the input
    :raise ValueError: if the shapes of the arguments are not identical
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if not sh0 == sh1:
        raise ValueError("inner: shape of arguments does not match")
    return generalTensorProduct(arg0, arg1, axis_offset=len(sh0))


def outer(arg0, arg1):
    """
    The outer product of the two arguments. The outer product is defined as:

    ``out[t,s]=arg0[t]*arg1[s]``

    where
        - s runs through ``arg0.Shape``
        - t runs through ``arg1.Shape``

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the outer product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    return generalTensorProduct(arg0, arg1, axis_offset=0)


def matrixmult(arg0, arg1):
    """
    See `matrix_mult`.
    """
    return matrix_mult(arg0, arg1)


def matrix_mult(arg0, arg1):
    """
    matrix-matrix or matrix-vector product of the two arguments.

    `out[s0]=Sigma_{r0} arg0[s0,r0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[r0,s1]`

    The second dimension of ``arg0`` and the first dimension of ``arg1`` must
    match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of at least rank 1
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the matrix-matrix or matrix-vector product of ``arg0`` and ``arg1``
             at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if not len(sh0) == 2:
        raise ValueError("first argument must have rank 2")
    if not len(sh1) == 2 and not len(sh1) == 1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTensorProduct(arg0, arg1, axis_offset=1)


def tensormult(arg0, arg1):
    """
    See `tensor_mult`.
    """
    return tensor_mult(arg0, arg1)


def tensor_mult(arg0, arg1):
    """
    The tensor product of the two arguments.

    For ``arg0`` of rank 2 this is

    `out[s0]=Sigma_{r0} arg0[s0,r0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[r0,s1]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2,s3]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2]`

    or

    `out[s0,s1]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1]`

    In the first case the second dimension of ``arg0`` and the last dimension of
    ``arg1`` must match and in the second case the two last dimensions of ``arg0``
    must match the two first dimensions of ``arg1``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater than 1 or 2 depending on the
                 rank of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if len(sh0) == 2 and (len(sh1) == 2 or len(sh1) == 1):
        return generalTensorProduct(arg0, arg1, axis_offset=1)
    elif len(sh0) == 4 and (len(sh1) == 2 or len(sh1) == 3 or len(sh1) == 4):
        return generalTensorProduct(arg0, arg1, axis_offset=2)
    else:
        raise ValueError("tensor_mult: first argument must have rank 2 or 4 and second rank must be in (1,2) or (2,3,4) respectively.")


def generalTensorProduct(arg0, arg1, axis_offset=0):
    """
    Generalized tensor product.

    `out[s,t]=Sigma_r arg0[s,r]*arg1[r,t]`

    where
        - s runs through ``arg0.Shape[:arg0.ndim-axis_offset]``
        - r runs through ``arg1.Shape[:axis_offset]``
        - t runs through ``arg1.Shape[axis_offset:]``

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``arg0`` and ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    if (isinstance(arg0, float) or isinstance(arg0, complex)) and (isinstance(arg1, float) or isinstance(arg1, complex)):
        return arg1 * arg0
    arg0, arg1 = util.matchType(arg0, arg1)
    if isinstance(arg0, sym.Symbol):
        sh0 = arg0.getShape()
        sh1 = util.getShape(arg1)
        if not sh0[arg0.getRank() - axis_offset:] == sh1[:axis_offset]:
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        if isinstance(arg1, float):
            return arg0 * arg1
        elif isinstance(arg1, numpy.ndarray) or isinstance(arg1, sym.Symbol):
            return arg0.tensorProduct(arg1, axis_offset)
        elif isinstance(arg1, escore.Data):
            raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0, numpy.ndarray):
        if not arg0.shape[arg0.ndim - axis_offset:] == arg1.shape[:axis_offset]:
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        arg0_c = arg0.copy()
        arg1_c = arg1.copy()
        sh0, sh1 = arg0.shape, arg1.shape
        d0, d1, d01 = 1, 1, 1
        for i in sh0[:arg0.ndim - axis_offset]:
            d0 *= i
        for i in sh1[axis_offset:]:
            d1 *= i
        for i in sh1[:axis_offset]:
            d01 *= i
        arg0_c.resize((d0, d01))
        arg1_c.resize((d01, d1))
        if arg0_c.dtype.kind == 'c':
            restype = arg0_c.dtype
        else:
            restype = arg1_c.dtype
        out = numpy.zeros((d0, d1), restype)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0, i1] = numpy.sum(arg0_c[i0, :] * arg1_c[:, i1])
        out.resize(sh0[:arg0.ndim - axis_offset] + sh1[axis_offset:])
        return out
    elif isinstance(arg0, escore.Data):
        if isinstance(arg1, sym.Symbol):
            raise TypeError("tensor product of Data and Symbol not supported yet")
        return escript_generalTensorProduct(arg0, arg1, axis_offset)
    raise TypeError("generalTensorProduct: Unsupported argument type")


def escript_generalTensorProduct(arg0, arg1, axis_offset, transpose=0):
    """
    Generalized tensor product for escript Data objects.

    :param arg0: first Data object
    :type arg0: `escript.Data`
    :param arg1: second Data object (not necessarily on the same function space)
    :type arg1: `escript.Data`
    :param axis_offset: axis offset for the tensor product
    :type axis_offset: ``int``
    :param transpose: transpose mode (0=none, 1=transpose arg0, 2=transpose arg1)
    :type transpose: ``int``
    :return: tensor product result
    :rtype: `escript.Data`
    """
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, transpose)


def transposed_matrix_mult(arg0, arg1):
    """
    transposed(matrix)-matrix or transposed(matrix)-vector product of the two
    arguments.

    `out[s0]=Sigma_{r0} arg0[r0,s0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[r0,s0]*arg1[r0,s1]`

    The function call ``transposed_matrix_mult(arg0,arg1)`` is equivalent to
    ``matrix_mult(transpose(arg0),arg1)``.

    The first dimension of ``arg0`` and ``arg1`` must match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of at least rank 1
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the product of the transpose of ``arg0`` and ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if not len(sh0) == 2:
        raise ValueError("first argument must have rank 2")
    if not len(sh1) == 2 and not len(sh1) == 1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTransposedTensorProduct(arg0, arg1, axis_offset=1)


def transposed_tensor_mult(arg0, arg1):
    """
    The tensor product of the transpose of the first and the second argument.

    For ``arg0`` of rank 2 this is

    `out[s0]=Sigma_{r0} arg0[r0,s0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[r0,s0]*arg1[r0,s1]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2,s3]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2]`

    or

    `out[s0,s1]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1]`

    In the first case the first dimension of ``arg0`` and the first dimension of
    ``arg1`` must match and in the second case the two first dimensions of
    ``arg0`` must match the two first dimensions of ``arg1``.

    The function call ``transposed_tensor_mult(arg0,arg1)`` is equivalent to
    ``tensor_mult(transpose(arg0),arg1)``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater of 1 or 2 depending on the
                 rank of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of transpose of arg0 and arg1 at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if len(sh0) == 2 and (len(sh1) == 2 or len(sh1) == 1):
        return generalTransposedTensorProduct(arg0, arg1, axis_offset=1)
    elif len(sh0) == 4 and (len(sh1) == 2 or len(sh1) == 3 or len(sh1) == 4):
        return generalTransposedTensorProduct(arg0, arg1, axis_offset=2)
    else:
        raise ValueError("first argument must have rank 2 or 4")


def generalTransposedTensorProduct(arg0, arg1, axis_offset=0):
    """
    Generalized tensor product of transposed of ``arg0`` and ``arg1``.

    `out[s,t]=Sigma_r arg0[r,s]*arg1[r,t]`

    where
        - s runs through ``arg0.Shape[axis_offset:]``
        - r runs through ``arg0.Shape[:axis_offset]``
        - t runs through ``arg1.Shape[axis_offset:]``

    The function call ``generalTransposedTensorProduct(arg0,arg1,axis_offset)``
    is equivalent to
    ``generalTensorProduct(transpose(arg0,arg0.ndim-axis_offset),arg1,axis_offset)``.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``transpose(arg0)`` and ``arg1`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    if (isinstance(arg0, float) and isinstance(arg1, float)) or (isinstance(arg0, complex) and isinstance(arg1, complex)):
        return arg1 * arg0
    arg0, arg1 = util.matchType(arg0, arg1)
    if isinstance(arg0, sym.Symbol):
        sh0 = arg0.getShape()
        sh1 = util.getShape(arg1)
        if not sh0[:axis_offset] == sh1[:axis_offset]:
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        if isinstance(arg1, float):
            return arg0 * arg1
        elif isinstance(arg1, numpy.ndarray) or isinstance(arg1, sym.Symbol):
            return arg0.transposedTensorProduct(arg1, axis_offset)
        elif isinstance(arg1, escore.Data):
            raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0, numpy.ndarray):
        if not arg0.shape[:axis_offset] == arg1.shape[:axis_offset]:
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        arg0_c = arg0.copy()
        arg1_c = arg1.copy()
        sh0, sh1 = arg0.shape, arg1.shape
        d0, d1, d01 = 1, 1, 1
        for i in sh0[axis_offset:]:
            d0 *= i
        for i in sh1[axis_offset:]:
            d1 *= i
        for i in sh0[:axis_offset]:
            d01 *= i
        arg0_c.resize((d01, d0))
        arg1_c.resize((d01, d1))
        target_type = arg0.dtype if arg0.dtype.kind == 'c' else arg1.dtype
        out = numpy.zeros((d0, d1), target_type)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0, i1] = numpy.sum(arg0_c[:, i0] * arg1_c[:, i1])
        out.resize(sh0[axis_offset:] + sh1[axis_offset:])
        return out
    elif isinstance(arg0, escore.Data):
        if isinstance(arg1, sym.Symbol):
            raise TypeError("tensor product of Data and Symbol not supported yet")
        return escript_generalTransposedTensorProduct(arg0, arg1, axis_offset)


def escript_generalTransposedTensorProduct(arg0, arg1, axis_offset):
    """
    Generalized tensor product of transposed arg0 and arg1 for escript Data objects.

    :param arg0: first Data object (will be transposed)
    :type arg0: `escript.Data`
    :param arg1: second Data object (not necessarily on the same function space)
    :type arg1: `escript.Data`
    :param axis_offset: axis offset for the tensor product
    :type axis_offset: ``int``
    :return: tensor product of transpose(arg0) and arg1
    :rtype: `escript.Data`
    """
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 1)


def matrix_transposed_mult(arg0, arg1):
    """
    matrix-transposed(matrix) product of the two arguments.

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[s1,r0]`

    The function call ``matrix_transposed_mult(arg0,arg1)`` is equivalent to
    ``matrix_mult(arg0,transpose(arg1))``.

    The last dimensions of ``arg0`` and ``arg1`` must match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of rank 1 or 2
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the product of ``arg0`` and the transposed of ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if not len(sh0) == 2:
        raise ValueError("first argument must have rank 2")
    if not len(sh1) == 2 and not len(sh1) == 1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTensorTransposedProduct(arg0, arg1, axis_offset=1)


def tensor_transposed_mult(arg0, arg1):
    """
    The tensor product of the first and the transpose of the second argument.

    For ``arg0`` of rank 2 this is

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[s1,r0]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,s3,r0,r1]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,r0,r1]`

    In the first case the second dimension of ``arg0`` and ``arg1`` must
    match and in the second case the two last dimensions of ``arg0`` must match
    the two last dimensions of ``arg1``.

    The function call ``tensor_transpose_mult(arg0,arg1)`` is equivalent to
    ``tensor_mult(arg0,transpose(arg1))``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater of 1 or 2 depending on rank
                 of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of ``arg0`` and the transposed of ``arg1`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    sh0 = util.getShape(arg0)
    sh1 = util.getShape(arg1)
    if len(sh0) == 2 and (len(sh1) == 2 or len(sh1) == 1):
        return generalTensorTransposedProduct(arg0, arg1, axis_offset=1)
    elif len(sh0) == 4 and (len(sh1) == 2 or len(sh1) == 3 or len(sh1) == 4):
        if len(sh1) == 2:
            return generalTensorTransposedProduct(arg0, transpose(arg1), axis_offset=2)
        else:
            return generalTensorTransposedProduct(arg0, arg1, axis_offset=2)
    else:
        raise ValueError("first argument must have rank 2 or 4")


def generalTensorTransposedProduct(arg0, arg1, axis_offset=0):
    """
    Generalized tensor product of ``arg0`` and transpose of ``arg1``.

    `out[s,t]=Sigma_r arg0[s,r]*arg1[t,r]`

    where
        - s runs through ``arg0.Shape[:arg0.ndim-axis_offset]``
        - r runs through ``arg0.Shape[arg1.ndim-axis_offset:]``
        - t runs through ``arg1.Shape[arg1.ndim-axis_offset:]``

    The function call ``generalTensorTransposedProduct(arg0,arg1,axis_offset)``
    is equivalent to
    ``generalTensorProduct(arg0,transpose(arg1,arg1.ndim-axis_offset),axis_offset)``.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``arg0`` and ``transpose(arg1)`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    from . import util
    if ((isinstance(arg0, float) or isinstance(arg0, complex)) and
            (isinstance(arg1, float) or isinstance(arg1, complex))):
        return arg1 * arg0
    arg0, arg1 = util.matchType(arg0, arg1)
    if isinstance(arg0, sym.Symbol):
        sh0 = arg0.getShape()
        sh1 = util.getShape(arg1)
        r1 = util.getRank(arg1)
        if not sh0[arg0.getRank() - axis_offset:] == sh1[r1 - axis_offset:]:
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        if isinstance(arg1, float) or isinstance(arg1, complex):
            return arg0 * arg1
        elif isinstance(arg1, numpy.ndarray) or isinstance(arg1, sym.Symbol):
            return arg0.tensorTransposedProduct(arg1, axis_offset)
        elif isinstance(arg1, escore.Data):
            raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0, numpy.ndarray):
        if not (arg0.shape[arg0.ndim - axis_offset:] == arg1.shape[arg1.ndim - axis_offset:] or
                arg0.shape[arg0.ndim - axis_offset:], tuple(reversed(arg1.shape[arg1.ndim - axis_offset:]))):
            raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument." % (axis_offset, axis_offset))
        arg0_c = arg0.copy()
        arg1_c = arg1.copy()
        sh0, sh1 = arg0.shape, arg1.shape
        d0, d1, d01 = 1, 1, 1
        for i in sh0[:arg0.ndim - axis_offset]:
            d0 *= i
        for i in sh1[:arg1.ndim - axis_offset]:
            d1 *= i
        for i in sh1[arg1.ndim - axis_offset:]:
            d01 *= i
        arg0_c.resize((d0, d01))
        arg1_c.resize((d1, d01))
        if arg0_c.dtype != numpy.double:
            out = numpy.zeros((d0, d1), arg0_c.dtype)
        else:
            out = numpy.zeros((d0, d1), numpy.double)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0, i1] = numpy.sum(arg0_c[i0, :] * arg1_c[i1, :])
        out.resize(sh0[:arg0.ndim - axis_offset] + sh1[:arg1.ndim - axis_offset])
        return out
    elif isinstance(arg0, escore.Data):
        if isinstance(arg1, sym.Symbol):
            raise TypeError("tensor product of Data and Symbol not supported yet")
        return escript_generalTensorTransposedProduct(arg0, arg1, axis_offset)


def escript_generalTensorTransposedProduct(arg0, arg1, axis_offset):
    """
    Generalized tensor product of arg0 and transposed arg1 for escript Data objects.

    :param arg0: first Data object
    :type arg0: `escript.Data`
    :param arg1: second Data object (will be transposed, not necessarily on the same function space)
    :type arg1: `escript.Data`
    :param axis_offset: axis offset for the tensor product
    :type axis_offset: ``int``
    :return: tensor product of arg0 and transpose(arg1)
    :rtype: `escript.Data`
    """
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 2)
