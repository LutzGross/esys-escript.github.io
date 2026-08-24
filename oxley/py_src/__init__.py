##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""A domain meshed with uniform rectangles or quadrilaterals
"""


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import warnings

import esys.escript       # This is just to ensure required libraries are loaded
from .oxleycpp import *
from .oxleycpp import Rectangle as _Rectangle, Brick as _Brick


# Parameters that never had any effect in oxley:
#   d0/d1/d2   - domain decomposition is handled internally by p4est, so these
#                subdivision counts were validated but never used.
#   periodic*  - periodicity was stored but never wired into the p4est
#                connectivity, so it had no effect on the mesh.
#   order      - stored in m_order but never read: the assembler always uses a
#                fixed 2-point Gauss rule (4 points in 2D, 8 in 3D, exact to
#                cubic), so the integration order is fixed regardless of this.
# They are accepted for backwards compatibility and ignored with a warning.
_LEGACY_RECTANGLE_ARGS = frozenset(['d0', 'd1', 'periodic0', 'periodic1',
                                    'order'])
_LEGACY_BRICK_ARGS = frozenset(['d0', 'd1', 'd2',
                                'periodic0', 'periodic1', 'periodic2',
                                'order'])


def _flatten_refine_level(refine_level, shape):
    """
    Normalise ``refine_level`` for the C++ factories: a scalar is passed through
    unchanged (uniform refinement), an array of shape ``shape`` (the block grid)
    is validated and flattened to a plain list in row-major (C) order, matching
    the per-block indexing used by the domain constructors.
    """
    try:
        return int(refine_level)
    except (TypeError, ValueError):
        pass
    import numpy as np
    arr = np.asarray(refine_level, dtype=int)
    if arr.shape != tuple(shape):
        raise ValueError("refine_level array shape %s does not match numBlocks %s"
                         % (arr.shape, tuple(shape)))
    return [int(v) for v in arr.reshape(-1)]


def _warn_legacy(factory, legacy, allowed):
    if not legacy:
        return
    unknown = set(legacy) - allowed
    if unknown:
        raise TypeError("%s() got an unexpected keyword argument '%s'"
                        % (factory, sorted(unknown)[0]))
    warnings.warn(
        "%s: argument(s) %s are no longer supported and are ignored. "
        "Domain decomposition is handled by p4est and periodic boundaries "
        "are not implemented." % (factory, ", ".join(sorted(legacy))),
        DeprecationWarning, stacklevel=3)


def Rectangle(n0=10, n1=10, l0=1.0, l1=1.0, refine_level=0,
              diracPoints=[], diracTags=[], comm=None, framework=None,
              **legacy):
    """
    Creates a rectangular p4est mesh of n0 x n1 blocks over the rectangle [0,l0] x [0,l1],
    each block uniformly subdivided ``refine_level`` times.

    :param n0: number of blocks in direction 0
    :param n1: number of blocks in direction 1
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param refine_level: refinement level applied to every block (a single int),
                         or an array of shape (n0, n1) giving a per-block level
                         (differing levels create hanging nodes at block seams)
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    _warn_legacy("Rectangle", legacy, _LEGACY_RECTANGLE_ARGS)
    refine_level = _flatten_refine_level(refine_level, (n0, n1))
    dom = _Rectangle(n0=n0, n1=n1, l0=l0, l1=l1, refine_level=refine_level,
                     diracPoints=diracPoints, diracTags=diracTags, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


def Brick(n0=10, n1=10, n2=10, l0=1.0, l1=1.0, l2=1.0, refine_level=0,
          diracPoints=[], diracTags=[], comm=None, framework=None,
          **legacy):
    """
    Creates a brick p8est mesh of n0 x n1 x n2 blocks over [0,l0] x [0,l1] x [0,l2],
    each block uniformly subdivided ``refine_level`` times.

    :param n0: number of blocks in direction 0
    :param n1: number of blocks in direction 1
    :param n2: number of blocks in direction 2
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param l2: length of side 2 or coordinate range of side 2
    :param refine_level: refinement level applied to every block (a single int),
                         or an array of shape (n0, n1, n2) giving a per-block level
                         (differing levels create hanging nodes at block seams)
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    _warn_legacy("Brick", legacy, _LEGACY_BRICK_ARGS)
    refine_level = _flatten_refine_level(refine_level, (n0, n1, n2))
    dom = _Brick(n0=n0, n1=n1, n2=n2, l0=l0, l1=l1, l2=l2, refine_level=refine_level,
                 diracPoints=diracPoints, diracTags=diracTags, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


def Block(numBlocks, length=None, origin=None, refine_level=0,
          diracPoints=[], diracTags=[], comm=None, framework=None):
    """
    Dimension-agnostic block-structured domain constructor (see the oxley design).

    Builds a coarse grid of ``numBlocks`` blocks (p4est/p8est trees) spanning the box
    ``[origin, origin+length]``, with every block uniformly subdivided ``refine_level``
    times. Returns a 2D (``Rectangle``) or 3D (``Brick``) domain from ``len(numBlocks)``.

    :param numBlocks: number of blocks per axis, e.g. ``(20, 20, 10)``
    :param length: physical extent per axis (default 1.0 per axis)
    :param origin: coordinate of the lower corner (default 0.0 per axis)
    :param refine_level: uniform refinement level applied to every block
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    numBlocks = tuple(numBlocks)
    dim = len(numBlocks)
    if dim not in (2, 3):
        raise ValueError("Block: numBlocks must have length 2 or 3, got %d" % dim)
    if length is None:
        length = (1.0,) * dim
    if origin is None:
        origin = (0.0,) * dim
    length = tuple(length)
    origin = tuple(origin)
    if len(length) != dim or len(origin) != dim:
        raise ValueError("Block: numBlocks, length and origin must have the same length")
    # translate origin+length into the (lower, upper) coordinate ranges the
    # Rectangle/Brick factories expect for each axis
    rng = tuple((origin[i], origin[i] + length[i]) for i in range(dim))
    common = dict(refine_level=refine_level, diracPoints=diracPoints,
                  diracTags=diracTags, comm=comm, framework=framework)
    if dim == 2:
        return Rectangle(n0=numBlocks[0], n1=numBlocks[1],
                         l0=rng[0], l1=rng[1], **common)
    return Brick(n0=numBlocks[0], n1=numBlocks[1], n2=numBlocks[2],
                 l0=rng[0], l1=rng[1], l2=rng[2], **common)


__nodocorecursion=['oxleycpp']
