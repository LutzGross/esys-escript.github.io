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

from esys.oxley.RefinementZone import *


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


def Rectangle(n0=10, n1=10, l0=1.0, l1=1.0,
              diracPoints=[], diracTags=[], comm=None, framework=None,
              **legacy):
    """
    Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    _warn_legacy("Rectangle", legacy, _LEGACY_RECTANGLE_ARGS)
    dom = _Rectangle(n0=n0, n1=n1, l0=l0, l1=l1,
                     diracPoints=diracPoints, diracTags=diracTags, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


def Brick(n0=10, n1=10, n2=10, l0=1.0, l1=1.0, l2=1.0,
          diracPoints=[], diracTags=[], comm=None, framework=None,
          **legacy):
    """
    Creates a brick p4est mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param n2: number of elements in direction 2
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param l2: length of side 2 or coordinate range of side 2
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    _warn_legacy("Brick", legacy, _LEGACY_BRICK_ARGS)
    dom = _Brick(n0=n0, n1=n1, n2=n2, l0=l0, l1=l1, l2=l2,
                 diracPoints=diracPoints, diracTags=diracTags, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


__nodocorecursion=['oxleycpp']
