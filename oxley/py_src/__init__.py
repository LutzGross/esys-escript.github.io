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

import esys.escript       # This is just to ensure required libraries are loaded
from .oxleycpp import *
from .oxleycpp import Rectangle as _Rectangle, Brick as _Brick

from esys.oxley.RefinementZone import *


def Rectangle(n0=10, n1=10, l0=1.0, l1=1.0, d0=-1, d1=-1,
              diracPoints=[], diracTags=[], periodic0=0, periodic1=0,
              comm=None, framework=None):
    """
    Creates a rectangular p4est mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param d0: number of subdivisions in direction 0
    :param d1: number of subdivisions in direction 1
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param periodic0: periodic boundary conditions in direction 0
    :param periodic1: periodic boundary conditions in direction 1
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    dom = _Rectangle(n0=n0, n1=n1, l0=l0, l1=l1, d0=d0, d1=d1,
                     diracPoints=diracPoints, diracTags=diracTags,
                     periodic0=periodic0, periodic1=periodic1, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


def Brick(n0=10, n1=10, n2=10, l0=1.0, l1=1.0, l2=1.0, d0=-1, d1=-1, d2=-1,
          diracPoints=[], diracTags=[], periodic0=0, periodic1=0, periodic2=0,
          comm=None, framework=None):
    """
    Creates a brick p4est mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param n2: number of elements in direction 2
    :param l0: length of side 0 or coordinate range of side 0
    :param l1: length of side 1 or coordinate range of side 1
    :param l2: length of side 2 or coordinate range of side 2
    :param d0: number of subdivisions in direction 0
    :param d1: number of subdivisions in direction 1
    :param d2: number of subdivisions in direction 2
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param periodic0: periodic boundary conditions in direction 0
    :param periodic1: periodic boundary conditions in direction 1
    :param periodic2: periodic boundary conditions in direction 2
    :param comm: MPI communicator (optional, from mpi4py)
    :param framework: solver framework to use (optional SolverFramework instance)
    :return: Domain object
    """
    dom = _Brick(n0=n0, n1=n1, n2=n2, l0=l0, l1=l1, l2=l2, d0=d0, d1=d1, d2=d2,
                 diracPoints=diracPoints, diracTags=diracTags,
                 periodic0=periodic0, periodic1=periodic1, periodic2=periodic2, comm=comm)
    if framework is not None:
        dom.setFramework(framework)
    return dom


__nodocorecursion=['oxleycpp']
