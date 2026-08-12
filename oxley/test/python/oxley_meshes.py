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

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
The meshes the oxley tests run on.

REGRESSION TESTS RUN ON A GRADED MESH. A uniform forest has face_code == 0 on
every element, so it never enters the hanging-node path at all - which is how
grad() came to be wrong by O(1) along every seam while every suite stayed
green. A uniform mesh is a development convenience, not a regression test.

One uniform CONTROL is kept, and only as a discriminator: when the graded
cases fail and the control passes, the fault is in the seams rather than in
oxley generally, and that is what made the grad bug diagnosable in minutes.

Meshes come from the CONSTRUCTOR - refine_level, an int or a per-block list -
and not from a refinement operation. A per-block list is what produces 2:1
seams, which is the property being tested; how a refinement gets applied is a
separate question with its own tests.

This module is not named run_* so the build does not collect it as a suite.
"""

from esys.oxley import Rectangle, Brick

# Per-block refinement levels. Named cases, from a survey of what p4est_balance
# actually produces: a big level jump does NOT give a cell many hanging edges -
# balance grades it gradually - and 3 or 4 hanging edges arise only from an
# isolated coarse cell ringed by finer ones, so both extremes are here.
MIXED_2D = [[3, 1, 2], [1, 2, 1], [2, 1, 3]]      # 1, 2 and 3 hanging edges
SEAM_2D = [[2], [3]]                              # the simplest 2:1 seam
ISOLATED_2D = [[1, 1, 1], [1, 0, 1], [1, 1, 1]]   # 4 hanging edges on one cell
PEAK_2D = [[0, 0, 1, 0, 0], [0, 1, 2, 1, 0], [1, 2, 3, 2, 1],
           [0, 1, 2, 1, 0], [0, 0, 1, 0, 0]]
MIXED_3D = [[[2, 1], [1, 2]], [[1, 2], [2, 1]]]

UNIFORM = 2                                       # the control

# what a suite that wants more than one mesh should iterate over
CASES_2D = [("mixed", MIXED_2D), ("seam", SEAM_2D), ("isolated", ISOLATED_2D),
            ("uniform_control", UNIFORM)]


def _blocks(levels):
    """block counts implied by a per-block level table"""
    if isinstance(levels, int):
        return 2, 2
    return len(levels), len(levels[0])


def graded(levels=None, length=1., **kwargs):
    """
    A 2D forest with 2:1 seams, on [0,length]^2.

    The domain is the unit square by default because escript's shared spatial
    suite asserts getX() lies in [0,1]^dim and that integrate(x_i**k) equals
    1/(k+1).
    """
    levels = MIXED_2D if levels is None else levels
    n0, n1 = _blocks(levels)
    return Rectangle(n0=n0, n1=n1, l0=length, l1=length,
                     refine_level=levels, **kwargs)


def uniform(level=UNIFORM, blocks=2, length=1., **kwargs):
    """the control: no seam anywhere, so the hanging path is not involved"""
    return Rectangle(n0=blocks, n1=blocks, l0=length, l1=length,
                     refine_level=level, **kwargs)


def graded3D(levels=None, length=1., **kwargs):
    """a 3D forest with 2:1 seams"""
    levels = MIXED_3D if levels is None else levels
    n0 = len(levels)
    n1 = len(levels[0])
    n2 = len(levels[0][0])
    return Brick(n0=n0, n1=n1, n2=n2, l0=length, l1=length, l2=length,
                 refine_level=levels, **kwargs)


def uniform3D(level=1, blocks=2, length=1., **kwargs):
    return Brick(n0=blocks, n1=blocks, n2=blocks, l0=length, l1=length,
                 l2=length, refine_level=level, **kwargs)
