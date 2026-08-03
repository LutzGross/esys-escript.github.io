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
The spatial function spaces on a GRADED (hanging-node) 2D forest.

Every other oxley suite builds a forest refined uniformly, where face_code is
zero on every element and no 2:1 seam exists. That leaves the whole hanging
node path untested: grad() on a graded forest was wrong by O(1) along every
seam - all 64 gradient tests here failed, plus 5 L2 norms - while every
existing suite stayed green, because a hanging corner's slot holds a master
rather than the value at the corner.

So these are the same shared suites run on a graded forest, plus the boundary
gradients that escript's Test_Util_Gradient_noBoundary deliberately omits.
The uniform case is kept as a control: it must pass whatever the seam code
does, and its failing would point somewhere else entirely.

The domain has to be the unit square, since the shared suite asserts getX()
lies in [0,1]^dim and that integrate(x_i**k) == 1/(k+1).
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle

from test_util_spatial_functions1 import \
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact

try:
    import esys.finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

# per-block refinement levels; both graded ones carry 2:1 seams, and the
# 3x3 case has elements with 1, 2 and 3 hanging edges
MIXED = [[3, 1, 2], [1, 2, 1], [2, 1, 3]]
PEAK = [[0, 0, 1, 0, 0], [0, 1, 2, 1, 0], [1, 2, 3, 2, 1],
        [0, 1, 2, 1, 0], [0, 0, 1, 0, 0]]
UNIFORM = 2


def graded(levels):
    """the forest for a per-block level table, on the unit square"""
    if isinstance(levels, int):
        n0 = n1 = 2
    else:
        n0, n1 = len(levels), len(levels[0])
    return Rectangle(n0=n0, n1=n1, l0=1., l1=1., refine_level=levels)


class Test_SpatialFunctionsOnGradedOxley2D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order = 1
        self.domain = graded(MIXED)

    def tearDown(self):
        del self.order
        del self.domain


class Test_SpatialFunctionsOnUniformOxley2D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """control: no seam anywhere, so the hanging path is not involved"""
    def setUp(self):
        self.order = 1
        self.domain = graded(UNIFORM)

    def tearDown(self):
        del self.order
        del self.domain


class Test_GradientOnBoundaryOnGradedOxley2D(unittest.TestCase):
    """
    Gradients on the boundary function spaces of a graded forest.

    Uncovered by the class above, which is the _noBoundary_ variant. It needs
    its own test because the boundary gradient combines the two values on the
    face with the element's other two corners to get the tangential
    derivative, and on the fine side of a seam one of those can be hanging.

    A globally linear field is the discriminator: it lies in the space
    exactly, so any departure from its constant gradient means a corner value
    was read wrong, not that the mesh is too coarse.
    """
    def _check(self, dom, fs, name):
        x = ContinuousFunction(dom).getX()
        u = 2. * x[0] + 3. * x[1] - 1.
        err = Lsup(grad(u, fs) - [2., 3.])
        self.assertLess(err, 1e-8, "%s: gradient off by %g" % (name, err))

    def _allSpaces(self, levels, label):
        dom = graded(levels)
        for fs, name in ((Function(dom), "Function"),
                         (ReducedFunction(dom), "ReducedFunction"),
                         (FunctionOnBoundary(dom), "FunctionOnBoundary"),
                         (ReducedFunctionOnBoundary(dom),
                          "ReducedFunctionOnBoundary")):
            self._check(dom, fs, "%s on %s" % (name, label))

    def test_gradient_mixed(self):
        self._allSpaces(MIXED, "3x3 mixed")

    def test_gradient_peak(self):
        self._allSpaces(PEAK, "centre peak")

    def test_gradient_uniform(self):
        self._allSpaces(UNIFORM, "uniform (control)")


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_SpatialFunctionsOnExportedGradedOxley2D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """
    The same suite on the finley mesh exported from the graded forest.

    The export resolves each seam by triangulating it, so unlike the forest
    itself it is a conforming P1 space; this is what fails if the split, the
    node numbering or the export's node ownership is wrong.
    """
    def setUp(self):
        self.order = 1
        self.domain = graded(MIXED).toFinley()

    def tearDown(self):
        del self.order
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
