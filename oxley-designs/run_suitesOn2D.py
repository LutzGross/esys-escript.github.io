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

"""
Runs escript's shared test suites against a COMPLEX (hanging-node) 2D oxley
forest and against the finley mesh exported from that same forest.

The suites in oxley/test/python only ever build a UNIFORM forest
(refine_level=1 applied to every block), so nothing there exercises a 2:1
seam. Here the same suites run on a graded forest, and then on its export -
the comparison is the point: the export is the only one of the two that is
expected to be a proper conforming P1 space.

No PDE suite is run on the oxley domain: oxley no longer assembles.

Usage:
    run-escript -n<P> -t1 run_suitesOn2D.py <oxley|finley> <group> [case]

    group: spatial | algebra | objects | all
    case:  mixed (default) | peak | uniform

The domain is the unit square in every case, because the spatial suite
asserts getX() lies in [0,1]^dim and that integrate(x_i**k) == 1/(k+1).
"""

import os
import sys

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle
# Required: without it boost::python has no registration for the concrete
# finley type and toFinley() hands back a bare Domain.
import esys.finley

from test_util import Test_util
from test_util_NaN_funcs import Test_util_NaN_funcs
from test_util_spatial_functions1 import \
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_objects import Test_Domain, Test_Dump, Test_SetDataPointValue, \
        Test_Lazy

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    print("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs(object):
        pass

# per-block refinement levels; the graded ones all carry 2:1 seams
CASES = {
    "mixed":   [[3, 1, 2], [1, 2, 1], [2, 1, 3]],
    "peak":    [[0, 0, 1, 0, 0], [0, 1, 2, 1, 0], [1, 2, 3, 2, 1],
                [0, 1, 2, 1, 0], [0, 0, 1, 0, 0]],
    "uniform": 2,
}

TARGET = sys.argv[1] if len(sys.argv) > 1 else "oxley"
GROUP = sys.argv[2] if len(sys.argv) > 2 else "all"
CASE = sys.argv[3] if len(sys.argv) > 3 else "mixed"

if TARGET not in ("oxley", "finley"):
    raise ValueError("target must be oxley or finley, not %s" % TARGET)
if CASE not in CASES:
    raise ValueError("case must be one of %s" % sorted(CASES))

WORKDIR = os.environ.get("OXLEY_WORKDIR", ".")

# a case guaranteed to differ from CASE, for the dump tests
OTHER = "peak" if CASE != "peak" else "mixed"


def forest(case):
    """the oxley domain for a case, on the unit square"""
    levels = CASES[case]
    if isinstance(levels, int):
        n0 = n1 = 2
    else:
        n0, n1 = len(levels), len(levels[0])
    return Rectangle(n0=n0, n1=n1, l0=1., l1=1., refine_level=levels)


def domain(case):
    """the domain under test: the forest itself, or its finley export"""
    d = forest(case)
    return d.toFinley() if TARGET == "finley" else d


class Test_SpatialFunctions2D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """integration, interpolation between function spaces, gradients, normals"""
    def setUp(self):
        self.order = 1
        self.domain = domain(CASE)

    def tearDown(self):
        del self.order
        del self.domain


class Test_Utils2D(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    """the large algebraic suite over the domain's function spaces"""
    def setUp(self):
        self.domain = domain(CASE)
        # escript needs a live reference to the domain
        self.functionspace = FunctionOnBoundary(self.domain)
        self.workdir = WORKDIR

    def tearDown(self):
        del self.functionspace
        del self.domain


class Test_DomainInterface2D(Test_Domain):
    def setUp(self):
        self.boundary_tag_list = [1, 2, 10, 20]
        self.domain = domain(CASE)
        self.rdomain = domain("uniform")

    def tearDown(self):
        del self.domain
        del self.rdomain
        del self.boundary_tag_list


class Test_DataOps2D(Test_Dump, Test_SetDataPointValue, Test_Lazy):
    if TARGET == "oxley":
        @unittest.skip("known oxley gap: HDF5 dump/load of expanded data "
                       "reports insufficient sample ids (see run_escriptOnOxley)")
        def test_DumpAndLoad_Expanded(self):
            pass

    def setUp(self):
        self.domain = domain(CASE)
        # the dump tests need domains that disagree about sample counts with
        # the one under test - load() is required to REFUSE them - so this must
        # be a different case from CASE, whatever CASE is.
        self.domain_with_different_number_of_samples = domain(OTHER)
        self.domain_with_different_number_of_data_points_per_sample = domain(OTHER)
        self.domain_with_different_sample_ordering = domain(CASE)
        self.filename_base = WORKDIR
        self.mainfs = Function(self.domain)
        self.otherfs = Solution(self.domain)

    def tearDown(self):
        del self.domain
        del self.domain_with_different_number_of_samples
        del self.domain_with_different_number_of_data_points_per_sample
        del self.domain_with_different_sample_ordering
        del self.mainfs
        del self.otherfs


class Test_GradOnBoundary2D(unittest.TestCase):
    """
    Gradients on the boundary function spaces.

    The shared suite is Test_Util_Gradient_noBoundary, so nothing above covers
    these. They need covering separately because the boundary gradient combines
    the two values on the face with the element's other two corners to get the
    tangential derivative, and on the fine side of a seam one of those corners
    can be hanging.

    A globally linear field is the discriminator: it lies in the space exactly,
    so any departure from its constant gradient is the discretisation getting a
    corner value wrong, not approximation error.
    """
    def setUp(self):
        self.domain = domain(CASE)

    def tearDown(self):
        del self.domain

    # Not a property of the export: finley's plain face elements (Line2 here,
    # Tri3 in 3D) carry only the nodes ON the face, and two nodes cannot
    # determine a 2D gradient - the normal derivative is simply not there, and
    # grad() returns a wrong answer rather than refusing. Measured on finley's
    # OWN Rectangle: useElementsOnFace=True gives 1.3e-15, False gives exactly
    # 3.0, the dropped y-component. toFinley() builds the plain kind by choice
    # (2026-07-28), so this is the price of that choice, not a defect here.
    if TARGET == "finley":
        _skip_boundary = unittest.skip(
                "finley's plain face elements hold only the on-face nodes, so "
                "no boundary gradient exists; needs useElementsOnFace-style "
                "parent-shaped faces in the export")
    else:
        _skip_boundary = lambda f: f

    def _check(self, fs, name):
        x = ContinuousFunction(self.domain).getX()
        u = 2. * x[0] + 3. * x[1] - 1.
        g = grad(u, fs)
        err = Lsup(g - [2., 3.])
        self.assertLess(err, 1e-8, "%s: gradient off by %g" % (name, err))

    @_skip_boundary
    def test_grad_onFunctionOnBoundary(self):
        self._check(FunctionOnBoundary(self.domain), "FunctionOnBoundary")

    @_skip_boundary
    def test_grad_onReducedFunctionOnBoundary(self):
        self._check(ReducedFunctionOnBoundary(self.domain),
                    "ReducedFunctionOnBoundary")

    def test_grad_onFunction(self):
        self._check(Function(self.domain), "Function")

    def test_grad_onReducedFunction(self):
        self._check(ReducedFunction(self.domain), "ReducedFunction")


GROUPS = {
    "spatial": [Test_SpatialFunctions2D],
    "algebra": [Test_Utils2D],
    "objects": [Test_DomainInterface2D, Test_DataOps2D],
    "boundary": [Test_GradOnBoundary2D],
}
GROUPS["all"] = (GROUPS["spatial"] + GROUPS["algebra"] + GROUPS["objects"]
                 + GROUPS["boundary"])

if __name__ == '__main__':
    if GROUP not in GROUPS:
        raise ValueError("group must be one of %s" % sorted(GROUPS))
    if getMPIRankWorld() == 0:
        print("=== %s domain, case '%s', group '%s', %d rank(s) ==="
              % (TARGET, CASE, GROUP, getMPISizeWorld()))
        sys.stdout.flush()
    result = run_tests(__name__, GROUPS[GROUP])
    if getMPIRankWorld() == 0:
        print("\n%s/%s/%s: %d run, %d failures, %d errors, %d skipped"
              % (TARGET, CASE, GROUP, result.testsRun, len(result.failures),
                 len(result.errors), len(result.skipped)))
        sys.stdout.flush()
