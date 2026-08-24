
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_util import Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
# from test_util_interpolation import Test_Util_Point_Data_Interpolation
from test_symfuncs import Test_symfuncs
from esys.escript import *
from esys.oxley import Rectangle, Brick

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    print("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

from test_util_NaN_funcs import Test_util_NaN_funcs

import oxley_meshes

# The meshes come from oxley_meshes: graded, because a uniform forest never
# enters the hanging-node path. There is no NX/NY here any more - it only fed
# d0/d1, which oxley ignores with a warning since p4est owns the partition.

class Test_UtilOnOxley(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=oxley_meshes.graded()
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """integration, interpolation and gradients on a graded mesh"""
    def setUp(self):
        self.order=1
        self.domain = oxley_meshes.graded()
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnUniformOxley2D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """
    The control. It must pass whatever the seam code does, so when the graded
    class above fails and this one does not, the fault is in the hanging-node
    path rather than in oxley generally - which is what made the O(1) error in
    grad() diagnosable in minutes instead of days.
    """
    def setUp(self):
        self.order=1
        self.domain = oxley_meshes.uniform()
    def tearDown(self):
        del self.order
        del self.domain


class Test_Util_SpatialFunctionsOnOxley3D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """
    3D, the uniform control: it must pass whatever the seam code does, so that
    a failure of the graded class below points at the hanging-node path.
    """
    def setUp(self):
        self.order=1
        self.domain = oxley_meshes.uniform3D()
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnGradedOxley3D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """the 3D counterpart of the graded 2D class - the seams are the point"""
    def setUp(self):
        self.order=1
        self.domain = oxley_meshes.graded3D()
    def tearDown(self):
        del self.order
        del self.domain

#TODO
# class Test_2D_Point_Data_Integration(Test_Util_Point_Data_Interpolation):
#     def setUp(self):
#         Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
#         StationsTags = ["A1", "A2", "A3", "A4" ]
#         self.domain=Rectangle(n0=5,n1=5, diracPoints=Stations, diracTags=StationsTags)
#     def tearDown(self):
#         del self.domain

# Todo
# class Test_3D_Point_Data_Integration(Test_Util_Point_Data_Interpolation):
#     def setUp(self):
#         Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.) ]
#         StationsTags = ["A1", "A2", "A3", "A4" ]
#         self.domain=Brick(n0=5,n1=5,n2=5,diracPoints=Stations,diracTags=StationsTags)
#     def tearDown(self):
#         del self.domain


class Test_GradientOnBoundaryOnOxley2D(unittest.TestCase):
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
        dom = (oxley_meshes.uniform(levels) if isinstance(levels, int)
               else oxley_meshes.graded(levels))
        for fs, name in ((Function(dom), "Function"),
                         (ReducedFunction(dom), "ReducedFunction"),
                         (FunctionOnBoundary(dom), "FunctionOnBoundary"),
                         (ReducedFunctionOnBoundary(dom),
                          "ReducedFunctionOnBoundary")):
            self._check(dom, fs, "%s on %s" % (name, label))

    def test_gradient_mixed(self):
        self._allSpaces(oxley_meshes.MIXED_2D, "3x3 mixed")

    def test_gradient_peak(self):
        self._allSpaces(oxley_meshes.PEAK_2D, "centre peak")

    def test_gradient_uniform(self):
        self._allSpaces(oxley_meshes.UNIFORM, "uniform (control)")


class Test_GradientOnOxley3D(unittest.TestCase):
    """
    Gradients on a graded 3D forest, in every function space.

    The 3D counterpart of the class above, and it needs its own because 3D has
    a kind of hanging node 2D does not: one in the middle of a coarse FACE,
    whose four masters include two that hang in turn on that face's edges. The
    isolated case is the one where those chains meet.

    Two fields. A LINEAR one has a constant gradient, so any departure means a
    corner value was read from the wrong place. A TRILINEAR one also lies in
    the space exactly but its gradient VARIES inside an element, so it pins the
    quadrature points themselves - a plain permutation of them passes the
    linear test and fails this one.
    """
    def _check(self, dom, fs, label):
        X = ContinuousFunction(dom).getX()
        err = Lsup(grad(2.*X[0] + 3.*X[1] - X[2] + 5., fs) - [2., 3., -1.])
        self.assertLess(err, 1e-8, "%s: linear gradient off by %g" % (label, err))

        x = fs.getX()
        g = grad(X[0]*X[1]*X[2], fs)
        exact = g * 0.
        exact[0] = x[1]*x[2]
        exact[1] = x[0]*x[2]
        exact[2] = x[0]*x[1]
        err = Lsup(g - exact)
        self.assertLess(err, 1e-8,
                        "%s: trilinear gradient off by %g" % (label, err))

    def _allSpaces(self, levels, label):
        dom = (oxley_meshes.uniform3D(levels) if isinstance(levels, int)
               else oxley_meshes.graded3D(levels))
        for fs, name in ((Function(dom), "Function"),
                         (ReducedFunction(dom), "ReducedFunction"),
                         (FunctionOnBoundary(dom), "FunctionOnBoundary"),
                         (ReducedFunctionOnBoundary(dom),
                          "ReducedFunctionOnBoundary")):
            self._check(dom, fs, "%s on %s" % (name, label))

    def test_gradient_mixed(self):
        self._allSpaces(oxley_meshes.MIXED_3D, "2x2x2 mixed")

    def test_gradient_seam(self):
        self._allSpaces(oxley_meshes.SEAM_3D, "single seam")

    def test_gradient_isolated(self):
        self._allSpaces(oxley_meshes.ISOLATED_3D, "isolated coarse cell")

    def test_gradient_uniform(self):
        self._allSpaces(1, "uniform (control)")


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
