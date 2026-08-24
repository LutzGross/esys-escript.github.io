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
escript's shared suites run on the mesh oxley EXPORTS to finley.

This is how the converter is tested as a MESH rather than as a set of arrays.
run_finleyExportOnOxley checks the export's geometry, its tags and the
transfer operators directly; here the exported mesh is simply handed to the
same suites any finley mesh would face, and has to behave like one.

What runs here is what depends on the mesh: integration, interpolation
between function spaces, gradients, normals, and PDE assembly and solves. The
algebraic util suite and the IO tests are deliberately absent - they are the
same code for any domain, so running them here would cost time and prove
nothing.

The source forests are GRADED. A uniform forest exports to a mesh with no 2:1
seam in it, so it would exercise none of the split.

KNOWN LIMIT OF THE EXPORT, not of this suite: toFinley() builds plain face
elements (Line2 in 2D), which carry only the nodes ON the face, so no
boundary gradient exists - grad(u, FunctionOnBoundary) returns the tangential
part alone. Measured on finley's OWN Rectangle: useElementsOnFace=True gives
1.3e-15, False gives exactly 3.0. Any test needing that is excluded by name
here, and the exclusion is the argument for switching the export to
parent-shaped faces.
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *

from test_util_spatial_functions1 import \
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from esys.escript.linearPDEs import Poisson
from test_linearPDEs import Test_Poisson, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1

import oxley_meshes

try:
    import esys.finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False


def converted(levels=None):
    """the finley mesh exported from a graded oxley forest"""
    return oxley_meshes.graded(levels).toFinley()


def converted3D(levels=None):
    """
    The same in 3D, where the split has more to do: a 2:1 seam puts a node at
    the centre of a coarse FACE and at the midpoints of its EDGES, and an edge
    can hang on its own where only a diagonal neighbour is finer. A coarse
    octant is then coned from its own centre - a node the exported mesh has and
    the forest does not - into up to forty-eight tetrahedra.
    """
    return oxley_meshes.graded3D(levels).toFinley()


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_SpatialFunctionsOnConvertedMesh(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """
    integration, interpolation between spaces, gradients and normals.

    The export resolves each 2:1 seam by triangulating it, so unlike the
    forest it came from it is a conforming P1 space - this is what fails if
    the split, the node numbering or the export's node ownership is wrong.
    """
    def setUp(self):
        self.order = 1
        self.domain = converted()

    def tearDown(self):
        del self.order
        del self.domain


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_PoissonOnConvertedMesh(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8

    def test_solve(self):
        """
        Replaces the shared test rather than skipping it.

        The shared version asks for 1e-5 relative on an exact solution that is
        QUADRATIC, and finley passes it with SECOND-ORDER elements -
        Rectangle(NE, NE, 2, useFullElementOrder=True) in
        run_linearPDEsOnFinley3. The export is P1 on triangles, so that
        tolerance is unreachable: what is left is discretisation error.

        What can be asked of a P1 mesh is that the error behave like one, so
        this solves the same problem on two refinements and checks the error
        falls at the expected rate. That is a statement about the exported
        discretisation; an absolute tolerance would only have been a statement
        about how fine the mesh happened to be.
        """
        def solveOn(levels):
            dom = oxley_meshes.graded(levels).toFinley()
            cf = ContinuousFunction(dom)
            x = cf.getX()
            u_ex = Scalar(1., cf)
            for i in range(dom.getDim()):
                u_ex *= x[i] * (2. - x[i])
            msk = Scalar(0., cf)
            for i in range(dom.getDim()):
                msk += whereZero(x[i])
            f = Scalar(0, cf)
            for i in range(dom.getDim()):
                f_i = Scalar(2., cf)
                for j in range(dom.getDim()):
                    if i != j:
                        f_i *= x[j] * (2. - x[j])
                f += f_i
            pde = Poisson(dom, debug=self.DEBUG)
            pde.setValue(f=f, q=msk)
            pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
            u = pde.getSolution()
            return Lsup(u - u_ex) / Lsup(u_ex)

        coarse = solveOn(oxley_meshes.MIXED_2D)
        fine = solveOn([[l + 1 for l in row] for row in oxley_meshes.MIXED_2D])

        self.assertLess(coarse, 5.e-2,
                        "the coarse solution is not even close: %g" % coarse)
        self.assertGreater(coarse / fine, 2.5,
                           "the error does not fall like a P1 discretisation: "
                           "%g -> %g is a factor of %g, expected about 4"
                           % (coarse, fine, coarse / fine))

    def setUp(self):
        self.domain = converted()

    def tearDown(self):
        del self.domain


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_LinearPDEOnConvertedMesh(Test_LinearPDE_noLumping,
                                    Test_assemblage_2Do1):
    """
    PDE assembly and solves on the exported mesh.

    The export is P1 on triangles, so the order-1 assemblage tests are the
    matching ones; anything expecting order 2 does not apply.
    """
    RES_TOL=1.e-7
    ABS_TOL=1.e-8

    def setUp(self):
        self.domain = converted()
        self.order = 1

    def tearDown(self):
        del self.domain


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_SpatialFunctionsOnConvertedMesh3D(
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    """
    The 3D export, on the same shared suite. What it adds over the 2D case is
    the tetrahedral split of a hanging octant: faces that are polygons rather
    than quads, an interior apex, and a boundary quad that becomes up to six
    triangles rather than two.
    """
    def setUp(self):
        self.order = 1
        self.domain = converted3D()

    def tearDown(self):
        del self.order
        del self.domain


@unittest.skipIf(not HAVE_FINLEY, "finley not available")
class Test_LinearPDEOnConvertedMesh3D(Test_LinearPDE_noLumping,
                                      Test_assemblage_3Do1):
    """
    PDE assembly and solves on the 3D export - P1 on tetrahedra, so the
    order-1 assemblage tests are the matching ones.
    """
    RES_TOL=1.e-7
    ABS_TOL=1.e-8

    def setUp(self):
        self.domain = converted3D()
        self.order = 1

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
