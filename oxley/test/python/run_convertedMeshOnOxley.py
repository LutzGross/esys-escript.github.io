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
from test_linearPDEs import Test_Poisson, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1

import oxley_meshes

try:
    import esys.finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False


def converted(levels=None):
    """the finley mesh exported from a graded oxley forest"""
    return oxley_meshes.graded(levels).toFinley()


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

    @unittest.skip("test_solve asks for 1e-5 relative on an exact solution "
                   "that is QUADRATIC, and finley passes it with second-order "
                   "elements - Rectangle(NE, NE, 2, useFullElementOrder=True) "
                   "in run_linearPDEsOnFinley3. The export is P1 on triangles, "
                   "so the error is discretisation error: measured 6.1e-3, "
                   "2.0e-3, 5.7e-4 as the forest is refined, clean O(h^2). "
                   "Reaching 1e-5 would take about 100k triangles.")
    def test_solve(self):
        pass

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


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
