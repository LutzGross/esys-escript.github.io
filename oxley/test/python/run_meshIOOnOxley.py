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
saveMesh / loadMesh.

loadMesh is a function, not a method, for the reason the domain has no refine
methods either: reading into an existing domain replaces its mesh, and every
Data already built over that domain is then the wrong size with nothing in the
types to say so. It returns a NEW domain instead.

That also means a saved mesh has to carry enough to rebuild a domain from
nothing. p4est's own files hold the connectivity and the quadrants but not
where the domain sits in space, so saveMesh writes that alongside as
<filename>.oxley; the geometry checks below are what would catch it going
missing.
"""

import glob
import os

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick, RefinementQueue2D, loadMesh

try:
    WORKDIR = os.environ['OXLEY_WORKDIR']
except KeyError:
    WORKDIR = '.'


def numElements(domain):
    """
    GLOBAL element count. The local one is no use here: saveMesh writes
    without the partition and loadMesh re-partitions, so the same mesh
    legitimately splits differently across ranks after a round trip.
    """
    perElement = 4 if domain.getDim() == 2 else 8
    local = Data(1., Function(domain)).getNumberOfDataPoints() // perElement
    return getMPIWorldSum(int(local))


class Test_MeshIOOnOxley(unittest.TestCase):
    def path(self, name):
        return os.path.join(WORKDIR, "_meshio_%s_%d" % (name, getMPISizeWorld()))

    def tearDown(self):
        """
        saveMesh writes several files per mesh - the header, the connectivity
        and the forest - and nothing else removes them, so a run used to leave
        them lying in whatever directory it started from. Cleared here rather
        than at the end of the module so that a test which fails still tidies
        up after itself.

        Rank 0 alone removes, behind a barrier: every rank writes the same
        paths, so removing on all of them is a race in which the losers find
        the file already gone.
        """
        MPIBarrierWorld()
        if getMPIRankWorld() == 0:
            for f in glob.glob(os.path.join(WORKDIR, "_meshio_*")):
                try:
                    os.remove(f)
                except OSError:
                    pass
        MPIBarrierWorld()

    def check_same_mesh(self, a, b, what):
        self.assertEqual(a.getDim(), b.getDim(), "%s: dimension" % what)
        self.assertEqual(numElements(a), numElements(b), "%s: element count" % what)
        self.assertAlmostEqual(integrate(Scalar(1., Function(a))),
                               integrate(Scalar(1., Function(b))), 10,
                               "%s: volume" % what)
        xa = ContinuousFunction(a).getX()
        xb = ContinuousFunction(b).getX()
        for i in range(a.getDim()):
            self.assertAlmostEqual(inf(xa[i]), inf(xb[i]), 10,
                                   "%s: lower bound of axis %d" % (what, i))
            self.assertAlmostEqual(sup(xa[i]), sup(xb[i]), 10,
                                   "%s: upper bound of axis %d" % (what, i))

    def test_roundtrip_2D(self):
        """graded, off-unit extent: the geometry has to survive too"""
        dom = Rectangle(n0=3, n1=2, l0=3., l1=2.,
                        refine_level=[[2, 1], [1, 2], [2, 0]])
        f = self.path("2d")
        dom.saveMesh(f)
        self.check_same_mesh(dom, loadMesh(f), "2D round trip")

    def test_roundtrip_2D_after_refinement(self):
        dom = Rectangle(n0=2, n1=2, l0=1., l1=1.)
        q = RefinementQueue2D()
        q.setRefinementLevel(3)
        q.refineRegion(x0=0., y0=0., x1=0.4, y1=0.4)
        dom = q.apply(dom)
        f = self.path("2dq")
        dom.saveMesh(f)
        self.check_same_mesh(dom, loadMesh(f), "2D refined round trip")

    def test_roundtrip_3D(self):
        dom = Brick(n0=2, n1=2, n2=2, l0=2., l1=2., l2=2., refine_level=1)
        f = self.path("3d")
        dom.saveMesh(f)
        self.check_same_mesh(dom, loadMesh(f), "3D round trip")

    def test_source_domain_is_untouched(self):
        """
        The point of loadMesh being a function. With the old method, reading a
        mesh into a domain left s below sized for a mesh that no longer
        existed.
        """
        dom = Rectangle(n0=2, n1=2, l0=1., l1=1., refine_level=2)
        s = Scalar(1., ContinuousFunction(dom))
        before = numElements(dom)

        other = Rectangle(n0=2, n1=2, l0=1., l1=1., refine_level=1)
        f = self.path("other")
        other.saveMesh(f)
        loaded = loadMesh(f)

        self.assertEqual(numElements(dom), before, "the domain changed")
        self.assertEqual(Lsup(s + s), 2.0, "Data over the domain broke")
        self.assertNotEqual(numElements(loaded), before,
                            "this test needs the two meshes to differ")

    def test_domain_has_no_loadMesh(self):
        dom = Rectangle(n0=2, n1=2, l0=1., l1=1.)
        self.assertFalse(hasattr(dom, "loadMesh"),
                         "loadMesh must not be a method on the domain: it "
                         "would replace the mesh under any Data already "
                         "defined over it")

    def test_missing_file_is_refused(self):
        self.assertRaises(RuntimeError, loadMesh,
                          os.path.join(WORKDIR, "_meshio_no_such_mesh"))


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
