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
The 2D oxley -> finley export of a graded (hanging-node) forest.

finley cannot represent a hanging node - one ReferenceElementSet per
ElementFile, and escript's q/r is pointwise Dirichlet rather than
u = (u_a+u_b)/2 - so a 2:1 seam is resolved by the TRIANGULATION: the hanging
position becomes an ordinary free node and the coarse element is split so
that it is a vertex on both sides. What comes out is a conforming P1 space on
a graded mesh.

Checks here, weakest to strongest:

  1. volume, boundary length and int(n.x)dS == dim*volume: geometry, and the
     winding of the face elements. The winding is also checked FACE BY FACE,
     because int(n.x)dS is a global sum in which a flipped face and a
     compensating one elsewhere cancel - and an inward normal silently
     reverses any Neumann or Robin term built on that face.
  2. a Poisson solve runs and stays bounded.
  3. the linear patch test: a linear field must be reproduced exactly, the
     standard check for a spurious constraint or kink.
  4. conformity - the one that matters. Every triangle edge must be shared by
     exactly two triangles or lie on the boundary; a quad split one way by one
     element and the other way by its neighbour leaves the mesh cracked, and
     NO physics test can see it, since a globally linear function lies in both
     triangulations and the patch test passes at 1e-16 on a cracked mesh.
     toFinley() checks this internally and throws, so getting a mesh back at
     all is that check passing.

Plus what does not travel with the mesh arrays and so was once dropped
silently: tag names, node tags and Dirac points.

Collective calls are made for every case BEFORE anything is asserted. An
assertion that fails on one rank only would otherwise leave that rank short
of a reduction the others still enter, and the run would hang somewhere far
from the cause.
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.oxley import Rectangle, toFinleyData, fromFinleyData

# REQUIRED: registers the concrete finley domain type with boost::python, or
# toFinley() hands back a base Domain with no getDescription
import esys.finley

# Cases chosen from a survey of what p4est_balance actually produces. Blocks
# are unit squares (l0=n0, l1=n1) so the elements stay square. Note the survey
# result: a big level jump does NOT give a cell many hanging edges - balance
# grades it gradually - and 3 or 4 hanging edges arise only from an ISOLATED
# COARSE CELL ringed by finer ones. Both extremes are here.
CASES = [
    # name, per-block levels, hanging edges per element seen in the survey
    ("one_seam",          [[2], [3]],                                 "1"),
    ("cascade",           [[3, 1], [1, 2]],                           "1,2"),
    ("mixed_3x3",         [[3,1,2],[1,2,1],[2,1,3]],                  "1,2,3"),
    ("isolated_coarse",   [[1,1,1],[1,0,1],[1,1,1]],                  "4"),
    ("big_jump",          [[0, 3], [3, 0]],                           "1,2"),
    ("diagonal_ramp",     [[0,1,2,3],[1,2,3,2],[2,3,2,1],[3,2,1,0]],  "1,2"),
    ("centre_peak",       [[0,0,1,0,0],[0,1,2,1,0],[1,2,3,2,1],
                           [0,1,2,1,0],[0,0,1,0,0]],                  "1,2"),
    ("checkerboard",      [[1,2,1,2],[2,1,2,1],[1,2,1,2],[2,1,2,1]],  "1,2"),
    ("corner_spike",      [[4,0,0],[0,0,0],[0,0,2]],                  "1,2"),
    ("stripes",           [[3,0,3],[0,3,0],[3,0,3],[0,3,0]],          "1,2"),
    ("uniform_control",   2,                                          "none"),
]


def blocks(levels):
    """the block counts, which are also the domain extent"""
    if isinstance(levels, int):
        return 2, 2
    return len(levels), len(levels[0])


def forest(levels, **kwargs):
    n0, n1 = blocks(levels)
    return Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                     refine_level=levels, **kwargs)


class Test_FinleyExportGeometry2D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # export once and share: toFinley() is the expensive part, and every
        # test below wants the same meshes
        cls.meshes = []
        for name, levels, hang in CASES:
            n0, n1 = blocks(levels)
            cls.meshes.append((name, forest(levels).toFinley(),
                               float(n0), float(n1)))

    @classmethod
    def tearDownClass(cls):
        del cls.meshes

    def test_volume(self):
        got = [(n, integrate(Scalar(1, Function(f))), lx * ly)
               for n, f, lx, ly in self.meshes]
        for name, vol, want in got:
            self.assertAlmostEqual(vol, want, 10, "volume of %s" % name)

    def test_boundary_length(self):
        got = [(n, integrate(Scalar(1, FunctionOnBoundary(f))), 2 * (lx + ly))
               for n, f, lx, ly in self.meshes]
        for name, length, want in got:
            self.assertAlmostEqual(length, want, 10,
                                   "boundary length of %s" % name)

    def test_face_winding_total(self):
        """int(n.x)dS == dim*volume"""
        got = []
        for name, f, lx, ly in self.meshes:
            fb = FunctionOnBoundary(f)
            got.append((name, integrate(inner(fb.getNormal(), fb.getX())),
                        2 * lx * ly))
        for name, val, want in got:
            self.assertAlmostEqual(val, want, 10, "int(n.x)dS on %s" % name)

    def test_face_winding_per_face(self):
        """every face normal points out, not just on average"""
        got = []
        for name, f, lx, ly in self.meshes:
            fb = FunctionOnBoundary(f)
            xb = fb.getX()
            outward = Vector(0, fb)
            outward[0] = whereZero(xb[0] - lx) - whereZero(xb[0])
            outward[1] = whereZero(xb[1] - ly) - whereZero(xb[1])
            got.append((name, inf(inner(fb.getNormal(), outward))))
        for name, worst in got:
            self.assertAlmostEqual(worst, 1.0, 10,
                                   "worst face normal on %s" % name)

    def test_boundary_tag_areas(self):
        got = []
        for name, f, lx, ly in self.meshes:
            for tag, want in (("left", ly), ("right", ly),
                              ("bottom", lx), ("top", lx)):
                mask = Scalar(0, FunctionOnBoundary(f))
                mask.setTaggedValue(tag, 1.0)
                got.append((name, tag, integrate(mask), want))
        for name, tag, area, want in got:
            self.assertAlmostEqual(area, want, 10,
                                   "tag '%s' area on %s" % (tag, name))

    def test_linear_patch(self):
        """a linear field must come back exactly: no spurious constraint"""
        got = []
        for name, f, lx, ly in self.meshes:
            x = ContinuousFunction(f).getX()
            exact = 1.0 + 2.0 * x[0] + 3.0 * x[1]
            pde = LinearPDE(f, numEquations=1)
            pde.setValue(A=kronecker(f))
            onbnd = (whereZero(x[0]) + whereZero(x[0] - lx)
                     + whereZero(x[1]) + whereZero(x[1] - ly))
            pde.setValue(q=onbnd, r=exact)
            pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
            pde.getSolverOptions().setTolerance(1e-12)
            u = pde.getSolution()
            got.append((name, Lsup(u - exact) / Lsup(exact)))
        for name, err in got:
            self.assertLess(err, 1e-9, "patch test on %s off by %g"
                            % (name, err))

    def test_poisson_runs(self):
        got = []
        for name, f, lx, ly in self.meshes:
            x = ContinuousFunction(f).getX()
            pde = LinearPDE(f, numEquations=1)
            pde.setSymmetryOn()
            pde.setValue(A=kronecker(f), Y=1.0,
                         q=whereZero(x[0]) + whereZero(x[0] - lx))
            pde.getSolverOptions().setTolerance(1e-12)
            got.append((name, sup(pde.getSolution())))
        for name, peak in got:
            self.assertTrue(peak == peak and abs(peak) < 1e30,
                            "Poisson solution on %s is not finite (%s)"
                            % (name, peak))


class Test_FinleyExportTagsAndDirac2D(unittest.TestCase):
    """
    What does NOT ride along with the mesh arrays, and so was once dropped
    without a word: tag names, node tags and Dirac points.
    """
    # asked-for points, deliberately off-node so that snapping is visible
    POINTS = [(0.3, 0.3), (0.7, 0.7), (0.1, 0.9)]
    TAGS = ["A", "B", "C"]
    LEVELS = [[3, 1, 2], [1, 2, 1], [2, 1, 3]]      # the 3x3 mixed forest

    def unitForest(self, **kwargs):
        """the graded forest on the UNIT square, so POINTS lie inside it"""
        n0, n1 = blocks(self.LEVELS)
        return Rectangle(n0=n0, n1=n1, l0=1., l1=1.,
                         refine_level=self.LEVELS, **kwargs)

    def hangingPosition(self, dom):
        """
        One hanging position, agreed by every rank.

        The adversarial place to ask for a Dirac point: it is NOT a node of the
        oxley domain, so oxley snaps to an lnodes node half an element away -
        but it IS a node of the export, at distance zero. A converter that
        located points on the exported mesh rather than reusing oxley's answer
        would put the point here instead, and nothing else would notice.

        Taken as the lexicographic maximum over ranks so all agree without a
        broadcast; fixed point because getMPIWorldMax only reduces ints, and
        the y only from ranks holding the winning x so the two coordinates
        cannot come from different points.
        """
        info = dom.getMeshInfo(True)
        xy = info["nodeCoords"].reshape(-1, 2)
        best = (-1, -1)
        for i in info["constrainedNodes"]:
            p = xy[int(i)]
            best = max(best, (int(round(p[0] * 1e6)), int(round(p[1] * 1e6))))
        gx = getMPIWorldMax(best[0])
        gy = getMPIWorldMax(best[1] if best[0] == gx else -1)
        return None if gx < 0 else (gx * 1e-6, gy * 1e-6)

    def moments(self, fs):
        """
        Count and coordinate moments of a Dirac space, summed over ranks.

        A fingerprint of the point SET that does not depend on which rank holds
        which point - ownership may legitimately differ between the domains,
        oxley assigning a point to the owner of the nearest node and finley to
        the owner of the DOF after prepare(). Fixed point again, since
        getMPIWorldSum reduces ints; the two domains hold the same doubles, so
        a 1e-6 grid is exact here rather than a tolerance.
        """
        x = fs.getX()
        s = [0] * 4
        for i in range(x.getNumberOfDataPoints()):
            p = x.getTupleForDataPoint(i)
            s[0] += int(round(p[0] * 1e6))
            s[1] += int(round(p[1] * 1e6))
            s[2] += int(round((p[0] * p[0] + p[1] * p[1]) * 1e6))
            s[3] += int(round(p[0] * p[1] * 1e6))
        return [getMPIWorldSum(x.getNumberOfDataPoints())] \
             + [getMPIWorldSum(v) for v in s]

    def test_tag_names_survive(self):
        dom = self.unitForest(diracPoints=self.POINTS, diracTags=self.TAGS)
        dom.setTagMap("myregion", 7)
        fin = dom.toFinley()
        names = [s.strip() for s in dom.showTagNames().split(",") if s.strip()]
        missing = [t for t in names if not fin.isValidTagName(t)]
        self.assertEqual(missing, [], "tag names lost by the export")
        wrong = [t for t in names if fin.getTag(t) != dom.getTag(t)]
        self.assertEqual(wrong, [], "tag values changed by the export")
        self.assertEqual(fin.getTag("myregion"), 7,
                         "a user-defined tag name must survive, not just the "
                         "boundary ones")

    def test_dirac_points_land_where_oxley_put_them(self):
        points, tags = list(self.POINTS), list(self.TAGS)
        hp = self.hangingPosition(self.unitForest())      # collective
        if hp is not None:
            points.append(hp)
            tags.append("H")
        dom = self.unitForest(diracPoints=points, diracTags=tags)
        fin = dom.toFinley()

        mo = self.moments(DiracDeltaFunctions(dom))
        mf = self.moments(DiracDeltaFunctions(fin))
        peaks = []
        for tag in tags:
            d = Data(0., DiracDeltaFunctions(fin))
            d.setTaggedValue(tag, 5.)
            peaks.append((tag, sup(d)))                   # collective

        self.assertEqual(mf[0], len(points), "the export lost Dirac points")
        self.assertEqual(mo[0], mf[0], "Dirac point count differs from oxley")
        self.assertEqual(mo[1:], mf[1:],
                         "Dirac points are not where oxley put them")
        for tag, peak in peaks:
            self.assertAlmostEqual(peak, 5., 10,
                                   "tag '%s' selects no point on the export"
                                   % tag)

    def test_node_tags_survive(self):
        dom = self.unitForest()
        cf = ContinuousFunction(dom)
        cf.setTags(5, whereNegative(cf.getX()[0] - 0.5))   # the left half
        fin = dom.toFinley()
        fcf = ContinuousFunction(fin)

        ind = Data(0., fcf)
        ind.setTaggedValue(5, 1.)
        stray = sup(ind * whereNonNegative(fcf.getX()[0] - 0.5))  # collective

        self.assertIn(5, fcf.getListOfTags(), "the node tag did not reach the "
                      "export")
        self.assertAlmostEqual(stray, 0., 10,
                               "a node outside the tagged region carries the "
                               "tag: the materialised nodes must inherit only "
                               "when their masters agree")



class Test_ContinuousFunctionTransfer2D(unittest.TestCase):
    """
    Carrying a ContinuousFunction between the forest and its export.

    Not an interpolation: the two meshes name their shared nodes identically,
    so the values are copied and the tests below expect EXACT equality, not
    closeness. The only nodes needing a rule are the ones the export has and
    the forest does not - the positions materialised at a 2:1 seam - which take
    the average of their masters.

    Works under MPI. finley redistributes the nodes when it prepares the
    domain, so rank r's finley nodes are not rank r's oxley nodes; the transfer
    passes a buffer around a ring over the global id range, which is finley's
    own idiom in NodeFile::gather_global and needs neither side to know the
    other's partition.
    """
    LEVELS = [[3, 1, 2], [1, 2, 1], [2, 1, 3]]

    def domains(self, levels=None):
        levels = self.LEVELS if levels is None else levels
        if isinstance(levels, int):
            n0 = n1 = 2
        else:
            n0, n1 = len(levels), len(levels[0])
        dom = Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                        refine_level=levels)
        return dom, dom.toFinley()

    def test_linear_field_is_exact_on_the_export(self):
        """
        A linear field lies in both spaces exactly, including at the seam
        positions - the average of two masters IS the value at their midpoint.
        So every node of the export must come out exact, not just the shared
        ones, and that is what checks the averaging rule.
        """
        for levels in (2, [[1], [2]], self.LEVELS):
            dom, fin = self.domains(levels)
            x = ContinuousFunction(dom).getX()
            u = 1. + 2. * x[0] + 3. * x[1]
            xf = ContinuousFunction(fin).getX()
            self.assertEqual(Lsup(toFinleyData(u, fin)
                                  - (1. + 2. * xf[0] + 3. * xf[1])), 0.,
                             "levels %s: the export is not exact" % (levels,))

    def test_round_trip_is_exact(self):
        for levels in (2, [[1], [2]], self.LEVELS):
            dom, fin = self.domains(levels)
            x = ContinuousFunction(dom).getX()
            for name, u in (("linear", 1. + 2. * x[0] + 3. * x[1]),
                            ("non-linear", sin(3. * x[0]) * cos(2. * x[1]))):
                back = fromFinleyData(toFinleyData(u, fin), dom)
                self.assertEqual(Lsup(back - u), 0.,
                                 "levels %s, %s field: the round trip lost "
                                 "something" % (levels, name))

    @unittest.skipIf(getMPISizeWorld() > 1,
                     "counts node samples, which double-count ghosts under MPI")
    def test_only_the_seam_nodes_are_averaged(self):
        """
        For a field that is NOT linear the seam values are averages rather than
        evaluations, so they must differ from the exact field - while every
        shared node still matches exactly. That the count of differing nodes is
        the count of materialised nodes is what says the averaging is confined
        to them.
        """
        dom, fin = self.domains()
        x = ContinuousFunction(dom).getX()
        u = sin(3. * x[0]) * cos(2. * x[1])
        xf = ContinuousFunction(fin).getX()
        diff = toFinleyData(u, fin) - (sin(3. * xf[0]) * cos(2. * xf[1]))

        nOx = u.getNumberOfDataPoints()
        nFin = diff.getNumberOfDataPoints()
        wrong = sum(1 for i in range(nFin)
                    if abs(diff.getTupleForDataPoint(i)[0]) > 1e-14)
        self.assertEqual(wrong, nFin - nOx,
                         "%d nodes differ from the exact field but only %d "
                         "exist solely on the export" % (wrong, nFin - nOx))

    def test_vector_data(self):
        dom, fin = self.domains()
        v = ContinuousFunction(dom).getX()
        self.assertEqual(Lsup(fromFinleyData(toFinleyData(v, fin), dom) - v), 0.)

    def test_wrong_function_space_is_refused(self):
        dom, fin = self.domains()
        self.assertRaises(RuntimeError, toFinleyData,
                          Data(1., Function(dom)), fin)

    def test_unrelated_domain_is_refused(self):
        """
        A finley mesh not built from this forest. The id ranges of two
        different forests overlap, so "every id I asked for was supplied" is
        NOT evidence of a match - the export carries a fingerprint of the
        forest it came from and that is what is checked.
        """
        dom, _ = self.domains()
        other = Rectangle(n0=2, n1=2, l0=2., l1=2., refine_level=1).toFinley()
        x = ContinuousFunction(dom).getX()
        self.assertRaises(RuntimeError, toFinleyData, x[0], other)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
