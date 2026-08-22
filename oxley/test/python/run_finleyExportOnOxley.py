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
from esys.oxley import (Rectangle, Brick, toFinleyData, fromFinleyData,
                       toFinleyReducedData, fromFinleyReducedData,
                       toFinleyBoundaryData, fromFinleyBoundaryData,
                       toFinleyFunctionData, fromFinleyFunctionData)

# REQUIRED: registers the concrete finley domain type with boost::python, or
# toFinley() hands back a base Domain with no getDescription
import oxley_meshes
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



class Test_ReducedFunctionTransfer2D(unittest.TestCase):
    """
    Carrying a ReducedFunction - one value per element - across.

    An octant becomes 2 to 6 triangles, so this one is not a copy: outbound the
    value is REPLICATED onto each simplex, inbound it is the AREA-WEIGHTED mean
    of them. The weights are recomputed on the oxley side from the split, which
    is a deterministic function of the hanging configuration, so the finley
    side never has to send areas.

    No map is stored either: the exported simplices carry ids that say which
    octant they came from, so the correspondence survives finley redistributing
    the mesh in prepare().

    Two invariants are worth more than the errors here. Replicate-then-average
    is the identity, because the weights of an octant sum to one. And the
    INTEGRAL is preserved in both directions, which is the property an error
    indicator needs - it is a density, and moving it between meshes must not
    create or destroy any of it.
    """
    CASES = [("conforming", 2), ("one_seam", [[1], [2]]),
             ("mixed_3x3", [[3, 1, 2], [1, 2, 1], [2, 1, 3]])]

    def domains(self, levels):
        n0, n1 = blocks(levels)
        dom = Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                        refine_level=levels)
        return dom, dom.toFinley()

    def test_replicate_then_average_is_the_identity(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x = ReducedFunction(dom).getX()
            u = 1. + x[0] * x[1]
            back = fromFinleyReducedData(toFinleyReducedData(u, fin), dom)
            self.assertEqual(Lsup(back - u), 0.,
                             "%s: the round trip changed the field" % name)

    def test_integral_is_preserved_outbound(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x = ReducedFunction(dom).getX()
            u = 1. + x[0] * x[1]
            a, b = integrate(u), integrate(toFinleyReducedData(u, fin))
            self.assertAlmostEqual(a, b, 10, "%s: integral changed on export "
                                   "(%.12g -> %.12g)" % (name, a, b))

    def test_integral_is_preserved_inbound(self):
        """
        The direction the adaptive loop needs: an indicator computed on the
        finley mesh coming home per octant.
        """
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            xf = ReducedFunction(fin).getX()
            ind = xf[0] * xf[0] + xf[1]
            home = fromFinleyReducedData(ind, dom)
            a, b = integrate(ind), integrate(home)
            self.assertAlmostEqual(a, b, 10, "%s: integral changed coming home "
                                   "(%.12g -> %.12g)" % (name, a, b))

    def test_constant_survives_any_split(self):
        """
        the weights of an octant sum to one, whatever pattern it was split by
        """
        dom, fin = self.domains(self.CASES[2][1])
        home = fromFinleyReducedData(Data(2.5, ReducedFunction(fin)), dom)
        self.assertEqual(Lsup(home - 2.5), 0.)

    def test_wrong_function_space_is_refused(self):
        dom, fin = self.domains(2)
        self.assertRaises(RuntimeError, toFinleyReducedData,
                          Data(1., ContinuousFunction(dom)), fin)

    def test_unrelated_domain_is_refused(self):
        dom, _ = self.domains(2)
        other = Rectangle(n0=3, n1=3, l0=3., l1=3., refine_level=1).toFinley()
        self.assertRaises(RuntimeError, toFinleyReducedData,
                          Data(1., ReducedFunction(dom)), other)



class Test_BoundaryTransfer2D(unittest.TestCase):
    """
    Carrying FunctionOnBoundary, and its reduced form, across.

    The boundary faces correspond one to one: the 2D split never subdivides a
    boundary edge, because a hanging node is the midpoint of a face that HAS a
    finer neighbour and is therefore interior. Measured: both meshes put the
    same physical quadrature points on the boundary.

    They do not agree on the ORDER though - neither of the faces nor of the
    points within a face - so this is a permutation. Faces are matched by an id
    that says which face of which octant they are; points within a face are
    matched by sorting each face's points by coordinate, so neither side has to
    send coordinates.

    A field that VARIES along the boundary is what makes a wrong permutation
    visible: a constant would survive any mismatch, and so would the integral.
    """
    CASES = [("conforming", 2), ("one_seam", [[1], [2]]),
             ("mixed_3x3", [[3, 1, 2], [1, 2, 1], [2, 1, 3]])]
    SPACES = [("FunctionOnBoundary", FunctionOnBoundary),
              ("ReducedFunctionOnBoundary", ReducedFunctionOnBoundary)]

    def domains(self, levels):
        n0, n1 = blocks(levels)
        dom = Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                        refine_level=levels)
        return dom, dom.toFinley()

    def test_varying_field_lands_on_the_right_faces(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            for label, fs in self.SPACES:
                x = fs(dom).getX()
                u = x[0] * x[0] + 3. * x[1]
                xf = fs(fin).getX()
                err = Lsup(toFinleyBoundaryData(u, fin)
                           - (xf[0] * xf[0] + 3. * xf[1]))
                self.assertLess(err, 1e-12,
                                "%s/%s: values landed on the wrong faces "
                                "(err %g)" % (name, label, err))

    def test_round_trip_is_exact(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            for label, fs in self.SPACES:
                x = fs(dom).getX()
                u = x[0] * x[0] + 3. * x[1]
                back = fromFinleyBoundaryData(toFinleyBoundaryData(u, fin), dom)
                self.assertEqual(Lsup(back - u), 0.,
                                 "%s/%s: the round trip changed the field"
                                 % (name, label))

    def test_surface_integral_is_preserved(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            for label, fs in self.SPACES:
                x = fs(dom).getX()
                u = x[0] * x[0] + 3. * x[1]
                a, b = integrate(u), integrate(toFinleyBoundaryData(u, fin))
                self.assertAlmostEqual(a, b, 10, "%s/%s: surface integral "
                                       "changed (%.12g -> %.12g)"
                                       % (name, label, a, b))

    def test_the_normal_survives(self):
        """
        the sharpest check available: a wrong face or a flipped edge shows up
        immediately, since the normal differs between neighbouring faces
        """
        dom, fin = self.domains(self.CASES[2][1])
        moved = toFinleyBoundaryData(FunctionOnBoundary(dom).getNormal(), fin)
        self.assertEqual(Lsup(moved - FunctionOnBoundary(fin).getNormal()), 0.)

    def test_wrong_function_space_is_refused(self):
        dom, fin = self.domains(2)
        self.assertRaises(RuntimeError, toFinleyBoundaryData,
                          Data(1., ContinuousFunction(dom)), fin)

    def test_unrelated_domain_is_refused(self):
        dom, _ = self.domains(2)
        other = Rectangle(n0=3, n1=3, l0=3., l1=3., refine_level=1).toFinley()
        self.assertRaises(RuntimeError, toFinleyBoundaryData,
                          Data(1., FunctionOnBoundary(dom)), other)



class Test_FunctionTransfer2D(unittest.TestCase):
    """
    Carrying Function - values at the quadrature points - across.

    The one transfer that is not a rearrangement: an octant carries the 2x2
    Gauss points, a Tri3 its three edge midpoints, and neither set contains
    the other. So values are EVALUATED, not moved, and what can be asked of it
    is exactness on the fields each side can represent:

      outbound  the four values of an octant are unisolvent for a bilinear
                function, so anything bilinear - every linear field included -
                crosses exactly.
      inbound   three edge midpoints are unisolvent for a linear function on
                the triangle, so a field linear on the split comes home
                exactly. That is what the export's own P1 space produces.

    The integral is preserved exactly OUTBOUND, and that is provable rather
    than lucky: 2x2 Gauss is exact for the bilinear interpolant, xy has total
    degree two, and the three-midpoint rule on a triangle is exact to degree
    two. Inbound it is not, since the octant's rule then samples a function
    that is only piecewise linear.
    """
    CASES = [("conforming", 2), ("one_seam", [[1], [2]]),
             ("mixed_3x3", [[3, 1, 2], [1, 2, 1], [2, 1, 3]])]
    TOL = 1e-12

    def domains(self, levels):
        n0, n1 = blocks(levels)
        dom = Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                        refine_level=levels)
        return dom, dom.toFinley()

    def test_linear_is_exact_outbound(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x, xf = Function(dom).getX(), Function(fin).getX()
            err = Lsup(toFinleyFunctionData(1. + 2.*x[0] + 3.*x[1], fin)
                       - (1. + 2.*xf[0] + 3.*xf[1]))
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_bilinear_is_exact_outbound(self):
        """the sharper claim: the octant's four values fix a bilinear field"""
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x, xf = Function(dom).getX(), Function(fin).getX()
            err = Lsup(toFinleyFunctionData(x[0]*x[1], fin) - xf[0]*xf[1])
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_linear_is_exact_inbound(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x, xf = Function(dom).getX(), Function(fin).getX()
            err = Lsup(fromFinleyFunctionData(1. + 2.*xf[0] + 3.*xf[1], dom)
                       - (1. + 2.*x[0] + 3.*x[1]))
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_round_trip_is_exact_for_a_linear_field(self):
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x = Function(dom).getX()
            u = 1. + 2.*x[0] + 3.*x[1]
            err = Lsup(fromFinleyFunctionData(toFinleyFunctionData(u, fin), dom)
                       - u)
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_integral_is_preserved_outbound(self):
        """holds for ANY field, not just the ones that cross exactly"""
        for name, levels in self.CASES:
            dom, fin = self.domains(levels)
            x = Function(dom).getX()
            w = sin(2.*x[0]) * cos(x[1])
            a, b = integrate(w), integrate(toFinleyFunctionData(w, fin))
            self.assertAlmostEqual(a, b, 10, "%s: %.12g -> %.12g"
                                   % (name, a, b))

    def test_wrong_function_space_is_refused(self):
        dom, fin = self.domains(2)
        self.assertRaises(RuntimeError, toFinleyFunctionData,
                          Data(1., ReducedFunction(dom)), fin)

    def test_unrelated_domain_is_refused(self):
        dom, _ = self.domains(2)
        other = Rectangle(n0=3, n1=3, l0=3., l1=3., refine_level=1).toFinley()
        self.assertRaises(RuntimeError, toFinleyFunctionData,
                          Data(1., Function(dom)), other)



class Test_ComplexTransfer2D(unittest.TestCase):
    """
    The same four transfers on complex data.

    There is no separate machinery: std::complex is two doubles in memory and
    every operation these transfers perform - copying a value, averaging
    masters, weighting by area, evaluating a shape function - is REAL-LINEAR,
    so it applies to the two parts independently. A complex field is a real one
    with twice as many components.

    Which is exactly why it needs testing: the real and imaginary parts must
    not be mixed or dropped, so the two parts here are different functions of
    position, and a field that lost its imaginary part - or copied the real one
    into it - would fail.
    """
    LEVELS = [[3, 1, 2], [1, 2, 1], [2, 1, 3]]

    def setUp(self):
        n0, n1 = blocks(self.LEVELS)
        self.dom = Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                             refine_level=self.LEVELS)
        self.fin = self.dom.toFinley()

    def tearDown(self):
        del self.dom
        del self.fin

    def field(self, fs, domain):
        x = fs(domain).getX()
        return (1. + 2.*x[0] + 3.*x[1]) + 1j * (0.5 - x[0] + 4.*x[1])

    def test_continuous_function(self):
        u = self.field(ContinuousFunction, self.dom)
        uf = toFinleyData(u, self.fin)
        self.assertTrue(uf.isComplex(), "the result lost its complexity")
        self.assertEqual(Lsup(uf - self.field(ContinuousFunction, self.fin)), 0.)
        self.assertEqual(Lsup(fromFinleyData(uf, self.dom) - u), 0.)

    def test_reduced_function(self):
        u = self.field(ReducedFunction, self.dom)
        uf = toFinleyReducedData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertEqual(Lsup(fromFinleyReducedData(uf, self.dom) - u), 0.)

    def test_function_on_boundary(self):
        u = self.field(FunctionOnBoundary, self.dom)
        uf = toFinleyBoundaryData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertLess(Lsup(uf - self.field(FunctionOnBoundary, self.fin)), 1e-12)
        self.assertEqual(Lsup(fromFinleyBoundaryData(uf, self.dom) - u), 0.)

    def test_function(self):
        u = self.field(Function, self.dom)
        uf = toFinleyFunctionData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertLess(Lsup(uf - self.field(Function, self.fin)), 1e-12)
        self.assertLess(Lsup(fromFinleyFunctionData(uf, self.dom) - u), 1e-12)

    def test_a_real_field_stays_real(self):
        r = Data(1., ContinuousFunction(self.dom))
        self.assertFalse(toFinleyData(r, self.fin).isComplex(),
                         "a real field must not come back complex")



# ---------------------------------------------------------------------------
# 3D. The forests are CONFORMING, because toFinley() still refuses a graded 3D
# one - the tetrahedral split has no hanging cases yet. What is being tested is
# the transfer, and it changes shape in 3D even without a seam:
#
#   - an octant becomes SIX tetrahedra and a boundary quad TWO triangles, so
#     only ContinuousFunction is still a copy. Everything else evaluates.
#   - the two sides' quadrature rules are unrelated: 2x2x2 Gauss on an octant
#     against a four-point rule on a Tet4, 2x2 on a boundary quad against three
#     edge midpoints on a Tri3. So values travel with the COORDINATES of the
#     points they were taken at and the receiver rebuilds the polynomial those
#     points determine; neither side has to know the other's rule.
#
# The exactness claims follow from what each set of points is unisolvent for,
# and that is what these tests check - not a tolerance.
# ---------------------------------------------------------------------------

CASES_3D = [
    ("unit_cube", dict(n0=2, n1=2, n2=2, l0=1., l1=1., l2=1., refine_level=1)),
    ("refined", dict(n0=1, n1=1, n2=1, l0=1., l1=1., l2=1., refine_level=2)),
    # not a cube, not the same number of blocks per axis, so nothing can hide
    # behind a symmetry of the split
    ("oblong", dict(n0=3, n1=2, n2=1, l0=1.5, l1=2., l2=3., refine_level=1)),
]

# The GRADED forests, from the shared table, are what the tetrahedral split has
# to cope with. A 2:1 seam in 3D puts a node at the centre of a coarse FACE and
# at the midpoints of its four EDGES, and an edge can hang on its own where only
# a diagonal neighbour is finer - so a coarse octant's faces are polygons of up
# to eight vertices, and it is coned from its own CENTRE rather than a corner,
# no corner of it serving once a face hangs.
#
# Everything asserted below holds on these too, with one exception that is its
# own class: a field the export cannot represent is no longer copied exactly,
# because the export has nodes the forest does not.
GRADED_3D = [("seam", oxley_meshes.SEAM_3D),
             ("mixed", oxley_meshes.MIXED_3D),
             ("isolated", oxley_meshes.ISOLATED_3D)]

CONFORMING_FORESTS_3D = [(n, (lambda kw=kw: Brick(**kw))) for n, kw in CASES_3D]
GRADED_FORESTS_3D = [(n, (lambda l=l: oxley_meshes.graded3D(l)))
                     for n, l in GRADED_3D]
FORESTS_3D = CONFORMING_FORESTS_3D + GRADED_FORESTS_3D


class Test_ContinuousFunctionTransfer3D(unittest.TestCase):
    """
    Still a copy in 3D: the export hands finley the forest's own node ids, so
    the two meshes name the same nodes, and a conforming forest has no
    materialised position needing a rule. Every field, however nonlinear, must
    therefore cross EXACTLY - which is a sharper statement than in 2D, where
    the seam nodes are averages.
    """
    def test_any_field_is_exact_both_ways(self):
        for name, make in CONFORMING_FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x, xf = ContinuousFunction(dom).getX(), ContinuousFunction(fin).getX()
            u = sin(3.*x[0]) * cos(2.*x[1]) * exp(x[2])
            uf = sin(3.*xf[0]) * cos(2.*xf[1]) * exp(xf[2])
            self.assertEqual(Lsup(toFinleyData(u, fin) - uf), 0.,
                             "%s: the export is not exact" % name)
            self.assertEqual(Lsup(fromFinleyData(uf, dom) - u), 0.,
                             "%s: the way home is not exact" % name)

    def test_vector_data(self):
        dom = CONFORMING_FORESTS_3D[0][1]()
        fin = dom.toFinley()
        v = ContinuousFunction(dom).getX()
        self.assertEqual(Lsup(fromFinleyData(toFinleyData(v, fin), dom) - v), 0.)


class Test_ReducedFunctionTransfer3D(unittest.TestCase):
    """
    One value per octant, replicated onto the six tets and volume-averaged on
    the way back. The weights are recomputed on the oxley side from the same
    cone split the export emitted, so the finley side never sends volumes.
    """
    def test_replicate_then_average_is_the_identity(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = ReducedFunction(dom).getX()
            u = 1. + x[0]*x[1]*x[2]
            back = fromFinleyReducedData(toFinleyReducedData(u, fin), dom)
            self.assertLess(Lsup(back - u), 1e-14,
                            "%s: the round trip changed the field" % name)

    def test_integral_is_preserved_both_ways(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = ReducedFunction(dom).getX()
            u = 1. + x[0]*x[1]*x[2]
            a, b = integrate(u), integrate(toFinleyReducedData(u, fin))
            self.assertAlmostEqual(a, b, 10, "%s: integral changed on export "
                                   "(%.12g -> %.12g)" % (name, a, b))
            xf = ReducedFunction(fin).getX()
            ind = xf[0]*xf[0] + xf[1]
            c, d = integrate(ind), integrate(fromFinleyReducedData(ind, dom))
            self.assertAlmostEqual(c, d, 10, "%s: integral changed coming home "
                                   "(%.12g -> %.12g)" % (name, c, d))

    def test_constant_survives_the_split(self):
        dom = CONFORMING_FORESTS_3D[2][1]()
        home = fromFinleyReducedData(Data(2.5, ReducedFunction(dom.toFinley())),
                                     dom)
        self.assertLess(Lsup(home - 2.5), 1e-14)


class Test_BoundaryTransfer3D(unittest.TestCase):
    """
    FunctionOnBoundary in 3D, where a boundary quad becomes two triangles.

    Not a permutation, unlike 2D. Outbound the quad's four points determine a
    function BILINEAR IN THE FACE'S OWN TWO AXES, so anything of that form
    crosses exactly; inbound each triangle's three points determine an affine
    function, so a linear field comes home exactly.

    The reduced space carries one value per face and can do no better than
    replicate it onto the two triangles, so unlike 2D it does NOT reproduce a
    varying field outbound. What it does keep is the surface integral and the
    round trip, and those are what is asserted.

    A field that VARIES along the boundary is what makes a wrong face pairing
    visible - the mesh view lists its boundary faces in one order and the
    domain its FunctionOnBoundary samples in another unless the two are built
    to agree, and a constant field would survive the mismatch. The normal is
    the sharpest form of that check.
    """
    def test_face_bilinear_is_exact_outbound(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x, xf = FunctionOnBoundary(dom).getX(), FunctionOnBoundary(fin).getX()
            u = x[0]*x[1] + x[1]*x[2] + x[2]*x[0]
            uf = xf[0]*xf[1] + xf[1]*xf[2] + xf[2]*xf[0]
            err = Lsup(toFinleyBoundaryData(u, fin) - uf)
            self.assertLess(err, 1e-12, "%s: %g" % (name, err))

    def test_linear_is_exact_inbound(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            for label, fs in (("FunctionOnBoundary", FunctionOnBoundary),
                              ("ReducedFunctionOnBoundary",
                               ReducedFunctionOnBoundary)):
                x, xf = fs(dom).getX(), fs(fin).getX()
                err = Lsup(fromFinleyBoundaryData(1.+2.*xf[0]+3.*xf[1]-xf[2],
                                                  dom)
                           - (1.+2.*x[0]+3.*x[1]-x[2]))
                self.assertLess(err, 1e-12, "%s/%s: %g" % (name, label, err))

    def test_round_trip_is_exact_for_a_linear_field(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            for label, fs in (("FunctionOnBoundary", FunctionOnBoundary),
                              ("ReducedFunctionOnBoundary",
                               ReducedFunctionOnBoundary)):
                x = fs(dom).getX()
                u = 1. + 2.*x[0] + 3.*x[1] - x[2]
                err = Lsup(fromFinleyBoundaryData(
                        toFinleyBoundaryData(u, fin), dom) - u)
                self.assertLess(err, 1e-12, "%s/%s: %g" % (name, label, err))

    def test_reduced_round_trip_is_the_identity(self):
        """
        replicate onto the two triangles, then average them by area: the
        weights sum to one, so ANY field survives - nonlinear included
        """
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = ReducedFunctionOnBoundary(dom).getX()
            u = sin(3.*x[0]) * cos(2.*x[1]) * x[2]
            err = Lsup(fromFinleyBoundaryData(
                    toFinleyBoundaryData(u, fin), dom) - u)
            self.assertLess(err, 1e-14, "%s: %g" % (name, err))

    def test_surface_integral_is_preserved(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            for label, fs in (("FunctionOnBoundary", FunctionOnBoundary),
                              ("ReducedFunctionOnBoundary",
                               ReducedFunctionOnBoundary)):
                x = fs(dom).getX()
                u = x[0]*x[1] + x[1]*x[2] + x[2]*x[0]
                a = integrate(u)
                b = integrate(toFinleyBoundaryData(u, fin))
                self.assertAlmostEqual(a, b, 10, "%s/%s: surface integral "
                                       "changed (%.12g -> %.12g)"
                                       % (name, label, a, b))

    def test_the_normal_survives(self):
        """
        The sharpest check that a value lands on the face it came from: the
        normal differs between neighbouring faces, and both triangles of a quad
        must receive the same one. It also says the split kept the winding, so
        the outward normal stayed outward.
        """
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            moved = toFinleyBoundaryData(FunctionOnBoundary(dom).getNormal(),
                                         fin)
            err = Lsup(moved - FunctionOnBoundary(fin).getNormal())
            self.assertLess(err, 1e-14, "%s: %g" % (name, err))


class Test_FunctionTransfer3D(unittest.TestCase):
    """
    Function in 3D: 2x2x2 Gauss points on an octant against a four-point rule
    on each of six Tet4s.

      outbound  the octant's eight values are unisolvent for a TRILINEAR
                function, so anything trilinear crosses exactly.
      inbound   a tet's four values are unisolvent for an AFFINE function, so a
                linear field comes home exactly. Which tet to read a point from
                is decided from the tets' vertices, since their quadrature
                points span only part of them.

    The volume integral is preserved outbound for any field. That is measured
    rather than derived: the four-point rule is not exact for the xyz term of
    the interpolant on a single tet, but the error cancels over the six of the
    cone. The oblong case above is there so this is not read off a symmetric
    mesh alone.
    """
    TOL = 1e-12

    def test_trilinear_is_exact_outbound(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x, xf = Function(dom).getX(), Function(fin).getX()
            err = Lsup(toFinleyFunctionData(x[0]*x[1]*x[2], fin)
                       - xf[0]*xf[1]*xf[2])
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_linear_is_exact_inbound(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x, xf = Function(dom).getX(), Function(fin).getX()
            err = Lsup(fromFinleyFunctionData(1.+2.*xf[0]+3.*xf[1]-xf[2], dom)
                       - (1.+2.*x[0]+3.*x[1]-x[2]))
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_round_trip_is_exact_for_a_linear_field(self):
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = Function(dom).getX()
            u = 1. + 2.*x[0] + 3.*x[1] - x[2]
            err = Lsup(fromFinleyFunctionData(toFinleyFunctionData(u, fin), dom)
                       - u)
            self.assertLess(err, self.TOL, "%s: %g" % (name, err))

    def test_integral_is_preserved_outbound(self):
        """
        On a CONFORMING octant only. The four-point rule is not exact for the
        xyz term of the trilinear interpolant on a single tetrahedron, and what
        makes the integral come out right is that the error cancels over the six
        of the cone. A hanging octant is coned from its centre into anything
        from twelve to forty-eight tetrahedra of unequal volume, and nothing
        makes those errors cancel - measured at a few times 1e-6 relative on the
        simplest seam. What survives there is the linear case below.
        """
        for name, make in CONFORMING_FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = Function(dom).getX()
            w = sin(2.*x[0]) * cos(x[1]) * exp(x[2])
            a, b = integrate(w), integrate(toFinleyFunctionData(w, fin))
            self.assertAlmostEqual(a, b, 10, "%s: %.12g -> %.12g" % (name, a, b))

    def test_a_linear_integral_is_preserved_on_any_forest(self):
        """
        A linear field is reproduced exactly at every point of the octant, and
        the four-point rule is exact to degree two, so the tetrahedra integrate
        it exactly however the octant was cut. This is the part of the claim
        above that does not depend on a cancellation.
        """
        for name, make in FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = Function(dom).getX()
            w = 1. + 2.*x[0] + 3.*x[1] - x[2]
            a, b = integrate(w), integrate(toFinleyFunctionData(w, fin))
            self.assertAlmostEqual(a, b, 12, "%s: %.12g -> %.12g" % (name, a, b))


class Test_ComplexTransfer3D(unittest.TestCase):
    """
    The same four transfers on complex data in 3D. As in 2D there is no
    separate machinery - a complex field is a real one with twice as many
    components - so what this checks is that the two parts are neither mixed
    nor dropped, which is why they are different functions of position here.
    """
    def setUp(self):
        self.dom = CONFORMING_FORESTS_3D[0][1]()
        self.fin = self.dom.toFinley()

    def tearDown(self):
        del self.dom
        del self.fin

    def field(self, fs, domain):
        x = fs(domain).getX()
        return (1. + 2.*x[0] + 3.*x[1] - x[2]) + 1j * (0.5 - x[0] + 4.*x[2])

    def test_continuous_function(self):
        u = self.field(ContinuousFunction, self.dom)
        uf = toFinleyData(u, self.fin)
        self.assertTrue(uf.isComplex(), "the result lost its complexity")
        self.assertEqual(Lsup(uf - self.field(ContinuousFunction, self.fin)), 0.)
        self.assertEqual(Lsup(fromFinleyData(uf, self.dom) - u), 0.)

    def test_reduced_function(self):
        u = self.field(ReducedFunction, self.dom)
        uf = toFinleyReducedData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertLess(Lsup(fromFinleyReducedData(uf, self.dom) - u), 1e-14)

    def test_function_on_boundary(self):
        u = self.field(FunctionOnBoundary, self.dom)
        uf = toFinleyBoundaryData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertLess(Lsup(uf - self.field(FunctionOnBoundary, self.fin)),
                        1e-12)
        self.assertLess(Lsup(fromFinleyBoundaryData(uf, self.dom) - u), 1e-12)

    def test_function(self):
        u = self.field(Function, self.dom)
        uf = toFinleyFunctionData(u, self.fin)
        self.assertTrue(uf.isComplex())
        self.assertLess(Lsup(uf - self.field(Function, self.fin)), 1e-12)
        self.assertLess(Lsup(fromFinleyFunctionData(uf, self.dom) - u), 1e-12)


class Test_GradedForestExports3D(unittest.TestCase):
    """
    The tetrahedral split on a forest with 2:1 seams.

    What is NOT tested here is conformity, and deliberately: a globally linear
    field lies in both triangulations of a planar quad, so a patch test passes
    at machine precision on a mesh full of cracks, and so do the volume and the
    surface area. Only the combinatorial face hash sees a mismatched diagonal,
    which is why that check lives inside the converter and runs on every export.
    These are the checks that catch the OTHER failures - a piece of the forest
    left uncovered, a face element that missed its triangles, a tetrahedron
    turned inside out.
    """
    def test_the_split_covers_the_forest(self):
        for name, make in GRADED_FORESTS_3D:
            fin = make().toFinley()
            v = integrate(Scalar(1., Function(fin)))
            self.assertAlmostEqual(v, 1., 10,
                                   "%s: the tetrahedra fill %.15g of the unit "
                                   "cube" % (name, v))

    def test_the_boundary_is_closed(self):
        """
        Six unit faces. A boundary quad whose edges carry hanging midpoints is
        a polygon of up to eight vertices and becomes up to six triangles; if
        any of them were dropped, or emitted twice, the area would say so.
        """
        for name, make in GRADED_FORESTS_3D:
            fin = make().toFinley()
            a = integrate(Scalar(1., FunctionOnBoundary(fin)))
            self.assertAlmostEqual(a, 6., 10, "%s: surface area %.15g"
                                   % (name, a))

    def test_a_linear_field_keeps_its_gradient(self):
        for name, make in GRADED_FORESTS_3D:
            fin = make().toFinley()
            x = fin.getX()
            g = grad(x[0] + 2.*x[1] - 3.*x[2], Function(fin))
            err = Lsup(g - [1., 2., -3.])
            self.assertLess(err, 1e-13, "%s: %g" % (name, err))

    def test_the_debug_hex_path_is_still_refused(self):
        """
        One element per octant cannot resolve a 2:1 seam, in any dimension.
        """
        dom = oxley_meshes.graded3D(oxley_meshes.SEAM_3D)
        self.assertFalse(dom.isConforming())
        self.assertRaises(RuntimeError, lambda: dom.toFinley(simplices=False))


class Test_GradedContinuousFunctionTransfer3D(unittest.TestCase):
    """
    The one claim that changes on a graded forest.

    The export has nodes the forest does not: a node at every seam position,
    and one at the centre of every octant that has a hanging node on it - the
    apex its tetrahedra are coned from, no corner of such an octant serving.
    Each of them is defined as the AVERAGE of octant corners, so:

      - a LINEAR field is reproduced at every one of them, and crosses exactly;
      - a field that is not linear differs at those nodes, and there alone;
      - coming home is exact for ANY field, because oxley's nodes are a subset
        of the export's and were copied, not interpolated.

    The last two together are what says the extra nodes are the only ones
    involved: if the copy were wrong anywhere else, the round trip would show it.
    """
    def test_linear_is_exact_outbound(self):
        for name, make in GRADED_FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x, xf = ContinuousFunction(dom).getX(), ContinuousFunction(fin).getX()
            err = Lsup(toFinleyData(1. + 2.*x[0] + 3.*x[1] - x[2], fin)
                       - (1. + 2.*xf[0] + 3.*xf[1] - xf[2]))
            self.assertLess(err, 1e-14, "%s: %g" % (name, err))

    def test_a_nonlinear_field_is_not_exact_outbound(self):
        """
        The counterpart, and the reason the test above is not vacuous: on a
        CONFORMING forest every field crosses exactly, so a linear one proves
        nothing about the averaging rule unless a nonlinear one is seen to fail.
        """
        dom = oxley_meshes.graded3D(oxley_meshes.MIXED_3D)
        fin = dom.toFinley()
        x, xf = ContinuousFunction(dom).getX(), ContinuousFunction(fin).getX()
        u = x[0]*x[0] + x[1]*x[2]
        uf = xf[0]*xf[0] + xf[1]*xf[2]
        self.assertGreater(Lsup(toFinleyData(u, fin) - uf), 1e-6)

    def test_coming_home_is_exact_for_any_field(self):
        for name, make in GRADED_FORESTS_3D:
            dom = make()
            fin = dom.toFinley()
            x = ContinuousFunction(dom).getX()
            u = sin(3.*x[0]) * cos(2.*x[1]) * exp(x[2])
            self.assertEqual(Lsup(fromFinleyData(toFinleyData(u, fin), dom) - u),
                             0., "%s: the way home is not a copy" % name)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
