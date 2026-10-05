########################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# Earth Systems Science Computational Center (ESSCC)
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
Earth Systems Science Computational Center (ESSCC)
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
RefinementQueue2D / RefinementQueue3D.

Refinement used to be a set of methods ON the domain that changed it in place,
so every Data already defined over that domain was silently stale afterwards
and nothing in the type system said so. A queue instead collects the
operations and hands back a NEW domain, leaving the source alone - which is
the property most of the tests below are about.

The refinement a domain is BORN with is a different thing and stays a
constructor argument: Rectangle/Brick take refine_level, an int or a per-block
list. A level has to mean the same in both, which is what
test_uniform_matches_the_constructor pins down.
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
# finley must be imported for toFinley() to hand back a usable domain
import esys.finley
from esys.oxley import Rectangle, Brick, RefinementQueue2D, RefinementQueue3D, \
                       fromFinleyData

N0 = 2
N1 = 2
N2 = 2
L0 = 1.
L1 = 1.
L2 = 1.


def refinedFraction(domain, coarse):
    """the fraction of the volume covered by elements finer than coarse's"""
    h = interpolate(domain.getSize(), Function(domain))
    vol = integrate(Scalar(1., Function(domain)))
    return integrate(whereNegative(h - 0.99*Lsup(coarse.getSize()))) / vol


def numElements(domain):
    """
    the global number of elements of domain. getMeshInfo() gives the elements
    this rank owns - each element belongs to exactly one rank, whereas the
    samples of a Function also count the overlap shared with neighbouring
    ranks - and they are summed over the communicator the domain lives on,
    which need not be all ranks. Without mpi4py there is no communicator to
    hand back, but then the domain lives on all ranks.
    """
    local = int(domain.getMeshInfo()["numElements"])
    comm = domain.getMPIComm()
    if comm is None:
        return getMPIWorldSum(local)
    return comm.allreduce(local)

class Test_RefinementQueue2D(unittest.TestCase):
    def setUp(self):
        self.domain = Rectangle(n0=N0, n1=N1, l0=L0, l1=L1)

    def tearDown(self):
        del self.domain

    def test_apply_returns_a_new_domain(self):
        before = numElements(self.domain)
        f = RefinementQueue2D()
        f.setRefinementLevel(1)
        f.refineUniform()
        refined = f.apply(self.domain)
        self.assertGreater(numElements(refined), before,
                           "apply() did not refine anything")
        self.assertEqual(numElements(self.domain), before,
                         "apply() modified the domain it was given; it must "
                         "work on a copy")

    def test_data_over_the_source_stays_valid(self):
        """
        The reason the domain has no refine methods at all.

        With in-place refinement this sequence was silently broken:

            dom = Rectangle(...)
            s  = Scalar(1., ContinuousFunction(dom))
            dom.refine()                       # dom now has more nodes
            sr = Scalar(1., ContinuousFunction(dom))
            s + sr                             # s is sized for the OLD dom

        s and sr claim the same function space on the same domain object while
        being different lengths, and nothing in the types says otherwise.
        Refining into a NEW domain makes them visibly different objects: s
        stays valid over the coarse mesh, sr belongs to the refined one, and
        mixing them is a plain error instead of a silent one.
        """
        s = Scalar(1., ContinuousFunction(self.domain))
        f = RefinementQueue2D()
        f.setRefinementLevel(1)
        f.refineUniform()
        refined = f.apply(self.domain)
        sr = Scalar(1., ContinuousFunction(refined))

        # the old Data is still usable and still sized for the coarse mesh
        self.assertEqual(Lsup(s + s), 2.0, "Data over the source domain broke")
        self.assertGreater(getMPIWorldSum(sr.getNumberOfDataPoints()),
                           getMPIWorldSum(s.getNumberOfDataPoints()),
                           "the refined domain should carry more nodes")
        # and the two cannot be mixed by accident
        self.assertRaises(RuntimeError, lambda: Lsup(s + sr))

    def test_uniform_matches_the_constructor(self):
        """a level must mean the same however it is asked for"""
        for level in (1, 2):
            f = RefinementQueue2D()
            f.setRefinementLevel(level)
            f.refineUniform()
            viaQueue = numElements(f.apply(self.domain))
            viaCtor = numElements(Rectangle(n0=N0, n1=N1, l0=L0, l1=L1,
                                            refine_level=level))
            self.assertEqual(viaQueue, viaCtor,
                             "level %d: queue gives %d elements, the "
                             "constructor %d" % (level, viaQueue, viaCtor))

    def test_level_is_absolute(self):
        """
        A level says what the elements should BE, not how many times to halve
        them: refine_uniform subdivides while quadrant->level < level. So
        re-applying the same queue changes nothing, and only asking for a
        higher level refines further. Worth pinning down, because "refine to
        level 1" and "refine once more" read alike and are not the same.
        """
        one = RefinementQueue2D()
        one.setRefinementLevel(1)
        one.refineUniform()
        once = one.apply(self.domain)
        self.assertEqual(numElements(one.apply(once)), numElements(once),
                         "a level is absolute, so re-applying it is a no-op")

        two = RefinementQueue2D()
        two.setRefinementLevel(2)
        two.refineUniform()
        self.assertGreater(numElements(two.apply(once)), numElements(once),
                           "a higher level must refine further")

    def test_refine_region_is_local(self):
        f = RefinementQueue2D()
        f.setRefinementLevel(2)
        # the region test is on the quadrant's ORIGIN corner, so the box has
        # to contain one: at level 0 the corners are (0,0), (.5,0), (0,.5),
        # (.5,.5), and only (0,0) is inside this one
        f.refineRegion(x0=0.0, y0=0.0, x1=0.4, y1=0.4)
        refined = f.apply(self.domain)
        self.assertGreater(numElements(refined), numElements(self.domain))
        # a region is not the whole domain, so it must cost less than refining
        # everything to the same level
        g = RefinementQueue2D()
        g.setRefinementLevel(2)
        g.refineUniform()
        self.assertLess(numElements(refined), numElements(g.apply(self.domain)))

    def test_refine_point(self):
        f = RefinementQueue2D()
        f.setRefinementLevel(2)
        f.refinePoint(x0=0.55, y0=0.55)
        self.assertGreater(numElements(f.apply(self.domain)),
                           numElements(self.domain))

    def test_refine_circle(self):
        """
        A circle refines every element it overlaps. It used to refine only an
        element whose centre lay exactly on the circle, i.e. none at all.
        """
        f = RefinementQueue2D()
        f.setRefinementLevel(2)
        f.refineCircle(x0=0.3, y0=0.6, r=0.1)
        refined = f.apply(self.domain)
        self.assertGreater(numElements(refined), numElements(self.domain),
                           "refineCircle refined nothing")
        g = RefinementQueue2D()
        g.setRefinementLevel(2)
        g.refineUniform()
        self.assertLess(numElements(refined), numElements(g.apply(self.domain)),
                        "a circle is not the whole domain")

    def test_refine_small_circle(self):
        """a circle inside one element, touching none of its corners"""
        f = RefinementQueue2D()
        f.setRefinementLevel(1)
        f.refineCircle(x0=0.25, y0=0.25, r=0.01)
        self.assertGreater(numElements(f.apply(self.domain)),
                           numElements(self.domain))

    def test_refine_border_by_name(self):
        """
        The Border enum is not exposed to python, so names are the API. On a
        2 x 1 domain of 8 x 4 blocks, a border strip thinner than a block
        refines exactly the layer of blocks along that border: a quarter of
        the domain along top or bottom, an eighth along left or right. This
        used to refine nothing, everything or an arbitrary part of the
        domain as soon as it was not the unit square.
        """
        coarse = Rectangle(n0=8, n1=4, l0=2., l1=1.)
        for name, fraction in (("top", 0.25), ("north", 0.25),
                               ("bottom", 0.25), ("south", 0.25),
                               ("left", 0.125), ("west", 0.125),
                               ("right", 0.125), ("east", 0.125)):
            f = RefinementQueue2D()
            f.refineBorder(border=name, dx=0.1, level=2)
            self.assertAlmostEqual(refinedFraction(f.apply(coarse), coarse),
                                   fraction, places=10,
                                   msg="border '%s' refined the wrong part" % name)

    def disc(self, domain, r=0.1):
        x = Function(domain).getX()
        return whereNegative(length(x - [0.3, 0.6]) - r)

    def test_refine_mask(self):
        """
        The queue holds only the name of a mask, a template not tied to any
        domain; the mask itself comes with apply()
        """
        f = RefinementQueue2D()
        f.refineMask("fault", level=3)
        refined = f.apply(self.domain, fault=self.disc(self.domain))
        self.assertGreater(numElements(refined), numElements(self.domain),
                           "refineMask refined nothing")
        g = RefinementQueue2D()
        g.refineUniform(level=3)
        self.assertLess(numElements(refined), numElements(g.apply(self.domain)),
                        "a mask is not the whole domain")
        # the same queue on another domain, with a mask on that domain
        other = Rectangle(n0=4, n1=4, l0=L0, l1=L1)
        self.assertGreater(numElements(f.apply(other, fault=self.disc(other))),
                           numElements(other))

    def test_zero_mask_refines_nothing(self):
        """a constant Data stores one value, not one per quadrature point"""
        f = RefinementQueue2D()
        f.refineMask("fault", level=3)
        self.assertEqual(numElements(f.apply(self.domain,
                                     fault=Scalar(0., Function(self.domain)))),
                         numElements(self.domain))

    def test_mask_anywhere_in_the_queue(self):
        """
        A mask is read on the domain passed to apply(), so it no longer has
        to come first: refining around it after a uniform refinement must
        give more than the uniform refinement alone
        """
        f = RefinementQueue2D()
        f.refineUniform(level=1)
        f.refineMask("fault", level=3)
        g = RefinementQueue2D()
        g.refineUniform(level=1)
        self.assertGreater(numElements(f.apply(self.domain,
                                       fault=self.disc(self.domain))),
                           numElements(g.apply(self.domain)))

    def test_mask_must_match_the_tags(self):
        f = RefinementQueue2D()
        f.refineMask("fault")
        mask = self.disc(self.domain)
        # no mask for a tag
        self.assertRaises(RuntimeError, f.apply, self.domain)
        # a mask for no tag: most likely a misspelt one
        self.assertRaises(RuntimeError, lambda: f.apply(self.domain,
                                                        fault=mask, fualt=mask))
        # a mask that is not Data, or not a scalar
        self.assertRaises(RuntimeError, lambda: f.apply(self.domain, fault=1.))
        self.assertRaises(RuntimeError, lambda: f.apply(self.domain,
                                        fault=Function(self.domain).getX()))
        # a mask on another domain
        other = Rectangle(n0=N0, n1=N1, l0=L0, l1=L1)
        self.assertRaises(RuntimeError, lambda: f.apply(self.domain,
                                        fault=Scalar(1., Function(other))))

    def test_mask_tag_must_be_a_name(self):
        """the tag is a keyword argument of apply()"""
        f = RefinementQueue2D()
        self.assertRaises(RuntimeError, f.refineMask, "my fault")
        self.assertRaises(RuntimeError, f.refineMask, "")
        self.assertRaises(RuntimeError, f.refineMask, "domain")

    def test_refined_domain_outlives_its_source(self):
        f = RefinementQueue2D()
        f.refineUniform(level=1)
        source = Rectangle(n0=N0, n1=N1, l0=L0, l1=L1)
        refined = f.apply(source)
        del source
        self.assertAlmostEqual(integrate(Scalar(1., Function(refined))),
                               L0*L1, places=10)

    def test_no_pde_on_a_refined_domain(self):
        """
        oxley's assembler ignores hanging nodes, so on a refined forest it
        used to return a wrong solution without complaint
        """
        f = RefinementQueue2D()
        f.refinePoint(x0=0.3, y0=0.6, level=2)
        self.assertRaises(NotImplementedError, LinearPDE, f.apply(self.domain))
        # a uniform forest has no hanging nodes and still assembles
        g = RefinementQueue2D()
        g.refineUniform(level=1)
        LinearPDE(g.apply(self.domain))

    def test_pde_on_the_finley_export(self):
        """
        the way to solve on a refined domain: the linear patch test is exact
        on the finley export, hanging nodes included
        """
        f = RefinementQueue2D()
        f.refinePoint(x0=0.3, y0=0.6, level=3)
        refined = f.apply(self.domain)
        fin = refined.toFinley()
        x = fin.getX()
        exact = 1. + 2*x[0] - x[1]
        pde = LinearPDE(fin)
        pde.setSymmetryOn()
        pde.setValue(A=kronecker(fin), r=exact,
                     q=whereZero(x[0]) + whereZero(x[0]-L0)
                     + whereZero(x[1]) + whereZero(x[1]-L1))
        pde.getSolverOptions().setTolerance(1e-12)
        u = fromFinleyData(pde.getSolution(), refined)
        xo = refined.getX()
        self.assertLess(Lsup(u - (1. + 2*xo[0] - xo[1])), 1e-8)

    def test_unknown_border_is_refused(self):
        f = RefinementQueue2D()
        self.assertRaises(RuntimeError, f.refineBorder, "sideways", 0.3)

    def test_empty_queue_changes_nothing(self):
        f = RefinementQueue2D()
        self.assertEqual(numElements(f.apply(self.domain)),
                         numElements(self.domain))

    def test_wrong_dimension_is_refused(self):
        brick = Brick(n0=N0, n1=N1, n2=N2, l0=L0, l1=L1, l2=L2)
        f = RefinementQueue2D()
        f.setRefinementLevel(1)
        f.refineUniform()
        self.assertRaises(RuntimeError, f.apply, brick)


class Test_RefinementQueue3D(unittest.TestCase):
    def setUp(self):
        self.domain = Brick(n0=N0, n1=N1, n2=N2, l0=L0, l1=L1, l2=L2)

    def tearDown(self):
        del self.domain

    def test_apply_returns_a_new_domain(self):
        before = numElements(self.domain)
        f = RefinementQueue3D()
        f.setRefinementLevel(1)
        f.refineUniform()
        refined = f.apply(self.domain)
        self.assertGreater(numElements(refined), before)
        self.assertEqual(numElements(self.domain), before,
                         "apply() modified the domain it was given")

    def test_uniform_matches_the_constructor(self):
        f = RefinementQueue3D()
        f.setRefinementLevel(1)
        f.refineUniform()
        self.assertEqual(numElements(f.apply(self.domain)),
                         numElements(Brick(n0=N0, n1=N1, n2=N2, l0=L0, l1=L1,
                                           l2=L2, refine_level=1)))

    def test_refine_region(self):
        f = RefinementQueue3D()
        f.setRefinementLevel(2)
        # see the 2D case: the box must contain a quadrant origin
        f.refineRegion(x0=0.0, y0=0.0, z0=0.0, x1=0.4, y1=0.4, z1=0.4)
        self.assertGreater(numElements(f.apply(self.domain)),
                           numElements(self.domain))

    def test_empty_queue_keeps_the_refinement(self):
        """apply() used to rebuild a Brick at level 0, dropping its refinement"""
        brick = Brick(n0=N0, n1=N1, n2=N2, l0=L0, l1=L1, l2=L2, refine_level=2)
        self.assertEqual(numElements(RefinementQueue3D().apply(brick)),
                         numElements(brick))

    def test_no_pde_on_a_refined_domain(self):
        f = RefinementQueue3D()
        f.refinePoint(x0=0.3, y0=0.6, z0=0.4, level=2)
        self.assertRaises(NotImplementedError, LinearPDE, f.apply(self.domain))

    def test_refine_mask(self):
        x = Function(self.domain).getX()
        f = RefinementQueue3D()
        f.refineMask("blob", level=2)
        refined = f.apply(self.domain,
                          blob=whereNegative(length(x - [0.3, 0.6, 0.4]) - 0.1))
        self.assertGreater(numElements(refined), numElements(self.domain))
        g = RefinementQueue3D()
        g.refineUniform(level=2)
        self.assertLess(numElements(refined), numElements(g.apply(self.domain)))

    def test_refine_border_by_name(self):
        """
        In 3D top and bottom are the faces normal to z, north (back) and
        south (front) those normal to y, and left (west) and right (east)
        those normal to x. See the 2D case for the fractions.
        """
        coarse = Brick(n0=8, n1=4, n2=4, l0=2., l1=1., l2=1.)
        for name, fraction in (("top", 0.25), ("bottom", 0.25),
                               ("north", 0.25), ("back", 0.25),
                               ("south", 0.25), ("front", 0.25),
                               ("left", 0.125), ("west", 0.125),
                               ("right", 0.125), ("east", 0.125)):
            f = RefinementQueue3D()
            f.refineBorder(border=name, dx=0.1, level=1)
            self.assertAlmostEqual(refinedFraction(f.apply(coarse), coarse),
                                   fraction, places=10,
                                   msg="border '%s' refined the wrong part" % name)
        self.assertRaises(RuntimeError, RefinementQueue3D().refineBorder,
                          "sideways", 0.1)

    def test_refine_sphere(self):
        """see the 2D circle: a sphere refines every element it overlaps"""
        f = RefinementQueue3D()
        f.setRefinementLevel(2)
        f.refineSphere(x0=0.3, y0=0.6, z0=0.4, r=0.1)
        refined = f.apply(self.domain)
        self.assertGreater(numElements(refined), numElements(self.domain),
                           "refineSphere refined nothing")
        g = RefinementQueue3D()
        g.setRefinementLevel(2)
        g.refineUniform()
        self.assertLess(numElements(refined), numElements(g.apply(self.domain)),
                        "a sphere is not the whole domain")

    def test_wrong_dimension_is_refused(self):
        rect = Rectangle(n0=N0, n1=N1, l0=L0, l1=L1)
        f = RefinementQueue3D()
        f.setRefinementLevel(1)
        f.refineUniform()
        self.assertRaises(RuntimeError, f.apply, rect)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
