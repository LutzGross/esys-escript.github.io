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
from esys.oxley import Rectangle, Brick, RefinementQueue2D, RefinementQueue3D

N0 = 2
N1 = 2
N2 = 2
L0 = 1.
L1 = 1.
L2 = 1.


def numElements(domain):
    """elements, not quadrature points: 4 per element in 2D, 8 in 3D"""
    perElement = 4 if domain.getDim() == 2 else 8
    return Data(1., Function(domain)).getNumberOfDataPoints() // perElement


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
        self.assertGreater(sr.getNumberOfDataPoints(),
                           s.getNumberOfDataPoints(),
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

    def test_refine_border_by_name(self):
        """the Border enum is not exposed to python, so names are the API"""
        for name in ("top", "bottom", "left", "right", "north", "south",
                     "east", "west"):
            f = RefinementQueue2D()
            f.setRefinementLevel(1)
            f.refineBorder(border=name, dx=0.3)
            self.assertGreater(numElements(f.apply(self.domain)),
                               numElements(self.domain),
                               "border '%s' refined nothing" % name)

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

    def test_wrong_dimension_is_refused(self):
        rect = Rectangle(n0=N0, n1=N1, l0=L0, l1=L1)
        f = RefinementQueue3D()
        f.setRefinementLevel(1)
        f.refineUniform()
        self.assertRaises(RuntimeError, f.apply, rect)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
