"""
D2 probe: does a degree-2 p4est_lnodes already number the hanging positions?

Background. The export needs a global id for every node of the simplex mesh,
including the hanging positions (a coarse octant's face centre or edge midpoint).
A degree-1 lnodes does not number those: they are not nodes at all, and the
corresponding element_nodes slot holds a far master instead. The plan of record
(D2) was therefore to number them ourselves - own the ones whose coarse octant we
own, canonical key (global quadrant id, local face/edge index), one MPI_Exscan for
the block offsets.

But at degree 2 those very positions ARE independent nodes: the 3x3(x3) Lobatto
layout puts a node at every face centre and edge midpoint. p4est documents
(p4est_lnodes.h:70-73) that a hanging slot "indexes the geometrically
corresponding independent node of a neighbor" - so the hanging corner of a fine
octant, which sits exactly at the coarse neighbour's face centre or edge midpoint,
should come back carrying that node's global id. If that holds, D2 is free: take
the corner slots of a degree-2 lnodes and p4est has done the numbering, the
ownership and the sharers for us.

The test. Map every corner slot of every octant to its geometric position and
check the id <-> position pairing both ways round. An id at two positions means
the slot held a far master (the degree-1 behaviour, hypothesis dead). A position
under two ids means the position got numbered twice, which would crack the mesh
between the elements that disagree. Both counts must be zero.

Run:
    ./bin/run-escript -n1 oxley-designs/probe_lnodes_degree2.py
    ./bin/run-escript -n2 -t1 oxley-designs/probe_lnodes_degree2.py

The per-rank check above is necessary but NOT sufficient under MPI: it cannot see
two ranks giving the same shared position different ids. Serial is the decisive
first test; the cross-rank check comes after.
"""
import esys.escript as esc
import esys.finley  # noqa: F401  - registers the type toFinley() returns
import esys.oxley as oxley

rank = esc.getMPIRankWorld()
failures = []


def probe(name, dom, expect_hanging):
    """
    Prints the counts and checks them against the ANSWER THIS PROBE ESTABLISHED
    (2026-07-29), so a change in p4est's behaviour shows up as a failure here:

      conforming forest     -> the corner slots are a clean bijection between ids
                               and positions;
      non-conforming forest -> they are not. A hanging corner slot holds a MASTER,
                               so the hanging position gets no id of its own from
                               the fine side, and the master's id turns up at two
                               positions. See probe_lnodes_slots.py for which id
                               lands where, and note that the position IS numbered
                               - on the COARSE element's non-corner slot.
    """
    r = dom.lnodesDegree2Report()
    (octants, slots, ids, clashes, hanging, localNodes,
     positions, splits) = r

    consistent = (clashes == 0 and splits == 0 and ids == positions)
    if expect_hanging:
        if hanging == 0:
            ok, note = False, "  *** forest is CONFORMING, nothing was tested ***"
        elif consistent:
            ok, note = False, "  *** slots are now CONSISTENT - p4est changed? ***"
        else:
            ok, note = True, "  (as established: masters at the hanging corners)"
    else:
        ok = consistent
        note = "" if ok else "  *** a CONFORMING forest must be consistent ***"

    if rank == 0:
        print("  %-34s octants=%-7d hangingOctants=%-6d nodes=%-7d"
              % (name, octants, hanging, localNodes))
        print("  %-34s slots=%-9d ids=%-8d positions=%-8d"
              % ("", slots, ids, positions))
        print("  %-34s idAtSeveralPositions=%-5d positionWithSeveralIds=%-5d  %s%s"
              % ("", clashes, splits, "as expected" if ok else "UNEXPECTED", note))
    if not ok:
        failures.append(name)
    return ok


if rank == 0:
    print("degree-2 lnodes probe: are the hanging positions numbered?  "
          "(%d rank(s))" % esc.getMPISizeWorld())
    print("2D")

# Non-conforming forests come from per-block refine_level: blocks at different
# levels leave a 2:1 seam along their shared boundary. (refineRegion/refinePoint
# are the older adaptive path and are not used here.)

# conforming control: no hanging anywhere, so this only checks that the corner
# slots of an ordinary forest behave as expected
probe("uniform L=1 (control)", oxley.Rectangle(n0=2, n1=2, refine_level=1),
      expect_hanging=False)

# one seam
probe("per-block [[2],[1]]",
      oxley.Rectangle(n0=2, n1=1, refine_level=[[2], [1]]), expect_hanging=True)

# a deep block, so the balance cascade puts seams of several depths next to each
# other and hanging nodes reach the outer boundary
probe("per-block [[3,1],[1,1]]",
      oxley.Rectangle(n0=2, n1=2, refine_level=[[3, 1], [1, 1]]),
      expect_hanging=True)

probe("per-block 3x3 mixed",
      oxley.Rectangle(n0=3, n1=3, refine_level=[[3, 1, 2], [1, 2, 1], [2, 1, 3]]),
      expect_hanging=True)

if rank == 0:
    print("3D")

probe("uniform L=1 (control)",
      oxley.Brick(n0=2, n1=2, n2=2, refine_level=1), expect_hanging=False)

# one seam, and the 3D-only cases with it: octants meeting the coarse neighbour
# along an EDGE or at a CORNER only, not across a face
probe("per-block [[[2]],[[1]]]",
      oxley.Brick(n0=2, n1=1, n2=1, refine_level=[[[2]], [[1]]]),
      expect_hanging=True)

probe("per-block 2x2x2 one deep block",
      oxley.Brick(n0=2, n1=2, n2=2,
                  refine_level=[[[3, 1], [1, 1]], [[1, 1], [1, 1]]]),
      expect_hanging=True)

if rank == 0:
    if failures:
        print("RESULT: UNEXPECTED  (%s)" % ", ".join(failures))
    else:
        print("RESULT: as established - the degree-2 corner slots are consistent "
              "only where nothing hangs, so the fine side cannot read the id of a "
              "hanging position and oxley must number those itself (D2)")
