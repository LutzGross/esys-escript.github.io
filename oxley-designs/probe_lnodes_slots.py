"""
Follow-up to probe_lnodes_degree2.py: WHAT does a hanging slot hold at degree 2?

The report says the degree-2 corner slots are not a consistent numbering of the
corner positions, so the hypothesis "a hanging corner slot carries the global id
of the coarse neighbour's face-centre / edge-midpoint node, which sits at exactly
that position" is wrong. This script says what the slot holds instead, because
the alternatives lead to different work:

  (a) the slot holds a MASTER of the hanging node (the degree-1 behaviour). Then
      the fine side cannot see the id at all and D2 stands as planned: number the
      hanging positions ourselves.
  (b) the slot holds the right node but our slot->position map is wrong. Then the
      probe is at fault, not p4est.

Method: build the smallest non-conforming forest (two blocks, levels 2 and 1, one
seam), then for each element with face_code != 0 print its corner slots, the id
each holds, and where that id's node actually is - "actually" taken from the slots
of CONFORMING elements only, which are the trustworthy ones.

Run: ./bin/run-escript -n1 oxley-designs/probe_lnodes_slots.py
"""
import esys.escript as esc
import esys.oxley as oxley

TOL = 1e-9


def canonical_positions(slots):
    """id -> position, taken only from elements where nothing hangs"""
    pos = {}
    for (e, slot, gid, x, y, z, fc, corner) in slots:
        if fc == 0:
            pos.setdefault(gid, (x, y, z))
    return pos


def dist(a, b):
    return max(abs(a[i] - b[i]) for i in range(3))


def report(name, dom, dim):
    slots = dom.lnodesDegree2Slots()
    pos = canonical_positions(slots)
    per_elem = {}
    for rec in slots:
        per_elem.setdefault(rec[0], []).append(rec)

    print("\n%s (%dD): %d elements" % (name, dim, len(per_elem)))

    shown = 0
    for e in sorted(per_elem):
        recs = per_elem[e]
        fc = recs[0][6]
        if fc == 0:
            continue
        bad = [r for r in recs
               if r[7] and r[2] in pos and dist(pos[r[2]], (r[3], r[4], r[5])) > TOL]
        if not bad:
            continue
        shown += 1
        if shown > 3:
            break
        print("  element %d, face_code %d" % (e, fc))
        for (_, slot, gid, x, y, z, _, corner) in recs:
            if not corner:
                continue
            here = (x, y, z)
            there = pos.get(gid)
            if there is None:
                where = "id not seen on any conforming element"
            elif dist(there, here) <= TOL:
                where = "id lives here"
            else:
                where = "id lives at (%g,%g,%g)  <-- MISMATCH" % there
            print("    slot %2d at (%g,%g,%g)  id %-6d  %s"
                  % (slot, x, y, z, gid, where))

        # Is the mismatched id one of this element's OTHER corners - i.e. a
        # master of the hanging node, the degree-1 behaviour?
        corner_ids = set(r[2] for r in recs if r[7])
        for r in bad:
            others = [o for o in recs
                      if o[7] and o[2] == r[2] and o[1] != r[1]]
            print("    -> slot %d's id also sits at this element's slot(s) %s"
                  % (r[1], [o[1] for o in others] if others else "none"))

    # Which corner positions never get an id of their own from a CORNER slot?
    # These are the hanging positions: the fine side only ever sees masters there.
    def rnd(x, y, z):
        return (round(x, 9), round(y, 9), round(z, 9))

    ids_at_corner = {}
    for (e, slot, gid, x, y, z, fc, corner) in slots:
        if corner:
            ids_at_corner.setdefault(rnd(x, y, z), set()).add(gid)
    unnumbered = [p for p, g in ids_at_corner.items()
                  if not any(p == rnd(*pos[gid]) for gid in g if gid in pos)]
    print("  corner positions with no id of their own: %d of %d"
          % (len(unnumbered), len(ids_at_corner)))

    # ...but the position IS the centre of a coarse face (or the midpoint of a
    # coarse edge), so it is a degree-2 node OF THE COARSE ELEMENT, which holds it
    # in a non-corner slot. If that id exists, the degree-2 numbering does cover
    # the hanging positions - it is only unreadable from the fine side.
    all_at = {}
    for rec in slots:
        all_at.setdefault(rnd(rec[3], rec[4], rec[5]), []).append(rec)

    covered = 0
    for p in unnumbered:
        owners = [r for r in all_at[p]
                  if not r[7] and r[2] in pos and rnd(*pos[r[2]]) == p]
        if owners:
            covered += 1
    print("  ...of which ARE numbered, on a non-corner slot elsewhere: %d of %d"
          % (covered, len(unnumbered)))
    for p in sorted(unnumbered)[:4]:
        owners = [r for r in all_at[p]
                  if not r[7] and r[2] in pos and rnd(*pos[r[2]]) == p]
        print("     %s  corner slots hold %s;  own id %s"
              % (p, sorted(ids_at_corner[p]),
                 ("%d (element %d slot %d)" % (owners[0][2], owners[0][0], owners[0][1]))
                 if owners else "NONE"))


report("two blocks, levels 2 and 1",
       oxley.Rectangle(n0=2, n1=1, refine_level=[[2], [1]]), 2)
report("two blocks, levels 2 and 1",
       oxley.Brick(n0=2, n1=1, n2=1, refine_level=[[[2]], [[1]]]), 3)
