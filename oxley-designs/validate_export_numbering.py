"""
Step 1 of the oxley -> finley export: the node numbering.

Checks the NODE half only - global ids and coordinates. Nothing here builds a
triangle; the simplex split is step 2. What must hold:

  1. AGREEMENT. One physical position carries one export id on every rank that
     sees it, and one export id belongs to one position. This is the property
     the whole design rests on: finley resolves the mesh by global id, so a
     hanging node the coarse rank and the fine rank number differently is two
     nodes, and the mesh is cracked along the seam.
  2. BLOCKS. Every id falls inside exactly one rank's block, and a rank's own
     lnodes nodes fall inside its own. finley's resolveNodeIds() allocates two
     dense arrays spanning the id range of the local elements, so ids that are
     not contiguous per rank cost O(global) memory per rank.
  3. COARSE-SIDE COVERAGE. Wherever an element face carries a hanging node, the
     rank holding that element knows it - via elementFaceHangingNode. This is
     the bug step 1 exists to fix: materialisation used to be driven by
     face_code, which p4est sets on the FINE side only, so a rank holding just
     the coarse octant of a seam never created the node at all.
  4. FINE-SIDE REDIRECT. The fine element's corner slot at a seam points at the
     hanging node rather than at a master lying outside the element.
  5. DECODE. Ids decode back to (owner, kind, oxley node), which is what the
     finley -> oxley transfer will use instead of a stored permutation.

Run at several rank counts, then compare:
    for n in 1 2 3 4; do ./bin/run-escript -n$n -t1 \
        oxley-designs/validate_export_numbering.py --dump $n; done
    python3 oxley-designs/validate_export_numbering.py --report
"""
import glob
import os
import sys

import numpy as np

OUT = os.environ.get("OXLEY_STEP1_DIR", "/tmp/oxley_step1")

# z-order corners of each face, p4est face order: -x, +x, -y, +y
FACE_CORNERS = [(0, 2), (1, 3), (0, 1), (2, 3)]

CASES = [
    ("one seam",            [[2], [3]]),
    ("cascade",             [[3, 1], [1, 2]]),
    ("3x3 mixed",           [[3, 1, 2], [1, 2, 1], [2, 1, 3]]),
    ("isolated coarse",     [[1, 1, 1], [1, 0, 1], [1, 1, 1]]),
    ("big jump",            [[0, 3], [3, 0]]),
    ("diagonal ramp",       [[0, 1, 2, 3], [1, 2, 3, 2], [2, 3, 2, 1], [3, 2, 1, 0]]),
    ("centre peak",         [[0, 0, 1, 0, 0], [0, 1, 2, 1, 0], [1, 2, 3, 2, 1],
                             [0, 1, 2, 1, 0], [0, 0, 1, 0, 0]]),
    ("checkerboard",        [[1, 2, 1, 2], [2, 1, 2, 1], [1, 2, 1, 2], [2, 1, 2, 1]]),
    ("uniform (control)",   2),
]


def key(p):
    """position -> exact hashable key. Test-side only; the numbering itself is
    combinatorial and never compares coordinates."""
    return (int(round(p[0] * 1e9)), int(round(p[1] * 1e9)))


def dump(nranks):
    import esys.escript as esc
    import esys.oxley as oxley

    rank = esc.getMPIRankWorld()
    os.makedirs(OUT, exist_ok=True)
    for ci, (name, levels) in enumerate(CASES):
        n0 = n1 = 2 if isinstance(levels, int) else 0
        if not isinstance(levels, int):
            n0, n1 = len(levels), len(levels[0])
        dom = oxley.Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                              refine_level=levels)
        info = dom.getMeshInfo(True)
        np.savez(os.path.join(OUT, "p%d_c%d_r%d.npz" % (nranks, ci, rank)),
                 coords=info["nodeCoords"].reshape(-1, 2),
                 exportId=info["nodeFinleyId"],
                 gid=info["nodeLnodesId"],
                 en=info["elementNodes"].reshape(-1, info["nodesPerElement"]),
                 efh=info["elementFaceHangingNode"],
                 dist=info["finleyDistribution"],
                 numOwned=info["numOwnedNodes"])
    if rank == 0:
        print("dumped %d cases at %d rank(s)" % (len(CASES), nranks))


class Checker:
    def __init__(self):
        self.failures = []

    def ok(self, cond, what):
        if not cond:
            self.failures.append(what)
        return cond


def check_case(c, label, ranks):
    dist = ranks[0]["dist"]
    owned = [int(d["numOwned"]) for d in ranks]
    nr = len(ranks)
    realOffset = np.concatenate([[0], np.cumsum(owned)])

    # 1. agreement, both directions
    idAt, posOf = {}, {}
    for d in ranks:
        for i, p in enumerate(d["coords"]):
            k, e = key(p), int(d["exportId"][i])
            idAt.setdefault(k, set()).add(e)
            posOf.setdefault(e, set()).add(k)
    split = {k: v for k, v in idAt.items() if len(v) > 1}
    merged = {e: v for e, v in posOf.items() if len(v) > 1}
    c.ok(not split, "%s: %d positions with >1 id" % (label, len(split)))
    c.ok(not merged, "%s: %d ids at >1 position" % (label, len(merged)))

    # 2. blocks
    for r, d in enumerate(ranks):
        e = d["exportId"]
        c.ok(bool(np.all(e >= 0)), "%s: rank %d has unset ids" % (label, r))
        c.ok(bool(np.all(e < dist[nr])),
             "%s: rank %d id past the last block" % (label, r))
        own = e[:owned[r]]
        c.ok(bool(np.all((own >= dist[r]) & (own < dist[r + 1]))),
             "%s: rank %d owned node outside its block" % (label, r))

    # 3. coarse-side coverage, and 4. fine-side redirect
    allPos = set(idAt.keys())
    missing = bad = 0
    for d in ranks:
        efh = d["efh"]
        for ei, en in enumerate(d["en"]):
            for f, (a, b) in enumerate(FACE_CORNERS):
                mid = 0.5 * (d["coords"][en[a]] + d["coords"][en[b]])
                k = key(mid)
                if k not in allPos:
                    continue                    # no hanging node on this face
                h = int(efh[ei * 4 + f]) if len(efh) else -1
                if h < 0:
                    missing += 1
                elif key(d["coords"][h]) != k:
                    bad += 1
    c.ok(missing == 0, "%s: %d faces whose hanging node the element cannot see"
         % (label, missing))
    c.ok(bad == 0, "%s: %d faces pointing at the wrong node" % (label, bad))

    # every element corner must sit at that corner, i.e. no slot still holding a
    # master from outside the element (the pre-materialisation state)
    stray = 0
    for d in ranks:
        for en in d["en"]:
            xy = d["coords"][en]
            if (xy[:, 0].max() - xy[:, 0].min() <= 0 or
                    xy[:, 1].max() - xy[:, 1].min() <= 0):
                stray += 1
    c.ok(stray == 0, "%s: %d degenerate elements" % (label, stray))

    # 5. decode
    wrong = 0
    for r, d in enumerate(ranks):
        for i in range(owned[r]):
            n = int(d["exportId"][i])
            rr = int(np.searchsorted(dist, n, side="right") - 1)
            k = n - dist[rr]
            if rr != r or k >= owned[r] or realOffset[rr] + k != int(d["gid"][i]):
                wrong += 1
    c.ok(wrong == 0, "%s: %d ids decoding to the wrong oxley node" % (label, wrong))

    return len(allPos)


def report():
    counts = sorted({int(os.path.basename(f).split("_")[0][1:])
                     for f in glob.glob(os.path.join(OUT, "p*_c*_r*.npz"))})
    if not counts:
        print("no dumps found in %s" % OUT)
        return 1
    c = Checker()
    print("rank counts found: %s\n" % counts)
    for ci, (name, _) in enumerate(CASES):
        positions = {}
        for p in counts:
            files = sorted(glob.glob(os.path.join(OUT, "p%d_c%d_r*.npz" % (p, ci))))
            ranks = [np.load(f) for f in files]
            if not ranks:
                continue
            positions[p] = check_case(c, "%s @%dr" % (name, p), ranks)
        same = len(set(positions.values())) <= 1
        c.ok(same, "%s: node count varies with rank count %s" % (name, positions))
        print("  %-20s nodes %-28s %s"
              % (name, positions, "OK" if same else "*** MISMATCH ***"))

    print("\nRESULT: %s" % ("ALL CHECKS PASSED" if not c.failures else
                            "FAILED\n  " + "\n  ".join(c.failures)))
    return 0 if not c.failures else 1


if __name__ == "__main__":
    if "--report" in sys.argv:
        sys.exit(report())
    n = int(sys.argv[sys.argv.index("--dump") + 1]) if "--dump" in sys.argv else 1
    dump(n)
