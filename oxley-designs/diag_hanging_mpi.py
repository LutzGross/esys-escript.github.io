"""
Why does the 2D hanging export fail under MPI but pass serially?

Each rank dumps its materialised mesh view to its own .npz - no collectives, so
this cannot deadlock - and a serial pass compares them. Two failure modes are
distinguished, and they need different fixes:

  (a) MISSING ON THE COARSE SIDE. Materialisation is driven by face_code, which
      p4est sets on the FINE element. A rank holding only the COARSE octant of a
      2:1 seam therefore never creates the midpoint node, so its element is split
      as a plain quad and the node is not a vertex of the coarse side. The row of
      the linear system for that node is then missing the coarse element's
      contribution.

  (b) DUPLICATED. The two fine elements of a seam can sit on different ranks.
      Each materialises the position independently and numbers it out of its OWN
      block in finaliseNodeNumbering(), so one physical point gets two global
      ids and the mesh is cracked there.

Run:
    ./bin/run-escript -n2 -t1 oxley-designs/diag_hanging_mpi.py
    python3 oxley-designs/diag_hanging_mpi.py --report
"""
import glob
import os
import sys

import numpy as np

OUT = "/tmp/claude-1000/-home-lgross-PycharmProjects-esys-escript-github-io/" \
      "e96cc7a5-2464-4e68-8ece-4733e0a85e05/scratchpad"
LEVELS = [[2], [3]]          # the one-seam case that fails on 2 ranks
N0, N1 = len(LEVELS), len(LEVELS[0])

# same quantisation the converter and addHangingNode() use to match positions
def key(p):
    return (int(round(p[0] * 1e9)), int(round(p[1] * 1e9)))


def dump():
    import esys.escript as esc
    import esys.oxley as oxley

    rank = esc.getMPIRankWorld()
    dom = oxley.Rectangle(n0=N0, n1=N1, l0=float(N0), l1=float(N1),
                          refine_level=LEVELS)
    info = dom.getMeshInfo(True)
    np.savez(os.path.join(OUT, "diag_rank%d.npz" % rank),
             coords=info["nodeCoords"].reshape(-1, 2),
             gid=info["nodeGlobalId"],
             en=info["elementNodes"].reshape(-1, info["nodesPerElement"]),
             numOwned=info["numOwnedNodes"])
    print("rank %d: %d nodes (%d owned), %d elements"
          % (rank, info["numNodes"], info["numOwnedNodes"], info["numElements"]))


def report():
    files = sorted(glob.glob(os.path.join(OUT, "diag_rank*.npz")))
    ranks = [np.load(f) for f in files]
    print("ranks: %d" % len(ranks))

    # every position anyone has a node at, and the global ids seen there
    idsAt = {}
    for r, d in enumerate(ranks):
        for i, p in enumerate(d["coords"]):
            idsAt.setdefault(key(p), set()).add(int(d["gid"][i]))

    dup = {k: v for k, v in idsAt.items() if len(v) > 1}
    print("\n(b) positions carrying MORE THAN ONE global id: %d" % len(dup))
    for k, v in sorted(dup.items())[:8]:
        print("      (%9.6f,%9.6f)  ids %s"
              % (k[0] / 1e9, k[1] / 1e9, sorted(v)))

    # (a): an element edge whose midpoint is a node SOMEWHERE, but not a node
    # this rank knows about - so this rank cannot insert it into the polygon
    zOrderEdges = [(0, 1), (1, 3), (3, 2), (2, 0)]   # c0 c1 c3 c2 boundary walk
    missing = 0
    for r, d in enumerate(ranks):
        local = {key(p) for p in d["coords"]}
        for e in d["en"]:
            for a, b in zOrderEdges:
                mid = 0.5 * (d["coords"][e[a]] + d["coords"][e[b]])
                k = key(mid)
                if k in idsAt and k not in local:
                    if missing < 8:
                        print("      rank %d element edge midpoint "
                              "(%9.6f,%9.6f) exists globally, not locally"
                              % (r, mid[0], mid[1]))
                    missing += 1
    print("\n(a) element edges whose hanging midpoint this rank cannot see: %d"
          % missing)

    print("\nverdict: %s"
          % ("both modes present" if dup and missing else
             "(b) duplicated ids only" if dup else
             "(a) coarse side blind only" if missing else
             "neither - look elsewhere"))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        dump()