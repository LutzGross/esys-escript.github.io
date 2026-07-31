"""
Which output slots does no rank fill?

weipa writes one shared point list: rank r emits the local nodes whose
nodeOutputIndex falls in [outputDistribution[r], outputDistribution[r+1]), and
every rank's connectivity indexes that list. A slot nobody claims is written as
a zero, which is what validate_output_hanging reports as an "unfilled sample".

Dumps per rank, then a serial pass works out, for every slot in the global list,
which ranks claim it - none, one, or several.

    ./bin/run-escript -n2 -t1 oxley-designs/diag_output_slots.py
    python3 oxley-designs/diag_output_slots.py --report
"""
import glob
import os
import sys

import numpy as np

OUT = os.environ.get("OXLEY_SLOT_DIR", "/tmp/oxley_slots")


def dump():
    import esys.escript as esc
    import esys.oxley as oxley

    rank = esc.getMPIRankWorld()
    if rank == 0:
        os.makedirs(OUT, exist_ok=True)
    esc.MPIBarrierWorld()

    dom = oxley.Rectangle(n0=2, n1=1, refine_level=[[2], [3]])
    info = dom.getMeshInfo(True)
    np.savez(os.path.join(OUT, "rank%d.npz" % rank),
             coords=info["nodeCoords"].reshape(-1, 2),
             gni=info["nodeOutputIndex"],
             dist=info["outputDistribution"],
             cn=info["constrainedNodes"],
             owner=info["constrainedOwner"],
             masters=info["constraintMasters"],
             numReal=info["numRealNodes"],
             numOwned=info["numOwnedNodes"])
    print("rank %d: %d nodes, %d real, %d owned, %d constrained"
          % (rank, len(info["nodeOutputIndex"]), info["numRealNodes"],
             info["numOwnedNodes"], len(info["constrainedNodes"])))


def report():
    files = sorted(glob.glob(os.path.join(OUT, "rank*.npz")))
    ranks = [np.load(f) for f in files]
    nr = len(ranks)
    dist = ranks[0]["dist"]
    total = int(dist[nr])
    print("ranks=%d  global slots=%d  distribution=%s" % (nr, total, list(dist)))

    claims = [[] for _ in range(total)]
    for r, d in enumerate(ranks):
        lo, hi = int(d["dist"][r]), int(d["dist"][r + 1])
        for i, g in enumerate(d["gni"]):
            g = int(g)
            if lo <= g < hi:                      # this rank writes it
                claims[g].append((r, i))

    unclaimed = [s for s in range(total) if not claims[s]]
    multi = [s for s in range(total) if len(claims[s]) > 1]
    print("slots nobody writes : %d" % len(unclaimed))
    print("slots >1 rank writes: %d" % len(multi))

    for s in unclaimed[:12]:
        owner = int(np.searchsorted(dist, s, side="right") - 1)
        # who HOLDS a node at that slot, even if outside its own range?
        holders = []
        for r, d in enumerate(ranks):
            for i, g in enumerate(d["gni"]):
                if int(g) == s:
                    kind = "real" if i < int(d["numReal"]) else "hanging"
                    holders.append("r%d[%d]%s%s" % (r, i, kind,
                                   "" if i < int(d["numReal"]) else
                                   " owner=%d" % int(d["owner"][
                                       list(d["cn"]).index(i)])))
        print("  slot %4d in block of rank %d, held by: %s"
              % (s, owner, ", ".join(holders) if holders else "NOBODY"))

    # is every local node's slot inside some block, and monotone within a block?
    for r, d in enumerate(ranks):
        lo, hi = int(d["dist"][r]), int(d["dist"][r + 1])
        own = [int(g) for g in d["gni"] if lo <= int(g) < hi]
        mono = all(own[k] < own[k + 1] for k in range(len(own) - 1))
        print("  rank %d writes %d slots, monotone in local order: %s"
              % (r, len(own), mono))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        dump()
