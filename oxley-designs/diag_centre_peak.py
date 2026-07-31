"""
Bisect the "centre peak 5x5" case, which stalls at 4 ranks but not at 2 or 3.

Each rank appends to its OWN file and flushes after every stage, so if this is a
collective reached by some ranks and not others, the last line of each file says
which ranks got where. Nothing is printed inside a `if rank == 0` guard - a
collective in there is itself a deadlock, and would hide the one being hunted.
"""
import os
import sys
import time

import esys.escript as esc
import esys.finley                     # registers the concrete domain type
import esys.oxley as oxley
from esys.escript import Function, FunctionOnBoundary, integrate

OUT = os.environ.get("OXLEY_DIAG_DIR", "/tmp/oxley_peak")
LEVELS = [[0, 0, 1, 0, 0],
          [0, 1, 2, 1, 0],
          [1, 2, 3, 2, 1],
          [0, 1, 2, 1, 0],
          [0, 0, 1, 0, 0]]

rank = esc.getMPIRankWorld()
size = esc.getMPISizeWorld()
if rank == 0:
    os.makedirs(OUT, exist_ok=True)
esc.MPIBarrierWorld()

log = open(os.path.join(OUT, "rank%d.log" % rank), "w")
t0 = time.time()


def stage(name):
    log.write("%8.2fs  %s\n" % (time.time() - t0, name))
    log.flush()
    os.fsync(log.fileno())


stage("start, size=%d" % size)

dom = oxley.Rectangle(n0=5, n1=5, l0=5.0, l1=5.0, refine_level=LEVELS)
stage("Rectangle built")

conf = dom.isConforming()               # collective
stage("isConforming -> %s" % conf)

info = dom.getMeshInfo(True)            # getMeshAccess: seams + numbering
stage("getMeshInfo: %d nodes, %d elements, %d hanging"
      % (info["numNodes"], info["numElements"],
         sum(1 for h in info["elementFaceHangingNode"] if h >= 0)))

fin = dom.toFinley()
stage("toFinley done")

v = integrate(esc.Scalar(1, Function(fin)))
stage("volume = %.12g" % v)

a = integrate(esc.Scalar(1, FunctionOnBoundary(fin)))
stage("boundary = %.12g" % a)

n = FunctionOnBoundary(fin).getNormal()
xb = FunctionOnBoundary(fin).getX()
stage("int(n.x) = %.12g" % integrate(esc.inner(n, xb)))

stage("ALL STAGES COMPLETE")
log.close()
sys.exit(0)
