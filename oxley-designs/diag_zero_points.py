"""
List the points a parallel .vtu writes as zero for a field that is nowhere zero.

    ./bin/run-escript -n2 -t1 oxley-designs/diag_zero_points.py
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

import esys.escript as esc
import esys.oxley as oxley
from esys.escript import ContinuousFunction
from esys.weipa import saveVTK

OUT = "/tmp/oxley_zero"
rank = esc.getMPIRankWorld()
if rank == 0:
    os.makedirs(OUT, exist_ok=True)
esc.MPIBarrierWorld()

dom = oxley.Rectangle(n0=2, n1=1, refine_level=[[2], [3]])
x = ContinuousFunction(dom).getX()
lin = 1.0 + 2.0 * x[0] + 3.0 * x[1]          # nowhere zero on this domain
path = os.path.join(OUT, "z.vtu")
saveVTK(path, lin=lin)
esc.MPIBarrierWorld()

if rank != 0:
    sys.exit(0)

root = ET.parse(path).getroot()
# the Points array carries no Name, so find it by its parent element
pts = val = None
for parent in root.iter():
    for da in parent:
        if not da.tag.endswith("DataArray"):
            continue
        if parent.tag.endswith("Points"):
            pts = np.array([float(v) for v in da.text.split()]).reshape(-1, 3)
        elif da.get("Name") == "lin":
            val = np.array([float(v) for v in da.text.split()])
if pts is None or val is None:
    raise SystemExit("could not find Points / lin in %s" % path)
want = 1.0 + 2.0 * pts[:, 0] + 3.0 * pts[:, 1]
bad = np.nonzero(np.abs(val - want) > 1e-5)[0]

print("points=%d  wrong=%d" % (len(val), len(bad)))
for i in bad:
    onseam = abs(pts[i, 0] - 0.5) < 1e-12
    ismid = (abs(pts[i, 1] * 8 - round(pts[i, 1] * 8)) < 1e-9)
    print("  point %3d at [%7.4f %7.4f]  got %8.4f want %8.4f  %s"
          % (i, pts[i, 0], pts[i, 1], val[i], want[i],
             "ON THE x=0.5 SEAM" if onseam else ""))
