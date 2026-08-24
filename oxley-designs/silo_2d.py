"""
Write Silo files for a 2D oxley forest and for the finley mesh it converts to,
for side-by-side inspection in VisIt.

    ./bin/run-escript -n1 -t1 oxley-designs/silo_2d.py
    ./bin/run-escript -n4 -t1 oxley-designs/silo_2d.py 2 _n4  # see the partition
    ./bin/run-escript -n1 -t1 oxley-designs/silo_2d.py 3      # refine level 3
    ./bin/run-escript -n1 -t1 oxley-designs/silo_2d.py "[[3,1],[1,2]]"

The last form gives each block its own level, which leaves a 2:1 seam - hanging
nodes - along the block boundaries; the block layout is taken from the list.

Produces, in the current directory (TAG is the second argument, and defaults to
"_hanging" for a per-block refinement so it does not overwrite the uniform run):
    oxley_2dTAG.silo    the forest as quadrilaterals
    finley_2dTAG.silo   the same mesh split into Tri3 (serial also for an
                        adaptive forest: a hanging node becomes an ordinary
                        vertex and the coarse element splits into 3-6 triangles)

Fields on both, so they can be compared directly:
    coords   node positions as a vector
    f        a smooth nodal field, sin(pi x) sin(pi y)
    rank     the MPI rank that owns each cell - shows the partition
and on the finley mesh only:
    sol      solution of -div(grad u) = 1, u = 0 on the left and right edges

In VisIt open both files, use Controls > Subset or the Mesh plot to see the
element shapes, and Pseudocolor on `rank` to see the decomposition.

FACE ELEMENTS. The .silo carries every mesh, so the boundary elements are always
in it. A .vtu holds ONE mesh, so weipa writes one file per mesh and only for
meshes that carry a variable: with volume/nodal fields alone the boundary is
simply not in the VTK output. Passing a field on FunctionOnBoundary (say
dom.getNormal()) writes it as <prefix>_FaceElements.vtu - but note that this also
renames the volume file to <prefix>_Elements.vtu. finley behaves identically.

Set REFINE to a nested list to build a forest with hanging nodes, e.g.
[[3,1],[1,2]] for a 2x2 block layout. Comparing oxley_2dTAG.silo with
finley_2dTAG.silo then shows the point of the converter: the quad carrying a
hanging node on its edge becomes a fan of triangles through that node, so the
T-junction is resolved by the triangulation rather than by a constraint.
"""
import ast
import math
import sys

try:
    import esys.escript as esc
    from esys.escript import ContinuousFunction, Function, sup, whereZero
    from esys.escript.linearPDEs import LinearPDE, SolverOptions
    import esys.finley as finley  # REQUIRED: registers the finley domain type,
                                  # otherwise toFinley() hands back a bare Domain
    import esys.oxley as oxley
    from esys.weipa import saveSilo, saveVTK
except ImportError as e:
    raise SystemExit(
        "cannot import esys (%s).\n\n"
        "Plain python3 works for a serial run, but needs the build on the path:\n"
        "    export PYTHONPATH=<escript root>\n"
        "    export LD_LIBRARY_PATH=<escript root>/lib/esys:<escript root>"
        "/esys.trilinos/lib\n"
        "Or run it through the launcher, which sets both:\n"
        "    <escript root>/bin/run-escript -n1 -t1 %s" % (e, __file__))

# ---------------------------------------------------------------- settings --
N0, N1 = 2, 2          # blocks (p4est trees)
REFINE = 2             # uniform level, or a nested list per block for hanging
# The domain extent follows the block layout so that every BLOCK is a unit
# square and hence every element is square. With a fixed 1x1 domain a 2x1 block
# layout would give 0.5x1.0 blocks and every element would inherit that 1:2
# aspect - correct, but hard to read in a viewer.
L0, L1 = float(N0), float(N1)

# argv[1]: the refinement, either an int for a uniform level or a nested list
#          giving one level per block, e.g. "[[3,1],[1,2]]" - blocks at different
#          levels leave a 2:1 seam, i.e. hanging nodes, along their boundary.
# argv[2]: a tag appended to the file names, so several runs can coexist
#          (defaults to "_hanging" for a per-block refinement, "" otherwise).
if len(sys.argv) > 1:
    REFINE = ast.literal_eval(sys.argv[1])
    if isinstance(REFINE, list):
        N0, N1 = len(REFINE), len(REFINE[0])
        L0, L1 = float(N0), float(N1)
TAG = sys.argv[2] if len(sys.argv) > 2 else ("_hanging" if isinstance(REFINE, list) else "")

rank = esc.getMPIRankWorld()


def say(msg):
    if rank == 0:
        print(msg, flush=True)


def common_fields(dom):
    """the fields written for both meshes, so they can be compared"""
    x = ContinuousFunction(dom).getX()
    owner = esc.Scalar(rank, Function(dom))
    owner.expand()                      # per cell, not a constant
    return {
        "coords": x,
        "f": esc.sin(math.pi * x[0]) * esc.sin(math.pi * x[1]),
        "rank": owner,
    }


def poisson(dom):
    pde = LinearPDE(dom, numEquations=1)
    pde.setSymmetryOn()
    x = dom.getX()
    pde.setValue(A=esc.kronecker(dom), Y=1.0,
                 q=whereZero(x[0]) + whereZero(x[0] - L0))
    pde.getSolverOptions().setTolerance(1e-12)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    return pde.getSolution()


say("building oxley Rectangle(%d, %d, refine_level=%s)" % (N0, N1, REFINE))
ox = oxley.Rectangle(n0=N0, n1=N1, l0=L0, l1=L1, refine_level=REFINE)

fields = common_fields(ox)
saveSilo("oxley_2d%s.silo" % TAG, **fields)
saveVTK("oxley_2d%s.vtu" % TAG, **fields)
say("wrote oxley_2d%s.silo and .vtu" % TAG)

# A non-conforming forest CAN be converted in 2D: the hanging position becomes an
# ordinary node of the finley mesh and the coarse element is split so that it is a
# vertex on both sides of the 2:1 seam. Under MPI that is not done yet - the node
# is materialised by the rank owning the COARSE octant of each seam, and both
# sides derive its global id from the same (octant, face) key, so the seam is
# resolved across ranks as well.
conforming = ox.isConforming()          # collective: every rank must call it

fin = ox.toFinley()
fields = common_fields(fin)
fields["sol"] = poisson(fin)
saveSilo("finley_2d%s.silo" % TAG, **fields)
saveVTK("finley_2d%s.vtu" % TAG, **fields)

# how the quads were split, so an odd-looking picture can be checked against the
# pattern table rather than guessed at: a quad with h hanging sides gives
# 2, 3, 4, 5 or 6 triangles depending on the configuration
nquad = esc.Scalar(1, Function(ox)).getNumberOfDataPoints() // 4    # 2x2 Gauss
ntri = esc.Scalar(1, Function(fin)).getNumberOfDataPoints() // 3    # Tri3
say("wrote finley_2d%s.silo  (Tri3, %d triangles from %d quads, %.2f per quad)"
    % (TAG, ntri, nquad, float(ntri) / max(nquad, 1)))
