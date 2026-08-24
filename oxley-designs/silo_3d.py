
"""
Write Silo files for a 3D oxley forest and for the finley mesh it converts to,
for side-by-side inspection in VisIt.

    ./bin/run-escript -n1 -t1 oxley-designs/silo_3d.py
    ./bin/run-escript -n4 -t1 oxley-designs/silo_3d.py 1 _n4  # see the partition
    ./bin/run-escript -n1 -t1 oxley-designs/silo_3d.py 2      # refine level 2
    ./bin/run-escript -n1 -t1 oxley-designs/silo_3d.py "[[[2,1],[1,1]],[[1,1],[1,2]]]"

The last form gives each block its own level, which leaves a 2:1 seam - hanging
nodes - along the block boundaries; the block layout is taken from the list.

Produces, in the current directory (TAG is the second argument, and defaults to
"_hanging" for a per-block refinement so it does not overwrite the uniform run):
    oxley_3dTAG.silo    the forest as hexahedra
    finley_3dTAG.silo   the same mesh split into Tet4, six per octant

Fields on both, so they can be compared directly:
    coords   node positions as a vector
    f        a smooth nodal field, sin(pi x) sin(pi y) sin(pi z)
    rank     the MPI rank that owns each cell - shows the partition
and on the finley mesh only:
    sol      solution of -div(grad u) = 1, u = 0 on the left and right faces

The tetrahedra are worth a look on their own: each hexahedron is coned from its
lowest-id corner over the three faces not containing it, so the six tets are not
all congruent and the pattern varies from cell to cell. Clip or use the Mesh
plot with "wireframe" to see it.

Set REFINE to a nested list to build a forest with hanging nodes, e.g.
[[[2,1],[1,1]],[[1,1],[1,2]]] for a 2x2x2 block layout. That still writes
oxley_3d.silo, but the conversion is refused for now - the simplex split does
not yet number the hanging positions.
"""
import ast
import math
import sys

try:
    import esys.escript as esc
    from esys.escript import ContinuousFunction, Function, whereZero
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
N0, N1, N2 = 2, 2, 2   # blocks (p8est trees)
REFINE = 1             # uniform level, or a nested list per block for hanging
# The domain extent follows the block layout so that every BLOCK is a unit cube
# and hence every element is a cube. With a fixed 1x1x1 domain a 2x1x1 block
# layout would give 0.5x1x1 blocks and every element would inherit that aspect -
# correct, but hard to read in a viewer.
L0, L1, L2 = float(N0), float(N1), float(N2)

# argv[1]: the refinement, either an int for a uniform level or a nested list
#          giving one level per block, e.g. "[[[2,1],[1,1]],[[1,1],[1,2]]]" -
#          blocks at different levels leave a 2:1 seam, i.e. hanging nodes, along
#          their boundary.
# argv[2]: a tag appended to the file names, so several runs can coexist
#          (defaults to "_hanging" for a per-block refinement, "" otherwise).
if len(sys.argv) > 1:
    REFINE = ast.literal_eval(sys.argv[1])
    if isinstance(REFINE, list):
        N0, N1, N2 = len(REFINE), len(REFINE[0]), len(REFINE[0][0])
        L0, L1, L2 = float(N0), float(N1), float(N2)
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
        "f": (esc.sin(math.pi * x[0]) * esc.sin(math.pi * x[1])
              * esc.sin(math.pi * x[2])),
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


say("building oxley Brick(%d, %d, %d, refine_level=%s)" % (N0, N1, N2, REFINE))
ox = oxley.Brick(n0=N0, n1=N1, n2=N2, l0=L0, l1=L1, l2=L2, refine_level=REFINE)

fields = common_fields(ox)
saveSilo("oxley_3d%s.silo" % TAG, **fields)
saveVTK("oxley_3d%s.vtu" % TAG, **fields)
say("wrote oxley_3d%s.silo and .vtu" % TAG)

conforming = ox.isConforming()          # collective: every rank must call it
if not conforming:
    say("forest has hanging nodes, so it cannot be converted yet - "
        "only oxley_3d%s.silo was written" % TAG)
    sys.exit(0)

fin = ox.toFinley()
fields = common_fields(fin)
fields["sol"] = poisson(fin)
saveSilo("finley_3d%s.silo" % TAG, **fields)
saveVTK("finley_3d%s.vtu" % TAG, **fields)
say("wrote finley_3d%s.silo  (Tet4, six tetrahedra per octant)" % TAG)
