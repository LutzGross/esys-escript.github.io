"""
Write Silo files for a 2D oxley forest and for the finley mesh it converts to,
for side-by-side inspection in VisIt.

    ./bin/run-escript -n1 -t1 oxley-designs/silo_2d.py
    ./bin/run-escript -n4 -t1 oxley-designs/silo_2d.py        # see the partition
    ./bin/run-escript -n1 -t1 oxley-designs/silo_2d.py 3      # refine level 3

Produces, in the current directory:
    oxley_2d.silo    the forest as quadrilaterals
    finley_2d.silo   the same mesh split into Tri3

Fields on both, so they can be compared directly:
    coords   node positions as a vector
    f        a smooth nodal field, sin(pi x) sin(pi y)
    rank     the MPI rank that owns each cell - shows the partition
and on the finley mesh only:
    sol      solution of -div(grad u) = 1, u = 0 on the left and right edges

In VisIt open both files, use Controls > Subset or the Mesh plot to see the
element shapes, and Pseudocolor on `rank` to see the decomposition.

Set REFINE to a nested list to build a forest with hanging nodes, e.g.
[[3,1],[1,2]] for a 2x2 block layout. That still writes oxley_2d.silo, but the
conversion is refused for now - the simplex split does not yet number the
hanging positions.
"""
import math
import sys

try:
    import esys.escript as esc
    from esys.escript import ContinuousFunction, Function, sup, whereZero
    from esys.escript.linearPDEs import LinearPDE, SolverOptions
    import esys.finley as finley  # REQUIRED: registers the finley domain type,
                                  # otherwise toFinley() hands back a bare Domain
    import esys.oxley as oxley
    from esys.weipa import saveSilo
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
L0, L1 = 1.0, 1.0

if len(sys.argv) > 1:
    REFINE = int(sys.argv[1])

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

saveSilo("oxley_2d.silo", **common_fields(ox))
say("wrote oxley_2d.silo")

conforming = ox.isConforming()          # collective: every rank must call it
if not conforming:
    say("forest has hanging nodes, so it cannot be converted yet - "
        "only oxley_2d.silo was written")
    sys.exit(0)

fin = ox.toFinley()
fields = common_fields(fin)
fields["sol"] = poisson(fin)
saveSilo("finley_2d.silo", **fields)
say("wrote finley_2d.silo  (Tri3, 2 triangles per quad)")
