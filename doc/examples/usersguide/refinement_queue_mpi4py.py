##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################
"""
Local mesh refinement on MPI sub-communicators.

The MPI ranks are split into groups with mpi4py. Each group builds the same
coarse oxley domain on its own sub-communicator, refines it around a source
to a different level with a RefinementQueue2D and solves a Poisson problem
on the result. The source itself is refined one level further through a
mask. The groups work at the same time; world rank 0 collects and
prints what they found. Each group saves its solution as
data/refinement_group<k>.vtu.

A RefinementQueue has no communicator of its own: apply() copies the domain
it is given, so the refined domain lives on the same communicator as that
domain. To refine on a sub-communicator, build the domain with comm=...
A mask refinement only names its mask in the queue; the mask is handed to
apply() and is defined on the domain being refined, so it lives on the same
sub-communicator. Each rank holds the mask on its own elements only: apply()
gathers the marked elements over the communicator, so every rank of the
group refines the same ones.

Usage:
    run-escript -n 4 refinement_queue_mpi4py.py
"""
__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os
from mpi4py import MPI
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.weipa import saveVTK
try:
    # finley must be imported for toFinley() to hand back a usable domain
    import esys.finley
    from esys.oxley import Rectangle, RefinementQueue2D, fromFinleyData
    HAVE_OXLEY = True
except ImportError:
    HAVE_OXLEY = False

# centre and radius of the source
XC, YC, R = 0.3, 0.6, 0.05
# the coarse mesh is at level 2; group k refines around the source to level
# 3+k, and the source itself to level 4+k
BASE_LEVEL = 2
MAX_GROUPS = 3


def sourceMask(domain):
    """positive where the source is, a mask for the elements to refine"""
    x = Function(domain).getX()
    return whereNegative(length(x - [XC, YC]) - R)


def solve(domain):
    """
    solves -div(grad(u)) = f, u = 0 on the boundary, f = 1 in the source, and
    returns u on domain. PDEs on an oxley domain are solved on its finley
    export: toFinley() turns the forest, hanging nodes included, into a
    conforming finley mesh, and fromFinleyData() brings the solution back.
    """
    fin = domain.toFinley()
    x = fin.getX()
    pde = LinearPDE(fin)
    pde.setSymmetryOn()
    gammaD = whereZero(x[0]) + whereZero(x[0]-1.) \
           + whereZero(x[1]) + whereZero(x[1]-1.)
    pde.setValue(A=kronecker(fin), Y=sourceMask(fin), q=gammaD)
    return fromFinleyData(pde.getSolution(), domain)


def numElements(domain):
    """
    the global number of elements of domain. getMeshInfo() gives the elements
    this rank owns - each element belongs to exactly one rank, whereas the
    samples of a Function also count the overlap shared with neighbouring
    ranks - and they are summed over the communicator the domain lives on,
    which need not be all ranks. Without mpi4py there is no communicator to
    hand back, but then the domain lives on all ranks.
    """
    local = int(domain.getMeshInfo()["numElements"])
    comm = domain.getMPIComm()
    if comm is None:
        return getMPIWorldSum(local)
    return comm.allreduce(local)

world = MPI.COMM_WORLD
if not HAVE_OXLEY:
    if world.rank == 0:
        print("This example requires the oxley domain.")
else:
    # create the output directory before the ranks go their separate ways:
    # mkDir is collective over all ranks
    mkDir("data")

    ngroups = min(world.size, MAX_GROUPS)
    color = world.rank % ngroups
    sub = world.Split(color, world.rank)
    level = BASE_LEVEL + 1 + color

    # the domain lives on the sub-communicator ...
    coarse = Rectangle(n0=4, n1=4, l0=1., l1=1., refine_level=BASE_LEVEL, comm=sub)

    # ... and so do the mask built on it and the refined domain apply() returns
    queue = RefinementQueue2D()
    queue.refineCircle(x0=XC, y0=YC, r=2*R, level=level)
    queue.refineMask("source", level=level+1)
    fine = queue.apply(coarse, source=sourceMask(coarse))

    u = solve(fine)
    # Lsup and integrate reduce over the domain's communicator, i.e. the group
    result = (color, sub.size, level, numElements(fine), Lsup(u),
              integrate(u))

    # every rank of the group takes part in writing the group's file
    saveVTK(os.path.join("data", "refinement_group%d.vtu" % color), u=u)

    # one report per group, collected on world rank 0
    results = world.gather(result if sub.rank == 0 else None, root=0)
    if world.rank == 0:
        print("%5s %5s %5s %9s %13s %13s" % ("group", "ranks", "level",
                                            "elements", "max(u)", "integral(u)"))
        for r in sorted(x for x in results if x is not None):
            print("%5d %5d %5d %9d %13.6e %13.6e" % r)

    sub.Free()
