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
Local mesh refinement on a 3D oxley domain with a RefinementQueue3D.

The 3D counterpart of refinement_queue.py. A queue collects refinement
operations and applies them to a Brick. apply() does NOT change the domain it
is given: it returns a NEW, refined domain.

The example solves, on the finley export of each mesh, the Poisson problem

    -div(grad(u)) = f   in [0,1]^3,   u = 0 on the boundary,

with a source f concentrated in a small ball, once on a coarse mesh and once
on a mesh refined around the source and along the top face. The solutions are
saved in the subdirectory "data".
"""
__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.weipa import saveVTK
try:
    # finley must be imported for toFinley() to hand back a usable domain
    import esys.finley
    from esys.oxley import Brick, RefinementQueue3D, fromFinleyData
    HAVE_OXLEY = True
except ImportError:
    HAVE_OXLEY = False

# centre and radius of the source
XC, YC, ZC, R = 0.3, 0.6, 0.4, 0.1


def sourceMask(domain):
    """positive where the source is, a mask for the elements to refine"""
    x = Function(domain).getX()
    return whereNegative(length(x - [XC, YC, ZC]) - R)


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
           + whereZero(x[1]) + whereZero(x[1]-1.) \
           + whereZero(x[2]) + whereZero(x[2]-1.)
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

if not HAVE_OXLEY:
    print("This example requires the oxley domain.")
else:
    # the coarse domain: 4 x 4 x 4 blocks, each subdivided once
    # -> 8 x 8 x 8 elements
    coarse = Brick(n0=4, n1=4, n2=4, l0=1., l1=1., l2=1., refine_level=1)

    # collect the refinements. A level is ABSOLUTE: an element is subdivided
    # until it reaches that level. In 3D, "top" and "bottom" are the faces
    # normal to x2, "north" ("back") and "south" ("front") those normal to x1,
    # and "east" ("right") and "west" ("left") those normal to x0.
    queue = RefinementQueue3D()
    queue.setRefinementLevel(3)
    queue.refineSphere(x0=XC, y0=YC, z0=ZC, r=2*R)      # around the source, level 3
    queue.refineMask("source", level=4)                 # the source itself, level 4
    queue.refineBorder(border="top", dx=0.1, level=2)   # a layer below x2=1
    queue.print()

    # apply() returns a NEW domain; coarse is not modified. The mask for the
    # tag "source" is passed as a keyword argument and must live on coarse.
    fine = queue.apply(coarse, source=sourceMask(coarse))

    u_coarse = solve(coarse)
    u_fine = solve(fine)

    print("coarse mesh: %d elements, max(u) = %e" % (numElements(coarse), Lsup(u_coarse)))
    print("refined mesh: %d elements, max(u) = %e" % (numElements(fine), Lsup(u_fine)))

    # save the solutions into the subdirectory "data" for visualisation, e.g.
    # with VisIt or ParaView. mkDir is MPI safe: only one rank creates it.
    mkDir("data")
    saveVTK(os.path.join("data", "refinement3d_coarse.vtu"), u=u_coarse)
    saveVTK(os.path.join("data", "refinement3d_fine.vtu"), u=u_fine)
