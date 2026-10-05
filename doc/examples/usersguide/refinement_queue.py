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
Local mesh refinement on an oxley domain with a RefinementQueue.

A queue collects refinement operations and applies them to a domain. apply()
does NOT change the domain it is given: it returns a NEW, refined domain. The
coarse domain, and any Data defined on it, stay valid.

A queue is a template which is not tied to a domain. A refinement driven by
data, refineMask, therefore only names its mask; the mask itself, defined on
the domain being refined, is handed to apply() under that name.

The example solves, on the finley export of each mesh, the Poisson problem

    -div(grad(u)) = f   in [0,1]^2,   u = 0 on the boundary,

with a source f concentrated in a small disc, once on a coarse mesh and once
on a mesh refined around the source and along the top border.
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
    from esys.oxley import Rectangle, RefinementQueue2D, fromFinleyData
    HAVE_OXLEY = True
except ImportError:
    HAVE_OXLEY = False

# centre and radius of the source
XC, YC, R = 0.3, 0.6, 0.05


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


def sourceMask(domain):
    """positive where the source is, a mask for the elements to refine"""
    x = Function(domain).getX()
    return whereNegative(length(x - [XC, YC]) - R)


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
    # the coarse domain: 4 x 4 blocks, each subdivided twice -> 16 x 16 elements
    coarse = Rectangle(n0=4, n1=4, l0=1., l1=1., refine_level=2)

    # collect the refinements. setRefinementLevel sets the default level for
    # operations which do not give one. A level is ABSOLUTE: an element is
    # subdivided until it reaches that level, so elements which are already
    # finer are left alone.
    queue = RefinementQueue2D()
    queue.setRefinementLevel(4)
    queue.refineCircle(x0=XC, y0=YC, r=2*R)             # around the source, level 4
    queue.refineMask("source", level=5)                 # the source itself, level 5
    queue.refineBorder(border="top", dx=0.1, level=3)   # a strip along y=1
    queue.refineRegion(x0=0.7, y0=0.0, x1=0.9, y1=0.2)  # a box, the default level
    queue.print()

    # apply() returns a NEW domain; coarse is not modified. The mask for the
    # tag "source" is passed as a keyword argument and must live on coarse.
    fine = queue.apply(coarse, source=sourceMask(coarse))

    u_coarse = solve(coarse)
    u_fine = solve(fine)

    print("coarse mesh: %d elements, max(u) = %e" % (numElements(coarse), Lsup(u_coarse)))
    print("refined mesh: %d elements, max(u) = %e" % (numElements(fine), Lsup(u_fine)))

    # u_coarse lives on the coarse domain and is still valid. Data on two
    # different domains cannot be combined: u_coarse + u_fine raises an error.

    # a queue can be applied again, to any 2D oxley domain, with the mask on
    # that domain:
    other = Rectangle(n0=8, n1=8, l0=1., l1=1., refine_level=1)
    finer = queue.apply(other, source=sourceMask(other))
    print("queue applied to an 8 x 8 block domain: %d elements" % numElements(finer))

    # save the solutions into the subdirectory "data" for visualisation, e.g.
    # with VisIt or ParaView. mkDir is MPI safe: only one rank creates it.
    mkDir("data")
    saveVTK(os.path.join("data", "refinement_coarse.vtu"), u=u_coarse)
    saveVTK(os.path.join("data", "refinement_fine.vtu"), u=u_fine)
