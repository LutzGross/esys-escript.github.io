/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#ifndef __OXLEY_FINLEYCONVERTER_H__
#define __OXLEY_FINLEYCONVERTER_H__

#include <oxley/OxleyDomain.h>

#include <escript/AbstractDomain.h>

namespace oxley {

/**
   \brief
   Builds a finley domain describing the same mesh as an oxley forest.

   The forest supplies the geometry and topology; finley owns the degrees of
   freedom, the parallel overlap and everything a PDE needs. Each rank hands
   over the part of the forest it holds, with nodes named by their lnodes global
   id, and finley resolves the references and distributes the mesh, so no
   communication happens here beyond the conformity check.

   Every octant is split into simplices - Tri3 in 2D, Tet4 in 3D - so that one
   element family serves both conforming and adaptive forests. A quad face is
   always cut by the diagonal through its lowest global node id, which both
   sides of the face compute independently and agree on, so the split needs no
   communication.

   Currently restricted to CONFORMING forests: the split of a hanging face
   needs global ids for the hanging positions, which lnodes does not give them.
   Converting a forest with hanging nodes throws.

   \param dom the forest to convert
   \param order integration order (1 or 2)
   \param reducedOrder reduced integration order (1 or 2)
   \param optimize whether to let finley repartition with ParMETIS. Pass false
                   to keep the partition the forest already has.
   \param simplices when false, emit one Rec4/Hex8 per octant instead of
                    splitting. Only for telling a fault in the finley handover
                    apart from a fault in the split - it is not a production
                    path, since a forest would then change element family the
                    moment refinement introduces a hanging node.
*/
escript::Domain_ptr toFinley(const OxleyDomain& dom, int order = -1,
                             int reducedOrder = -1, bool optimize = false,
                             bool simplices = true);

} // namespace oxley

#endif // __OXLEY_FINLEYCONVERTER_H__
