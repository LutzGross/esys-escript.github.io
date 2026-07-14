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

#ifndef __OXLEY_MESHACCESS_H__
#define __OXLEY_MESHACCESS_H__

#include <vector>

namespace oxley {

/**
   \brief
   A p4est-independent view of the mesh, consumed by output (weipa) and, later,
   the assembler. It is the single public description of the mesh topology, so
   the node numbering scheme stays entirely inside the domain.

   Node numbering is lnodes-based: the local nodes are laid out as
   [ owned | ghost ]. An owned local node i has global id globalNodeOffset+i;
   ghost local nodes carry their explicit global id in nodeGlobalId. In serial,
   numOwnedNodes == numNodes, globalNodeOffset == 0 and global id == local id.

   Element corners are listed in p4est z-order (x fastest, then y, then z).
*/
struct MeshAccess
{
    int numDim = 0;             ///< spatial dimension (2 or 3)
    int nodesPerElement = 0;    ///< corners per element (4 in 2D, 8 in 3D)
    long numNodes = 0;          ///< local nodes (owned + ghost)
    long numOwnedNodes = 0;     ///< nodes owned by this rank
    long numElements = 0;       ///< local leaf elements
    long globalNodeOffset = 0;  ///< global id of this rank's first owned node

    /// node coordinates, size numNodes*numDim (node i at [i*numDim + d])
    std::vector<double> nodeCoords;
    /// global id of each local node, size numNodes
    std::vector<long> nodeGlobalId;
    /// per-element local node indices, size numElements*nodesPerElement (z-order)
    std::vector<long> elementNodes;
    /// per-element tag, size numElements
    std::vector<long> elementTags;
};

} // namespace oxley

#endif // __OXLEY_MESHACCESS_H__
