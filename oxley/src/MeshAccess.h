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

#include <array>
#include <cmath>
#include <map>
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

    /// corners per boundary face (2 in 2D, 4 in 3D)
    int nodesPerFace = 0;
    /// local boundary faces on this rank
    long numFaces = 0;
    /// per-face local node indices, size numFaces*nodesPerFace. Unlike
    /// elementNodes these are NOT in z-order: they are wound so that the
    /// right-hand rule yields the OUTWARD normal (in 2D, the domain lies to the
    /// left of the directed edge). Consumers that compute normals depend on it.
    std::vector<long> faceNodes;
    /// per-face tag, size numFaces (left/right/bottom/top[/front/back])
    std::vector<long> faceTags;
    /// per-face index of the element the face belongs to, size numFaces
    std::vector<long> faceElements;

    // ------------------------------------------------------------------------
    // Materialised hanging positions (getMeshAccess(materializeHanging=true)).
    //
    // A degree-1 lnodes does not number a hanging position: the corresponding
    // element_nodes slot holds a MASTER instead, so an element whose corner hangs
    // has no node of its own there. That is fine for assembly, which resolves the
    // constraint, but not for a consumer that needs one node per corner - weipa
    // draws a cell from its corner list, and a master in that list is a corner in
    // the wrong place.
    //
    // With materializeHanging the hanging positions become real nodes, appended
    // after the lnodes ones: indices [numRealNodes, numNodes) and the tail of
    // nodeCoords/nodeGlobalId. They carry NO sample of the domain's nodal
    // function space, so a consumer with per-node values must fill them from the
    // constraint below - a hanging node's value is the average of its masters.
    // ------------------------------------------------------------------------

    /// nodes that come from lnodes; equals numNodes unless hanging positions
    /// were materialised, in which case the rest are the materialised ones
    long numRealNodes = 0;
    /// masters per materialised node: 2 for an edge midpoint, 4 for a 3D face
    /// centre. Shorter lists are padded with master -1 and weight 0.
    int mastersPerConstrainedNode = 0;
    /// local index of each materialised node, size = numNodes - numRealNodes
    std::vector<long> constrainedNodes;
    /// masters of each materialised node as local node indices, size
    /// constrainedNodes.size()*mastersPerConstrainedNode
    std::vector<long> constraintMasters;
    /// weight of each master, same size and layout as constraintWeights
    std::vector<double> constraintWeights;

    // ------------------------------------------------------------------------
    // Output numbering.
    //
    // nodeGlobalId is the lnodes numbering: unique, but a rank's owned ids are
    // NOT one contiguous block once the materialised nodes are added, and the
    // materialised ones sit far above the rest. A writer that emits one shared
    // list of points (VTK does; Silo writes a block per rank and does not care)
    // needs a numbering where each rank owns exactly one contiguous range, so
    // that every node is written by exactly one rank and the connectivity of
    // every rank indexes the same list.
    //
    // nodeOutputIndex is that numbering: rank r owns
    // [outputDistribution[r], outputDistribution[r+1]), its lnodes-owned nodes
    // first, then its materialised ones.
    // ------------------------------------------------------------------------

    /// contiguous global output index of each local node, size numNodes
    std::vector<long> nodeOutputIndex;
    /// first output index of each rank, size mpiSize+1; the last entry is the
    /// global number of output nodes
    std::vector<long> outputDistribution;
};

/// Materialised hanging nodes already created, keyed by position quantised to
/// 1e-9, so the fine elements meeting at one hanging position share a node.
typedef std::map<std::array<long,3>, long> HangingNodeMap;

/**
   \brief
   Appends a node at a hanging position, or returns the one already there.

   \param m the mesh being built; numNodes, nodeCoords, nodeGlobalId and the
            constraint arrays all grow by one when a node is created. Global ids
            are left at -1 for OxleyDomain::assignHangingNodeIds().
   \param seen positions already materialised
   \param xyz the position (three components; the third is ignored in 2D)
   \param masters local node indices this position is the average of
   \param weights their weights, summing to 1
   \param n how many masters, at most m.mastersPerConstrainedNode
   \return the local index of the node at that position
*/
inline long addHangingNode(MeshAccess& m, HangingNodeMap& seen,
                           const double xyz[3], const long* masters,
                           const double* weights, int n)
{
    std::array<long,3> key;
    for (int d = 0; d < 3; ++d)
        key[d] = (long) std::llround(xyz[d] * 1e9);
    HangingNodeMap::const_iterator it = seen.find(key);
    if (it != seen.end())
        return it->second;

    const long idx = m.numNodes++;
    for (int d = 0; d < m.numDim; ++d)
        m.nodeCoords.push_back(xyz[d]);
    m.nodeGlobalId.push_back(-1);
    m.constrainedNodes.push_back(idx);
    for (int k = 0; k < m.mastersPerConstrainedNode; ++k) {
        m.constraintMasters.push_back(k < n ? masters[k] : -1);
        m.constraintWeights.push_back(k < n ? weights[k] : 0.);
    }
    seen[key] = idx;
    return idx;
}

} // namespace oxley

#endif // __OXLEY_MESHACCESS_H__
