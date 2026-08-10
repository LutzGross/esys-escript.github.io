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

   Local nodes are laid out as [ owned | ghost ], and everything here that names
   a node - elementNodes, faceNodes, the constraint arrays - uses that LOCAL
   index. On top of it sit three global numberings, which cannot be collapsed
   into one because their requirements conflict:

     nodeLnodesId   p4est's own numbering, with a block above it for the
                    materialised hanging positions, which lnodes does not name.
                    This is the id space escript's node samples live in.

     nodeDenseIndex DENSE and contiguous per rank. A writer emitting one shared
                    point list (VTK) needs every node written by exactly one
                    rank and every rank's connectivity indexing the same list,
                    so it cannot tolerate holes.

     nodeFinleyId   DERIVED, so that both ranks either side of a 2:1 seam
                    compute the same id for a shared hanging node without
                    communicating. Deriving it costs reserved slots, hence
                    holes - the exact thing nodeDenseIndex cannot have.

   Dense and derivable are mutually exclusive here; that is the whole reason
   there are three. They also differ in who OWNS a hanging node: the dense
   numbering gives it to a rank that has it as a cell corner (a fine-side rank),
   the finley numbering to the coarse octant's rank, whose simplices use it.

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
    /// global index of this rank's first element, in the forest's own element
    /// order. With it a local element names itself the same way on every rank,
    /// which is what lets the export give its simplices derivable ids.
    long globalElementOffset = 0;

    /// node coordinates, size numNodes*numDim (node i at [i*numDim + d])
    std::vector<double> nodeCoords;
    /// global id of each local node, size numNodes
    std::vector<long> nodeLnodesId;
    /// tag of each local node, size numNodes. A node materialised at a hanging
    /// position has no tag of its own - it is no node of the domain - so it
    /// inherits from its masters, see inheritedTag().
    std::vector<long> nodeTags;
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
    /// which face of that element this is, in p4est's own face order
    /// (-x,+x,-y,+y[,-z,+z]), size numFaces. Together with the element's global
    /// index it names a boundary face the same way on every rank, which is what
    /// lets the export give its face elements derivable ids.
    std::vector<long> faceDirections;

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
    // nodeCoords/nodeLnodesId. They carry NO sample of the domain's nodal
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
    // The DENSE numbering - for a writer emitting one shared point list.
    //
    // nodeLnodesId cannot serve: a rank's owned ids are not one contiguous block
    // once the materialised nodes are added, and those sit far above the rest.
    // Here rank r owns exactly [denseDistribution[r], denseDistribution[r+1]),
    // its lnodes-owned nodes first, then the hanging ones it writes.
    //
    // ORDER IS LOAD-BEARING: within a rank's range the index must INCREASE with
    // the local node index, because the writers walk local nodes, emit those
    // falling in their own range, and index them afterwards by this array
    // (weipa's OxleyNodes::writeCoordinatesVTK and DataVar::writeToVTK). Order
    // them any other way and points are written in one order and referenced in
    // another - the cells stop being cells.
    // ------------------------------------------------------------------------

    /// dense contiguous index of each local node, size numNodes
    std::vector<long> nodeDenseIndex;
    /// first dense index of each rank, size mpiSize+1; the last entry is the
    /// global number of nodes written
    std::vector<long> denseDistribution;

    // ------------------------------------------------------------------------
    // The FINLEY numbering: the ids handed over on export.
    //
    // Neither of the other two can do this job. lnodes does not name the hanging
    // positions at all, and the dense numbering is assigned in each rank's own
    // order, so the rank holding only the COARSE side of a 2:1 seam cannot work
    // out what the rank holding the fine side called the node. finley resolves
    // the mesh by global id, so the two must agree or the seam is a crack.
    //
    // This one is therefore DERIVED rather than assigned, from data p4est
    // already replicates everywhere. Rank r owns
    //     [ offset_r, offset_r + owned_r + 4*quads_r )
    // with owned_r its lnodes node count and quads_r its octant count, and
    //     lnodes node, lnodes global id g  ->  offset_r + (g - realOffset_r)
    //     hanging node, octant Q, face f   ->  offset_r + owned_r
    //                                           + 4*(Q - firstQuad_r) + f
    // A hanging position is the midpoint of exactly one coarse octant's face,
    // so (Q, f) names it uniquely and both sides of a seam form the same key.
    //
    // Contiguous per rank because finley's resolveNodeIds() allocates two dense
    // arrays spanning the id range of the local elements; the four reserved
    // slots per octant leave holes, which cost a constant factor on finley's
    // temporary labelling buffer and are packed away by createDenseDOFLabeling.
    // Being closed-form in both directions, it also decodes back to the oxley
    // node without a stored permutation.
    // ------------------------------------------------------------------------

    /// finley id of each local node, size numNodes
    std::vector<long> nodeFinleyId;
    /// first finley id of each rank, size mpiSize+1
    std::vector<long> finleyDistribution;
    /// Which rank writes each materialised node in the DENSE numbering, parallel
    /// to constrainedNodes. NOT the coarse-side owner the finley id is derived
    /// from - this is the lowest-numbered rank holding a FINE octant of the seam,
    /// i.e. one that actually has the node as a corner of one of its cells. The
    /// coarse side does not, since a hanging node is no corner of the coarse
    /// quad, so a coarse-side writer emits a point it can give no value to.
    std::vector<int> hangingWriterRank;
    /// local index of the hanging node on each element face, else -1; size
    /// numElements*2*numDim, face in p4est order. The coarse side of a seam has
    /// no other way to reach it - the node is a corner of the finer neighbour,
    /// so it is absent from this element's own corner list.
    std::vector<long> elementFaceHangingNode;
};

/// Materialised hanging nodes already created, keyed by position quantised to
/// 1e-9, so the fine elements meeting at one hanging position share a node.
typedef std::map<std::array<long,3>, long> HangingNodeMap;

/**
   \brief
   The tag a node materialised at a hanging position takes from its masters.

   The position is not a node of the oxley domain, so nothing tagged it and
   there is no prior answer to reuse - unlike a Dirac point, which addPoints()
   has already resolved. The rule: inherit when every master agrees, otherwise
   0 (untagged). Agreement means the whole edge the node sits on carries one
   tag, so a region tagged through its nodes keeps a consistent boundary; where
   the masters disagree the node is on the border between two tagged regions
   and picking either would be arbitrary.
*/
inline long inheritedTag(const MeshAccess& m, const long* masters, int n)
{
    long tag = 0;
    for (int k = 0; k < n; ++k) {
        const long mi = masters[k];
        if (mi < 0 || mi >= (long) m.nodeTags.size())
            return 0;
        if (k == 0)
            tag = m.nodeTags[mi];
        else if (m.nodeTags[mi] != tag)
            return 0;
    }
    return tag;
}

/**
   \brief
   Appends a node at a hanging position, or returns the one already there.

   \param m the mesh being built; numNodes, nodeCoords, nodeLnodesId and the
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
    m.nodeLnodesId.push_back(-1);
    if (!m.nodeTags.empty())
        m.nodeTags.push_back(inheritedTag(m, masters, n));
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
