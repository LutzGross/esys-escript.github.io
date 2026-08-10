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

#include <oxley/FinleyConverter.h>
#include <oxley/Brick.h>
#include <oxley/OxleyException.h>
#include <oxley/Rectangle.h>

#include <escript/FunctionSpaceFactory.h>
#include <finley/FinleyDomain.h>

#include <algorithm>
#include <limits>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <vector>

namespace oxley {

namespace {

/// The four corners of each octant face, in z-order corner indices, listed
/// counter-clockwise as seen from OUTSIDE the octant. Faces run -x,+x,-y,+y,-z,+z.
/// Same table getMeshAccess() winds its boundary faces with.
const int brickFaceCorners[6][4] = {
    {0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6} };

/// The boundary polygon of a 2D element, counter-clockwise, in z-order indices.
const int rectPolygon[4] = {0, 1, 3, 2};

/// is z-order corner c on face f of an octant?
inline bool cornerOnFace(int c, int f)
{
    switch (f) {
        case 0:  return !(c & 1);
        case 1:  return  (c & 1) != 0;
        case 2:  return !(c & 2);
        case 3:  return  (c & 2) != 0;
        case 4:  return !(c & 4);
        default: return  (c & 4) != 0;
    }
}

/// six times the signed volume of a tetrahedron, from local node indices
inline double signedVolume(const std::vector<double>& X, long a, long b,
                           long c, long d)
{
    const double u0 = X[b*3  ] - X[a*3  ], v0 = X[c*3  ] - X[a*3  ], w0 = X[d*3  ] - X[a*3  ];
    const double u1 = X[b*3+1] - X[a*3+1], v1 = X[c*3+1] - X[a*3+1], w1 = X[d*3+1] - X[a*3+1];
    const double u2 = X[b*3+2] - X[a*3+2], v2 = X[c*3+2] - X[a*3+2], w2 = X[d*3+2] - X[a*3+2];
    return u0*(v1*w2 - v2*w1) - u1*(v0*w2 - v2*w0) + u2*(v0*w1 - v1*w0);
}

/// twice the signed area of a triangle, from local node indices
inline double signedArea(const std::vector<double>& X, long a, long b, long c)
{
    return (X[b*2] - X[a*2]) * (X[c*2+1] - X[a*2+1])
         - (X[c*2] - X[a*2]) * (X[b*2+1] - X[a*2+1]);
}

// ---------------------------------------------------------------------------
// The 2D split, as a lookup rather than a decision.
//
// Work in the quad's counter-clockwise frame: corners V0..V3 (z-order indices
// rectPolygon), side s running Vs -> V(s+1), and Hs the hanging node on side s.
// The configuration is a 4-bit mask over sides, and rotating it collapses the
// 16 masks onto SIX patterns.
//
// Choosing the split by table rather than by a rule keeps neighbouring elements
// in the same configuration looking the same. The earlier version picked an apex
// by smallest global node id, which made the triangulation a function of the
// NUMBERING: two quads in identical situations came out split differently, which
// is what made the meshes look shredded. It is also the shape that generalises -
// in 3D a decision rule has to cope with faces and edges hanging in combination,
// where a table is just longer.
//
// Conformity does not depend on any of this in 2D: neighbours share an EDGE, an
// edge has no diagonal to disagree about, and which boundary edges the polygon
// has is fixed by the mask alone.
// ---------------------------------------------------------------------------

/// Ids reserved per octant for the simplices it is split into. The 2D split
/// emits at most 6 (the four-hanging pattern) and the conforming 3D cone 6, so
/// this covers both; it must grow if the 3D split ever handles hanging faces.
const int MAX_SIMPLICES = 6;

/// side of the CCW frame -> p4est face index
const int sideToFace[4] = {2, 1, 3, 0};

/// Symbolic vertices in a pattern: 0..3 are V0..V3, 4..7 are H0..H3.
struct SplitPattern
{
    int mask;               ///< the canonical configuration
    int numTriangles;
    int tri[6][3];
};

const SplitPattern splitPatterns[6] = {
    // no hanging node: one diagonal, always the same one in the octant's frame
    { 0x0, 2, { {0,1,2}, {0,2,3} } },
    // one: fan from the hanging node, which beats fanning past it
    { 0x1, 3, { {4,1,2}, {4,2,3}, {4,3,0} } },
    // two adjacent: cut the corner between them, fan the rest
    { 0x3, 4, { {4,1,5}, {4,5,2}, {4,2,3}, {4,3,0} } },
    // two opposite: halve the quad between the two hanging nodes, split each half
    { 0x5, 4, { {0,4,6}, {0,6,3}, {4,1,2}, {4,2,6} } },
    // three: cut both corners that lie between two hanging sides
    { 0x7, 5, { {4,1,5}, {5,2,6}, {5,6,3}, {5,3,0}, {5,0,4} } },
    // four: cut all four corners and halve what is left - symmetric, and no
    // interior node, which would need a fifth reserved id per octant
    { 0xf, 6, { {7,0,4}, {4,1,5}, {5,2,6}, {6,3,7}, {4,5,6}, {4,6,7} } }
};

/// Finds the rotation bringing `mask` onto one of the six canonical patterns.
/// Rotating by r means side s of the pattern is side (s+r)%4 of the element.
inline const SplitPattern* matchPattern(int mask, int& rotation)
{
    for (int r = 0; r < 4; ++r) {
        int cm = 0;
        for (int s = 0; s < 4; ++s)
            cm |= ((mask >> ((s + r) & 3)) & 1) << s;
        for (int p = 0; p < 6; ++p) {
            if (splitPatterns[p].mask == cm) {
                rotation = r;
                return &splitPatterns[p];
            }
        }
    }
    return NULL;                        // unreachable: the six cover all 16
}

/// Appends one triangle, orienting it by MEASURING the signed area rather than
/// reasoning about winding conventions.
inline void emitTriangle(std::vector<index_t>& elems, std::vector<int>& tags,
                         const std::vector<double>& X, long a, long b, long c,
                         int tag)
{
    if (signedArea(X, a, b, c) < 0.)
        std::swap(b, c);
    elems.push_back((index_t) a);
    elems.push_back((index_t) b);
    elems.push_back((index_t) c);
    tags.push_back(tag);
}

/// position within `n` of the entry with the smallest global node id
template <typename GetId>
inline int lowestIdPos(const int* corners, int n, GetId gid)
{
    int best = 0;
    for (int k = 1; k < n; ++k)
        if (gid(corners[k]) < gid(corners[best]))
            best = k;
    return best;
}

/**
   Verifies that the simplices tile the domain without cracks: every face of
   every simplex must be shared with exactly one other simplex, or lie on the
   boundary. A quad face split one way by one element and the other way by its
   neighbour produces faces with a count of one that are not on the boundary,
   which is exactly what this catches - and which no physics test would, since
   a globally linear function lies in both triangulations of a planar quad.

   In serial the check is exact. Under MPI a face on the rank interface also has
   a local count of one, so only the "never more than twice" half is checked;
   the split rule is a deterministic function of global node ids, so agreeing in
   serial implies agreeing across ranks.
*/
template <int NF, int FN>
void checkConformity(const std::vector<index_t>& elementNodes, int nodesPerElem,
                     const int (*localFaces)[FN],
                     const std::vector<index_t>& boundaryNodes,
                     bool exact, const char* what)
{
    typedef std::array<index_t, FN> Key;
    std::map<Key, int> count;
    const size_t ne = elementNodes.size() / nodesPerElem;
    for (size_t e = 0; e < ne; ++e) {
        for (int f = 0; f < NF; ++f) {
            Key k;
            for (int i = 0; i < FN; ++i)
                k[i] = elementNodes[e * nodesPerElem + localFaces[f][i]];
            std::sort(k.begin(), k.end());
            count[k]++;
        }
    }

    std::set<Key> boundary;
    const size_t nb = boundaryNodes.size() / FN;
    for (size_t b = 0; b < nb; ++b) {
        Key k;
        for (int i = 0; i < FN; ++i)
            k[i] = boundaryNodes[b * FN + i];
        std::sort(k.begin(), k.end());
        boundary.insert(k);
    }

    long tooMany = 0, danglingInterior = 0;
    for (typename std::map<Key, int>::const_iterator it = count.begin();
         it != count.end(); ++it) {
        if (it->second > 2)
            ++tooMany;
        else if (it->second == 1 && exact && boundary.find(it->first) == boundary.end())
            ++danglingInterior;
    }

    if (tooMany || danglingInterior) {
        std::stringstream ss;
        ss << "toFinley: the " << what << " mesh is not conforming - "
           << tooMany << " face(s) shared by more than two elements, "
           << danglingInterior << " interior face(s) with only one element. "
              "The face split rule disagreed between neighbours.";
        throw OxleyException(ss.str());
    }
}

const int tetFaces[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} };
const int triFaces[3][2] = { {0,1}, {1,2}, {0,2} };

} // anonymous namespace


namespace {

/**
   \brief
   Carries the domain's Dirac points across to the exported mesh.

   oxley has already done the hard part: addPoints() resolved every point to the
   nearest OWNED lnodes node - never a hanging position, which is no node of the
   oxley domain - and settled collectively which single rank keeps it. So there
   is nothing to locate here. The point belongs at that node's position, which
   is a node of the export too, the export's node set being a SUPERSET of the
   lnodes one (it materialises the hanging positions as well).

   Re-searching on the exported mesh would be worse than redundant: a point
   could then snap to a materialised hanging position, which oxley would never
   have chosen, and the two domains would disagree about where the source sits.

   The placement goes through finley's own addDiracPoints because
   createFromArrays() ends in prepare(), which redistributes nodes between
   ranks: any local node index computed before that is meaningless afterwards,
   and only the coordinates survive the move. finley then finds the node lying
   at distance zero, which is the one oxley picked.

   finley's addDiracPoints is COLLECTIVE and reduces over the whole point list,
   so every rank has to pass the SAME list, while oxley's is partitioned one
   point per owning rank. Hence the gather.
*/
void transferDiracPoints(const OxleyDomain& dom, const MeshAccess& m,
                         const escript::JMPI& mpiInfo,
                         escript::Domain_ptr& target)
{
    const std::vector<DiracPoint>& pts = dom.getDiracPoints();
    const int dim = m.numDim;

    std::vector<double> coords;
    std::vector<int> tags;
    coords.reserve(pts.size() * dim);
    tags.reserve(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        const long n = (long) pts[i].node;
        if (n < 0 || n >= m.numNodes)
            throw OxleyException("toFinley: a Dirac point names a node which "
                    "is not in the mesh view.");
        for (int d = 0; d < dim; ++d)
            coords.push_back(m.nodeCoords[(size_t) n * dim + d]);
        tags.push_back(pts[i].tag);
    }

#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        const int size = mpiInfo->size;
        int localCount = (int) tags.size();
        std::vector<int> perRank(size, 0);
        MPI_Allgather(&localCount, 1, MPI_INT, &perRank[0], 1, MPI_INT,
                      mpiInfo->comm);

        std::vector<int> tagDispl(size, 0), coordCount(size, 0),
                         coordDispl(size, 0);
        int total = 0;
        for (int r = 0; r < size; ++r) {
            tagDispl[r] = total;
            coordDispl[r] = total * dim;
            coordCount[r] = perRank[r] * dim;
            total += perRank[r];
        }

        std::vector<int> allTags(total, 0);
        std::vector<double> allCoords((size_t) total * dim, 0.);
        MPI_Allgatherv(tags.empty() ? NULL : &tags[0], localCount, MPI_INT,
                       allTags.empty() ? NULL : &allTags[0], &perRank[0],
                       &tagDispl[0], MPI_INT, mpiInfo->comm);
        MPI_Allgatherv(coords.empty() ? NULL : &coords[0], localCount * dim,
                       MPI_DOUBLE, allCoords.empty() ? NULL : &allCoords[0],
                       &coordCount[0], &coordDispl[0], MPI_DOUBLE,
                       mpiInfo->comm);
        tags.swap(allTags);
        coords.swap(allCoords);
    }
#endif

    // the list is now identical on every rank, so this test is too - and
    // addDiracPoints returns before its own reduction when there is nothing to
    // place, which would deadlock if only some ranks reached it.
    if (tags.empty())
        return;

    finley::FinleyDomain* fd = dynamic_cast<finley::FinleyDomain*>(target.get());
    if (!fd)
        throw OxleyException("toFinley: the exported mesh is not a finley "
                "domain, so its Dirac points cannot be set.");
    fd->addDiracPoints(coords, tags);
    // addDiracPoints fills the Point table but leaves the tags-in-use list
    // stale, so every caller inside finley follows it with this; without it
    // DiracDeltaFunctions has the points but reports no tags, and setting a
    // value by tag name finds nothing to set.
    fd->getPoints()->updateTagList();
}


/**
   For each local element, how many simplices the split emits and what fraction
   of the element each one covers.

   Recomputed rather than recorded: the split is a deterministic function of the
   hanging configuration and the element geometry, both of which MeshAccess
   already carries, so the weights can be worked out on the oxley side alone.
   That is what lets the inbound transfer weight by area without the finley side
   ever having to send areas along.

   2D only. The fractions sum to one per element, so a field that is constant
   over an element comes back unchanged whatever the split.
*/
void simplexWeights(const MeshAccess& m, std::vector<int>& count,
                    std::vector<std::vector<double> >& weight)
{
    if (m.numDim != 2)
        throw OxleyException("simplexWeights: 2D only so far.");
    const int V = m.nodesPerElement;
    count.assign(m.numElements, 0);
    weight.assign(m.numElements, std::vector<double>());

    for (long e = 0; e < m.numElements; ++e) {
        const long* en = &m.elementNodes[(size_t) e * V];
        long V4[4], H4[4];
        int mask = 0;
        for (int side = 0; side < 4; ++side) {
            V4[side] = en[rectPolygon[side]];
            H4[side] = m.elementFaceHangingNode.empty() ? -1
                     : m.elementFaceHangingNode[(size_t) e * 4 + sideToFace[side]];
            if (H4[side] >= 0)
                mask |= 1 << side;
        }
        int rot = 0;
        const SplitPattern* pat = matchPattern(mask, rot);
        count[e] = pat->numTriangles;

        double total = 0.;
        weight[e].resize(pat->numTriangles);
        for (int t = 0; t < pat->numTriangles; ++t) {
            long v[3];
            for (int k = 0; k < 3; ++k) {
                const int sym = pat->tri[t][k];
                v[k] = (sym < 4) ? V4[(sym + rot) & 3] : H4[(sym - 4 + rot) & 3];
            }
            const double a = std::fabs(signedArea(m.nodeCoords, v[0], v[1], v[2]));
            weight[e][t] = a;
            total += a;
        }
        if (total <= 0.)
            throw OxleyException("simplexWeights: an element has no area.");
        for (int t = 0; t < pat->numTriangles; ++t)
            weight[e][t] /= total;
    }
}

/// how an exported mesh says which forest it came from
std::string exportTag(const OxleyDomain& dom)
{
    std::stringstream ss;
    ss << "[forest " << dom.forestChecksum() << "]";     // collective
    return ss.str();
}

} // anonymous namespace

escript::Domain_ptr toFinley(const OxleyDomain& dom, int order,
                             int reducedOrder, bool optimize, bool simplices)
{
    const bool conforming = dom.isConforming();     // collective
    if (!conforming && (!simplices || dom.getDim() != 2))
        throw OxleyException("toFinley: the forest has hanging nodes. Only the "
                "2D simplex split handles them so far; the 3D split and the "
                "debug Rec4/Hex8 path still need a conforming forest.");

    // With hanging nodes present the mesh view must materialise them: a hanging
    // position is then a real node with a global id, so it can be a vertex of
    // the triangles on BOTH sides of a 2:1 seam. finley cannot represent a
    // hanging node (one ReferenceElementSet per ElementFile, and escript's q/r
    // is pointwise Dirichlet, not u = (u_a+u_b)/2), so the seam has to be
    // resolved by the triangulation instead - the node becomes an ordinary free
    // degree of freedom.
    const MeshAccess m = dom.getMeshAccess(!conforming);
    if (m.numDim != 2 && m.numDim != 3) {
        std::stringstream ss;
        ss << "toFinley: unsupported dimension " << m.numDim;
        throw OxleyException(ss.str());
    }

    const int V = m.nodesPerElement;
    const int FV = m.nodesPerFace;
    // The ids handed to finley. The export numbering is the one both sides of a
    // 2:1 seam derive independently, so a hanging node is one node across ranks;
    // nodeLnodesId cannot do that, since it numbers materialised nodes in each
    // rank's own creation order. Brick does not build it yet, hence the fallback.
    const std::vector<long>& gidOf =
            m.nodeFinleyId.empty() ? m.nodeLnodesId : m.nodeFinleyId;

    finley::MeshArrays out;
    out.numDim = m.numDim;

    // nodes: hand over every local node, owned or ghost. finley tolerates a
    // node being supplied by more than one rank and fetches any it still needs.
    out.nodeId.resize(m.numNodes);
    for (long i = 0; i < m.numNodes; ++i)
        out.nodeId[i] = (index_t) gidOf[i];
    out.nodeCoords = m.nodeCoords;

    // node tags. Unlike the element and face tags these do not ride along with
    // an entity of the mesh, so leaving them out lost them silently: a region
    // tagged through its nodes on the oxley side came out untagged here.
    if ((long) m.nodeTags.size() == m.numNodes) {
        out.nodeTag.resize(m.numNodes);
        for (long i = 0; i < m.numNodes; ++i)
            out.nodeTag[i] = (int) m.nodeTags[i];
    }

    // element node tables are built in LOCAL indices first, so the conformity
    // check and the orientation fix can use the coordinates, then translated
    // to global ids at the end.
    std::vector<index_t> elems, faces;

    if (!simplices) {
        // ---- debug path: one finley element per octant --------------------
        // Kept only for isolating problems in the finley handover from problems
        // in the split. NOT a production path: a forest would then change
        // element family the moment refinement introduces a hanging node.
        static const int zToFinley2D[4] = {0, 1, 3, 2};
        static const int zToFinley3D[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        const int* zToFinley = (m.numDim == 2) ? zToFinley2D : zToFinley3D;

        out.elementType = (m.numDim == 2) ? finley::Rec4 : finley::Hex8;
        out.faceElementType = (m.numDim == 2) ? finley::Line2 : finley::Rec4;

        elems.resize((size_t) m.numElements * V);
        out.elementTag.resize(m.numElements);
        for (long e = 0; e < m.numElements; ++e) {
            for (int c = 0; c < V; ++c)
                elems[(size_t) e * V + zToFinley[c]] =
                        (index_t) m.elementNodes[(size_t) e * V + c];
            out.elementTag[e] = (int) m.elementTags[e];
            out.elementId.push_back(
                    (index_t)((m.globalElementOffset + e) * MAX_SIMPLICES));
        }
        faces.resize((size_t) m.numFaces * FV);
        out.faceTag.resize(m.numFaces);
        for (long f = 0; f < m.numFaces; ++f) {
            for (int c = 0; c < FV; ++c)
                faces[(size_t) f * FV + c] =
                        (index_t) m.faceNodes[(size_t) f * FV + c];
            out.faceTag[f] = (int) m.faceTags[f];
        }

    } else if (m.numDim == 2) {
        // ---- 2D: each quad is split by the pattern its configuration names ---
        out.elementType = finley::Tri3;
        out.faceElementType = finley::Line2;

        if (!conforming && m.elementFaceHangingNode.empty())
            throw OxleyException("toFinley: the mesh view has no hanging node "
                    "table; getMeshAccess was not asked to materialise them.");

        elems.reserve((size_t) m.numElements * 3 * 3);
        out.elementTag.reserve(m.numElements * 3);
        for (long e = 0; e < m.numElements; ++e) {
            const long* en = &m.elementNodes[(size_t) e * V];

            // the element's own vertices, counter-clockwise, and the hanging
            // node on each side (-1 where the neighbour is not finer)
            long V4[4], H4[4];
            int mask = 0;
            for (int s = 0; s < 4; ++s) {
                V4[s] = en[rectPolygon[s]];
                H4[s] = m.elementFaceHangingNode.empty() ? -1
                      : m.elementFaceHangingNode[(size_t) e * 4 + sideToFace[s]];
                if (H4[s] >= 0)
                    mask |= 1 << s;
            }

            int rot = 0;
            const SplitPattern* pat = matchPattern(mask, rot);
            for (int t = 0; t < pat->numTriangles; ++t) {
                long v[3];
                for (int k = 0; k < 3; ++k) {
                    const int sym = pat->tri[t][k];
                    v[k] = (sym < 4) ? V4[(sym + rot) & 3]
                                     : H4[(sym - 4 + rot) & 3];
                }
                emitTriangle(elems, out.elementTag, m.nodeCoords,
                             v[0], v[1], v[2], (int) m.elementTags[e]);
                // an id that says which octant this triangle came from, so a
                // transfer can find its way back without a stored map. Sparse -
                // an octant reserves MAX_SIMPLICES ids and uses 2 to 6 of them -
                // which finley does not mind: element ids are only stored and
                // handed back, never used to size anything.
                out.elementId.push_back(
                        (index_t)((m.globalElementOffset + e) * MAX_SIMPLICES + t));
            }
        }
        // boundary edges are unchanged by the split: each is an edge of exactly
        // one of the fan triangles
        faces.resize((size_t) m.numFaces * FV);
        out.faceTag.resize(m.numFaces);
        for (long f = 0; f < m.numFaces; ++f) {
            for (int c = 0; c < FV; ++c)
                faces[(size_t) f * FV + c] =
                        (index_t) m.faceNodes[(size_t) f * FV + c];
            out.faceTag[f] = (int) m.faceTags[f];
        }

    } else {
        // ---- 3D: cone each octant from its lowest-id corner ----------------
        // over the three faces that do not contain it, each split by the
        // diagonal through ITS own lowest-id corner. Six tetrahedra. The three
        // faces that do contain the apex inherit a fan through it, which is the
        // same rule, because the apex is their lowest corner too.
        out.elementType = finley::Tet4;
        out.faceElementType = finley::Tri3;

        elems.reserve((size_t) m.numElements * 6 * 4);
        out.elementTag.reserve(m.numElements * 6);
        for (long e = 0; e < m.numElements; ++e) {
            const long* en = &m.elementNodes[(size_t) e * V];
            int apex = 0;
            for (int c = 1; c < 8; ++c)
                if (gidOf[en[c]] < gidOf[en[apex]])
                    apex = c;

            int tet = 0;                // numbers this octant's tets, for its ids
            for (int f = 0; f < 6; ++f) {
                if (cornerOnFace(apex, f))
                    continue;                    // near face, the cone covers it
                const int* q = brickFaceCorners[f];
                const int mpos = lowestIdPos(q, 4,
                        [&](int c) { return gidOf[en[c]]; });
                for (int t = 0; t < 2; ++t) {
                    long a = en[apex];
                    long b = en[q[mpos]];
                    long c = en[q[(mpos + 1 + t) % 4]];
                    long d = en[q[(mpos + 2 + t) % 4]];
                    if (signedVolume(m.nodeCoords, a, b, c, d) < 0.)
                        std::swap(c, d);
                    elems.push_back((index_t) a);
                    elems.push_back((index_t) b);
                    elems.push_back((index_t) c);
                    elems.push_back((index_t) d);
                    out.elementTag.push_back((int) m.elementTags[e]);
                    out.elementId.push_back((index_t)(
                            (m.globalElementOffset + e) * MAX_SIMPLICES
                            + (long) tet++));
                }
            }
        }

        // boundary quads split by the same rule, so they match the tet faces.
        // getMeshAccess() already wound them counter-clockwise seen from
        // outside, and fanning preserves that, so the normals stay outward.
        faces.reserve((size_t) m.numFaces * 2 * 3);
        out.faceTag.reserve(m.numFaces * 2);
        for (long f = 0; f < m.numFaces; ++f) {
            const long* fn = &m.faceNodes[(size_t) f * FV];
            static const int ident[4] = {0, 1, 2, 3};
            const int mpos = lowestIdPos(ident, 4,
                    [&](int c) { return gidOf[fn[c]]; });
            for (int t = 0; t < 2; ++t) {
                faces.push_back((index_t) fn[mpos]);
                faces.push_back((index_t) fn[(mpos + 1 + t) % 4]);
                faces.push_back((index_t) fn[(mpos + 2 + t) % 4]);
                out.faceTag.push_back((int) m.faceTags[f]);
            }
        }
    }

    // the split must leave no cracks; see checkConformity()
    const bool exact = (dom.getMPI()->size == 1);
    if (simplices) {
        if (m.numDim == 3)
            checkConformity<4, 3>(elems, 4, tetFaces, faces, exact, "tetrahedral");
        else
            checkConformity<3, 2>(elems, 3, triFaces, faces, exact, "triangular");
    }

    // local node indices -> global node ids
    out.elementNodes.resize(elems.size());
    for (size_t i = 0; i < elems.size(); ++i)
        out.elementNodes[i] = (index_t) gidOf[elems[i]];
    out.faceNodes.resize(faces.size());
    for (size_t i = 0; i < faces.size(); ++i)
        out.faceNodes[i] = (index_t) gidOf[faces[i]];

    // Element ids are left for finley to assign, which numbers them
    // consecutively across ranks. Face ids follow on from the elements, as
    // finley's own generators arrange them, so the two do not overlap.
    escript::JMPI mpiInfo = dom.getMPI();
    const long numElements = (long) (out.elementNodes.size() /
            finley::ReferenceElement::getInfo(out.elementType)->numNodes);
    const long numFaces = (long) (out.faceNodes.size() /
            finley::ReferenceElement::getInfo(out.faceElementType)->numNodes);

    index_t globalNumElements = (index_t) numElements;
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        index_t local = (index_t) numElements;
        MPI_Allreduce(&local, &globalNumElements, 1, MPI_DIM_T, MPI_SUM,
                      mpiInfo->comm);
    }
#endif
    // clear of the element ids, which now reserve MAX_SIMPLICES per octant
    index_t faceIdOffset = (index_t) globalNumElements * MAX_SIMPLICES;
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        index_t local = (index_t) numFaces;
        index_t scan = 0;
        MPI_Exscan(&local, &scan, 1, MPI_DIM_T, MPI_SUM, mpiInfo->comm);
        if (mpiInfo->rank == 0)
            scan = 0;
        faceIdOffset += scan;
    }
#endif
    out.faceId.resize(numFaces);
    for (long f = 0; f < numFaces; ++f)
        out.faceId[f] = faceIdOffset + (index_t) f;

    // Tag NAMES. The values travel with the elements and faces, but a name is
    // domain state, so it has to be copied or it is lost - and a user's own
    // name (setTagMap) or a Dirac tag would be silently unknown on the export.
    // Copy the source domain's map first: writing the boundary names literally
    // here, as this used to, does not merely lose names, it can CONTRADICT the
    // source if the user has remapped one of them. The defaults below are a
    // fallback for the names the source domain happens not to define.
    out.tagMap = dom.getTagMap();
    const std::pair<const char*, int> defaultTags[] = {
        { "left", 1 }, { "right", 2 }, { "bottom", 10 }, { "top", 20 },
        { "front", 100 }, { "back", 200 }
    };
    const int numDefaults = (m.numDim == 3) ? 6 : 4;
    for (int i = 0; i < numDefaults; ++i) {
        if (out.tagMap.find(defaultTags[i].first) == out.tagMap.end())
            out.tagMap[defaultTags[i].first] = defaultTags[i].second;
    }

    // The name carries a fingerprint of the forest so that a later transfer can
    // tell this mesh came from THIS forest. Without it two unrelated meshes
    // with overlapping id ranges would exchange values without complaint.
    std::stringstream name;
    name << "finley mesh from oxley " << (m.numDim == 2 ? "Rectangle" : "Brick")
         << " " << exportTag(dom);

    escript::Domain_ptr result = finley::FinleyDomain::createFromArrays(
            out, name.str(), order, reducedOrder, optimize, mpiInfo);

    transferDiracPoints(dom, m, mpiInfo, result);
    return result;
}

namespace {

/// the mesh view toFinley() would build for this domain, so the node ids and
/// the constraint arrays match the export exactly
MeshAccess viewMatchingExport(const OxleyDomain& dom)
{
    return dom.getMeshAccess(!dom.isConforming());      // collective
}


/// refuses a target that is not the export of this forest
void checkSameForest(const OxleyDomain& dom, const escript::AbstractDomain& other,
                     const char* what)
{
    // getDescription() answers "FinleyMesh" for every finley mesh, so the name
    // is what distinguishes them - toFinley() stamps the forest into it.
    const finley::FinleyDomain* fd =
            dynamic_cast<const finley::FinleyDomain*>(&other);
    const std::string tag = exportTag(dom);              // collective
    if (fd == NULL || fd->getName().find(tag) == std::string::npos)
        throw OxleyException(std::string(what) + ": that domain was not "
                "exported from this forest. Pass the domain toFinley() "
                "returned for it.");
}

/// global id -> local index over a domain's nodes
std::map<long,long> nodeIndexById(const escript::AbstractDomain& dom,
                                  int fsType, long n)
{
    const dim_t* ids = dom.borrowSampleReferenceIDs(fsType);
    std::map<long,long> byId;
    for (long i = 0; i < n; ++i)
        byId[(long) ids[i]] = i;
    return byId;
}

/// the ids the export gave the nodes: derived where Rectangle builds them,
/// the plain lnodes ids where it does not (Brick, conforming only)
const std::vector<long>& exportIds(const MeshAccess& m)
{
    return m.nodeFinleyId.empty() ? m.nodeLnodesId : m.nodeFinleyId;
}

void checkNodalSpace(const escript::Data& d, const char* what)
{
    const int fs = d.getFunctionSpace().getTypeCode();
    if (fs != Nodes && fs != DegreesOfFreedom && fs != ReducedNodes
            && fs != ReducedDegreesOfFreedom)
        throw OxleyException(std::string(what) + ": the data must live on "
                "ContinuousFunction or Solution.");
}


/**
   Moves values keyed by GLOBAL ID between two meshes that do not share a
   partition.

   This is finley's own idiom, from NodeFile::gather_global: cut the global id
   range evenly (JMPI::setDistribution), then pass ONE buffer around a ring,
   each rank writing the ids it can supply as the buffer passes over its block,
   and on a second lap reading out the ids it wants. The point is that neither
   side ever learns the other's partition - the ring visits every block, so all
   that has to agree is the id space itself, which is what the export numbering
   guarantees. That is exactly the difficulty here: finley redistributes its
   nodes in prepare(), so rank r's finley nodes are not rank r's oxley nodes.

   The buffer is the global id RANGE divided by the rank count, and the export
   numbering is sparse - it reserves four slots per octant for the positions
   that may hang - so the buffer is a few times larger than the node count. The
   cost stops there; nothing downstream sees the holes.

   \param mpi the communicator both meshes were built on
   \param numComp components per value
   \param haveId ids this rank can supply, haveVal their values
   \param wantId ids this rank needs, wantVal filled with their values
*/
void exchangeByGlobalId(const escript::JMPI& mpi, int numComp,
                        const std::vector<long>& haveId,
                        const std::vector<double>& haveVal,
                        const std::vector<long>& wantId,
                        std::vector<double>& wantVal)
{
    wantVal.assign(wantId.size() * numComp, 0.);
    std::vector<char> got(wantId.size(), 0);

    // the id range, over every rank: a rank may want an id it cannot supply
    long lo = std::numeric_limits<long>::max(), hi = std::numeric_limits<long>::min();
    for (size_t k = 0; k < haveId.size(); ++k) {
        lo = std::min(lo, haveId[k]); hi = std::max(hi, haveId[k]);
    }
    for (size_t k = 0; k < wantId.size(); ++k) {
        lo = std::min(lo, wantId[k]); hi = std::max(hi, wantId[k]);
    }
    if (haveId.empty() && wantId.empty()) { lo = 0; hi = 0; }
#ifdef ESYS_MPI
    if (mpi->size > 1) {
        long gl = lo, gh = hi;
        MPI_Allreduce(&gl, &lo, 1, MPI_LONG, MPI_MIN, mpi->comm);
        MPI_Allreduce(&gh, &hi, 1, MPI_LONG, MPI_MAX, mpi->comm);
    }
#endif

    std::vector<index_t> dist(mpi->size + 1, 0);
    const dim_t bufLen = mpi->setDistribution((index_t) lo, (index_t) hi, &dist[0]);
    const long UNSET = lo - 1;

    std::vector<long> idBuf((size_t) bufLen, UNSET);
    std::vector<double> valBuf((size_t) bufLen * numComp, 0.);

    int bufferRank = mpi->rank;
#ifdef ESYS_MPI
    const int dest = mpi->mod_rank(mpi->rank + 1);
    const int source = mpi->mod_rank(mpi->rank - 1);
    MPI_Status status;
#endif

    // lap one: fill. Every rank writes what it can supply for whichever block
    // the buffer is currently carrying; after size steps the buffer is back
    // with its owner, holding every value anyone had for that block.
    for (int p = 0; p < mpi->size; ++p) {
#ifdef ESYS_MPI
        if (p > 0) {
            MPI_Sendrecv_replace(&idBuf[0], (int) bufLen, MPI_LONG, dest,
                    mpi->counter(), source, mpi->counter(), mpi->comm, &status);
            MPI_Sendrecv_replace(&valBuf[0], (int) (bufLen * numComp),
                    MPI_DOUBLE, dest, mpi->counter()+1, source,
                    mpi->counter()+1, mpi->comm, &status);
            mpi->incCounter(2);
        }
#endif
        bufferRank = mpi->mod_rank(bufferRank - 1);
        const long first = (long) dist[bufferRank], last = (long) dist[bufferRank+1];
        for (size_t k = 0; k < haveId.size(); ++k) {
            const long id = haveId[k];
            if (id < first || id >= last)
                continue;
            const size_t pos = (size_t)(id - first);
            idBuf[pos] = id;
            for (int c = 0; c < numComp; ++c)
                valBuf[pos*numComp + c] = haveVal[k*numComp + c];
        }
    }

    // lap two: read. The buffer starts on its owner and visits every rank, so
    // each rank sees every block exactly once.
    bufferRank = mpi->rank;
    for (int p = 0; p < mpi->size; ++p) {
        const long first = (long) dist[bufferRank], last = (long) dist[bufferRank+1];
        for (size_t k = 0; k < wantId.size(); ++k) {
            const long id = wantId[k];
            if (id < first || id >= last)
                continue;
            const size_t pos = (size_t)(id - first);
            if (idBuf[pos] == UNSET)
                continue;               // nobody supplied it; reported below
            for (int c = 0; c < numComp; ++c)
                wantVal[k*numComp + c] = valBuf[pos*numComp + c];
            got[k] = 1;
        }
#ifdef ESYS_MPI
        if (p < mpi->size - 1) {
            MPI_Sendrecv_replace(&idBuf[0], (int) bufLen, MPI_LONG, dest,
                    mpi->counter(), source, mpi->counter(), mpi->comm, &status);
            MPI_Sendrecv_replace(&valBuf[0], (int) (bufLen * numComp),
                    MPI_DOUBLE, dest, mpi->counter()+1, source,
                    mpi->counter()+1, mpi->comm, &status);
            mpi->incCounter(2);
        }
#endif
        bufferRank = mpi->mod_rank(bufferRank - 1);
    }

    for (size_t k = 0; k < wantId.size(); ++k) {
        if (!got[k]) {
            std::stringstream ss;
            ss << "transfer between the forest and its export: no value was "
                  "supplied for node id " << wantId[k] << ". Is the target the "
                  "domain toFinley() built from this forest?";
            throw OxleyException(ss.str());
        }
    }
}

} // anonymous namespace

escript::Data toFinleyData(const escript::Data& source, escript::Domain_ptr target)
{
    checkNodalSpace(source, "toFinleyData");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(
            source.getFunctionSpace().getDomain().get());
    if (dom == NULL)
        throw OxleyException("toFinleyData: the source must live on an oxley "
                "domain.");
    if (target.get() == NULL)
        throw OxleyException("toFinleyData: no target domain given.");
    if (source.isComplex())
        throw OxleyException("toFinleyData: complex data is not supported yet.");
    checkSameForest(*dom, *target, "toFinleyData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const std::vector<long>& gid = exportIds(m);
    const int numComp = source.getDataPointSize();

    // what this rank can supply: its own nodes, plus the seam positions it
    // materialised, which are no nodes of the forest and so take the average
    // of their masters
    std::vector<long> haveId;
    std::vector<double> haveVal;
    haveId.reserve(m.numNodes);
    haveVal.reserve((size_t) m.numNodes * numComp);
    for (long i = 0; i < m.numRealNodes; ++i) {
        const double* in = source.getSampleDataRO(i, (double) 0);
        haveId.push_back(gid[i]);
        for (int c = 0; c < numComp; ++c)
            haveVal.push_back(in[c]);
    }
    const int mpc = m.mastersPerConstrainedNode;
    for (size_t k = 0; k < m.constrainedNodes.size(); ++k) {
        const long node = m.constrainedNodes[k];
        std::vector<double> v(numComp, 0.);
        for (int j = 0; j < mpc; ++j) {
            const long master = m.constraintMasters[k*mpc + j];
            const double w = m.constraintWeights[k*mpc + j];
            if (master < 0 || w == 0.)
                continue;
            const double* in = source.getSampleDataRO(master, (double) 0);
            for (int c = 0; c < numComp; ++c)
                v[c] += w * in[c];
        }
        haveId.push_back(gid[node]);
        for (int c = 0; c < numComp; ++c)
            haveVal.push_back(v[c]);
    }

    escript::Data result(0., source.getDataPointShape(),
                         escript::continuousFunction(*target), true);
    result.requireWrite();
    const escript::FunctionSpace targetFS = escript::continuousFunction(*target);
    const long n = (long) result.getNumSamples();
    const dim_t* ids = target->borrowSampleReferenceIDs(targetFS.getTypeCode());
    std::vector<long> wantId(n);
    for (long j = 0; j < n; ++j)
        wantId[j] = (long) ids[j];

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numComp, haveId, haveVal, wantId, wantVal);

    for (long j = 0; j < n; ++j) {
        double* out = result.getSampleDataRW(j, (double) 0);
        for (int c = 0; c < numComp; ++c)
            out[c] = wantVal[(size_t) j*numComp + c];
    }
    return result;
}

escript::Data fromFinleyData(const escript::Data& source, escript::Domain_ptr target)
{
    checkNodalSpace(source, "fromFinleyData");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(target.get());
    if (dom == NULL)
        throw OxleyException("fromFinleyData: the target must be an oxley "
                "domain.");
    if (source.isComplex())
        throw OxleyException("fromFinleyData: complex data is not supported "
                "yet.");
    checkSameForest(*dom, *(source.getFunctionSpace().getDomain()),
                    "fromFinleyData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const std::vector<long>& gid = exportIds(m);
    const int numComp = source.getDataPointSize();

    const escript::FunctionSpace sourceFS = source.getFunctionSpace();
    const long ns = (long) source.getNumSamples();
    const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
            sourceFS.getTypeCode());
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * numComp);
    for (long j = 0; j < ns; ++j) {
        haveId[j] = (long) ids[j];
        const double* in = source.getSampleDataRO(j, (double) 0);
        for (int c = 0; c < numComp; ++c)
            haveVal[(size_t) j*numComp + c] = in[c];
    }

    // only the shared nodes come back; a materialised seam position is no node
    // of the forest, so its value has nowhere to go
    std::vector<long> wantId(m.numRealNodes);
    for (long i = 0; i < m.numRealNodes; ++i)
        wantId[i] = gid[i];

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numComp, haveId, haveVal, wantId, wantVal);

    escript::Data result(0., source.getDataPointShape(),
                         escript::continuousFunction(*target), true);
    result.requireWrite();
    for (long i = 0; i < m.numRealNodes; ++i) {
        double* out = result.getSampleDataRW(i, (double) 0);
        for (int c = 0; c < numComp; ++c)
            out[c] = wantVal[(size_t) i*numComp + c];
    }
    return result;
}


namespace {

/// accumulates the global id <-> position pairing seen at octant corner slots
/// and counts disagreements in BOTH directions. Either one breaks a mesh built
/// from these ids: one id at two positions tears the two apart, one position
/// under two ids leaves a crack between the elements that disagree.
struct CornerIdCheck
{
    typedef std::array<double,3> Pos;
    typedef std::array<long,3> Key;          // position quantised to 1e-9

    std::map<long, Pos> posOf;               // global id -> position
    std::map<Key, long> idAt;                // position -> global id
    long slots = 0;                          // corner slots examined
    long clashes = 0;                        // one id, several positions
    long splits = 0;                         // one position, several ids

    static Key key(const Pos& p)
    {
        Key k;
        for (int d = 0; d < 3; ++d)
            k[d] = (long) std::llround(p[d] * 1e9);
        return k;
    }

    void add(long gid, double x, double y, double z)
    {
        ++slots;
        const Pos p = {{x, y, z}};

        std::map<long, Pos>::iterator it = posOf.find(gid);
        if (it == posOf.end()) {
            posOf[gid] = p;
        } else {
            for (int d = 0; d < 3; ++d) {
                if (std::abs(it->second[d] - p[d]) > 1e-9) { ++clashes; break; }
            }
        }

        std::map<Key, long>::iterator jt = idAt.find(key(p));
        if (jt == idAt.end())
            idAt[key(p)] = gid;
        else if (jt->second != gid)
            ++splits;
    }
};

} // anonymous namespace

/// One traversal of a degree-2 lnodes: emits every element_nodes slot with the
/// global id it holds and the geometric position of the slot itself. Both the
/// report and the raw dump read this; the point of the exercise is exactly the
/// comparison between the two, so they must come from the same walk.
static void collectDegree2Slots(const OxleyDomain& dom,
                                std::vector<SlotRecord>& out,
                                long& numLocalNodes)
{
    const Rectangle* rect = dynamic_cast<const Rectangle*>(&dom);
    const Brick* brick = dynamic_cast<const Brick*>(&dom);
    out.clear();
    numLocalNodes = 0;

    if (rect) {
        p4est_t* p4 = rect->p4est;
        p4est_ghost_t* gh = p4est_ghost_new(p4, P4EST_CONNECT_FULL);
        p4est_lnodes_t* ln = p4est_lnodes_new(p4, gh, 2);
        numLocalNodes = ln->num_local_nodes;
        const int vn = ln->vnodes;               // 9 for degree 2 in 2D
        out.reserve((size_t) ln->num_local_elements * vn);
        long e = 0;
        for (p4est_topidx_t t = p4->first_local_tree; t <= p4->last_local_tree; ++t) {
            p4est_tree_t* tree = p4est_tree_array_index(p4->trees, t);
            sc_array_t* quads = &tree->quadrants;
            const p4est_locidx_t Q = (p4est_locidx_t) quads->elem_count;
            for (p4est_locidx_t q = 0; q < Q; ++q, ++e) {
                p4est_quadrant_t* quad = p4est_quadrant_array_index(quads, q);
                const p4est_qcoord_t half = P4EST_QUADRANT_LEN(quad->level) / 2;
                for (int slot = 0; slot < vn; ++slot) {
                    // lexicographic yx order: slot = i + 3*j, i,j in {0,1,2}
                    const int i = slot % 3, j = (slot / 3) % 3;
                    const p4est_locidx_t lid =
                            ln->element_nodes[(size_t) e * vn + slot];
                    SlotRecord r;
                    r.element = e;
                    r.slot = slot;
                    r.gid = (long) p4est_lnodes_global_index(ln, lid);
                    r.faceCode = (int) ln->face_code[e];
                    r.corner = (i != 1 && j != 1);
                    double xy[3] = {0., 0., 0.};
                    p4est_qcoord_to_vertex(p4->connectivity, t,
                            quad->x + i * half, quad->y + j * half, xy);
                    r.x = xy[0]; r.y = xy[1]; r.z = 0.;
                    out.push_back(r);
                }
            }
        }
        p4est_lnodes_destroy(ln);
        p4est_ghost_destroy(gh);

    } else if (brick) {
        p8est_t* p8 = brick->p8est;
        p8est_ghost_t* gh = p8est_ghost_new(p8, P8EST_CONNECT_FULL);
        p8est_lnodes_t* ln = p8est_lnodes_new(p8, gh, 2);
        numLocalNodes = ln->num_local_nodes;
        const int vn = ln->vnodes;                // 27 for degree 2 in 3D
        out.reserve((size_t) ln->num_local_elements * vn);
        long e = 0;
        for (p4est_topidx_t t = p8->first_local_tree; t <= p8->last_local_tree; ++t) {
            p8est_tree_t* tree = p8est_tree_array_index(p8->trees, t);
            sc_array_t* octs = &tree->quadrants;
            const p4est_locidx_t Q = (p4est_locidx_t) octs->elem_count;
            for (p4est_locidx_t q = 0; q < Q; ++q, ++e) {
                p8est_quadrant_t* oct = p8est_quadrant_array_index(octs, q);
                const p4est_qcoord_t half = P8EST_QUADRANT_LEN(oct->level) / 2;
                for (int slot = 0; slot < vn; ++slot) {
                    // lexicographic zyx order: slot = i + 3*j + 9*k
                    const int i = slot % 3, j = (slot / 3) % 3, k = slot / 9;
                    const p4est_locidx_t lid =
                            ln->element_nodes[(size_t) e * vn + slot];
                    SlotRecord r;
                    r.element = e;
                    r.slot = slot;
                    r.gid = (long) p8est_lnodes_global_index(ln, lid);
                    r.faceCode = (int) ln->face_code[e];
                    r.corner = (i != 1 && j != 1 && k != 1);
                    double xyz[3] = {0., 0., 0.};
                    p8est_qcoord_to_vertex(p8->connectivity, t,
                            oct->x + i * half, oct->y + j * half, oct->z + k * half,
                            xyz);
                    r.x = xyz[0]; r.y = xyz[1]; r.z = xyz[2];
                    out.push_back(r);
                }
            }
        }
        p8est_lnodes_destroy(ln);
        p8est_ghost_destroy(gh);

    } else {
        throw OxleyException("degree-2 lnodes probe: domain is neither a "
                             "Rectangle nor a Brick");
    }
}

std::vector<long> lnodesDegree2Report(const OxleyDomain& dom)
{
    std::vector<long> out(8, 0);
    std::vector<SlotRecord> slots;
    long numLocalNodes = 0;
    collectDegree2Slots(dom, slots, numLocalNodes);

    CornerIdCheck chk;
    long numOctants = 0, hangingOctants = 0, lastElement = -1;
    for (size_t s = 0; s < slots.size(); ++s) {
        const SlotRecord& r = slots[s];
        if (r.element != lastElement) {
            lastElement = r.element;
            ++numOctants;
            if (r.faceCode != 0)
                ++hangingOctants;
        }
        if (r.corner)
            chk.add(r.gid, r.x, r.y, r.z);
    }

    out[0] = numOctants;
    out[1] = chk.slots;
    out[2] = (long) chk.posOf.size();
    out[3] = chk.clashes;
    out[4] = hangingOctants;
    out[5] = numLocalNodes;
    out[6] = (long) chk.idAt.size();
    out[7] = chk.splits;
    return out;
}

std::vector<SlotRecord> lnodesDegree2Slots(const OxleyDomain& dom)
{
    std::vector<SlotRecord> slots;
    long numLocalNodes = 0;
    collectDegree2Slots(dom, slots, numLocalNodes);
    return slots;
}


escript::Data toFinleyReducedData(const escript::Data& source,
                                  escript::Domain_ptr target)
{
    if (source.getFunctionSpace().getTypeCode() != ReducedElements)
        throw OxleyException("toFinleyReducedData: the data must live on "
                "ReducedFunction.");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(
            source.getFunctionSpace().getDomain().get());
    if (dom == NULL)
        throw OxleyException("toFinleyReducedData: the source must live on an "
                "oxley domain.");
    if (target.get() == NULL)
        throw OxleyException("toFinleyReducedData: no target domain given.");
    if (source.isComplex())
        throw OxleyException("toFinleyReducedData: complex data is not "
                "supported yet.");
    checkSameForest(*dom, *target, "toFinleyReducedData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = source.getDataPointSize();
    std::vector<int> count;
    std::vector<std::vector<double> > weight;
    simplexWeights(m, count, weight);

    // one octant value, repeated onto each simplex it was split into
    std::vector<long> haveId;
    std::vector<double> haveVal;
    for (long e = 0; e < m.numElements; ++e) {
        const double* in = source.getSampleDataRO(e, (double) 0);
        for (int t = 0; t < count[e]; ++t) {
            haveId.push_back((m.globalElementOffset + e) * MAX_SIMPLICES + t);
            for (int c = 0; c < numComp; ++c)
                haveVal.push_back(in[c]);
        }
    }

    escript::Data result(0., source.getDataPointShape(),
                         escript::reducedFunction(*target), true);
    result.requireWrite();
    const escript::FunctionSpace targetFS = escript::reducedFunction(*target);
    const long n = (long) result.getNumSamples();
    const dim_t* ids = target->borrowSampleReferenceIDs(targetFS.getTypeCode());
    std::vector<long> wantId(n);
    for (long j = 0; j < n; ++j)
        wantId[j] = (long) ids[j];

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numComp, haveId, haveVal, wantId, wantVal);

    for (long j = 0; j < n; ++j) {
        double* out = result.getSampleDataRW(j, (double) 0);
        for (int c = 0; c < numComp; ++c)
            out[c] = wantVal[(size_t) j*numComp + c];
    }
    return result;
}

escript::Data fromFinleyReducedData(const escript::Data& source,
                                    escript::Domain_ptr target)
{
    if (source.getFunctionSpace().getTypeCode() != ReducedElements)
        throw OxleyException("fromFinleyReducedData: the data must live on "
                "ReducedFunction.");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(target.get());
    if (dom == NULL)
        throw OxleyException("fromFinleyReducedData: the target must be an "
                "oxley domain.");
    if (source.isComplex())
        throw OxleyException("fromFinleyReducedData: complex data is not "
                "supported yet.");
    checkSameForest(*dom, *(source.getFunctionSpace().getDomain()),
                    "fromFinleyReducedData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = source.getDataPointSize();
    std::vector<int> count;
    std::vector<std::vector<double> > weight;
    simplexWeights(m, count, weight);

    const escript::FunctionSpace sourceFS = source.getFunctionSpace();
    const long ns = (long) source.getNumSamples();
    const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
            sourceFS.getTypeCode());
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * numComp);
    for (long j = 0; j < ns; ++j) {
        haveId[j] = (long) ids[j];
        const double* in = source.getSampleDataRO(j, (double) 0);
        for (int c = 0; c < numComp; ++c)
            haveVal[(size_t) j*numComp + c] = in[c];
    }

    // ask for every simplex of every octant this rank owns
    std::vector<long> wantId;
    for (long e = 0; e < m.numElements; ++e)
        for (int t = 0; t < count[e]; ++t)
            wantId.push_back((m.globalElementOffset + e) * MAX_SIMPLICES + t);

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numComp, haveId, haveVal, wantId, wantVal);

    // and average them by the area each covers, so a field that is constant
    // over the octant comes back unchanged whatever the split
    escript::Data result(0., source.getDataPointShape(),
                         escript::reducedFunction(*target), true);
    result.requireWrite();
    size_t pos = 0;
    for (long e = 0; e < m.numElements; ++e) {
        double* out = result.getSampleDataRW(e, (double) 0);
        for (int c = 0; c < numComp; ++c)
            out[c] = 0.;
        for (int t = 0; t < count[e]; ++t, ++pos)
            for (int c = 0; c < numComp; ++c)
                out[c] += weight[e][t] * wantVal[pos*numComp + c];
    }
    return result;
}

} // namespace oxley
