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
#include <array>
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

/**
   Ids reserved per element for the simplices it is split into, and per boundary
   face for the ones IT is split into.

   2D emits at most six triangles, the four-hanging pattern. A conforming octant
   is six tetrahedra, but a hanging one is coned from its centre over a boundary
   whose six faces carry up to eight triangles each - a hanging face is four
   sub-quads of two - so forty-eight. A 2D boundary edge is never split; a 3D
   boundary quad is two triangles, or up to six when all four of its edges carry
   a hanging midpoint and it has become an octagon.

   The reservation is what makes the ids DERIVABLE: simplex t of global element
   Q is Q*maxSimplices + t on every rank, with no map to store and nothing to
   agree. The unused ids are holes, which finley does not mind - it only stores
   element ids and hands them back.
*/
inline int maxSimplices(int numDim) { return (numDim == 2) ? 6 : 48; }
inline int maxFaceSimplices(int numDim) { return (numDim == 2) ? 2 : 6; }

/// z-order corner pairs of the twelve octant edges, p8est's own numbering
const int brickEdgeCorners[12][2] = {
    {0,1}, {2,3}, {4,5}, {6,7}, {0,2}, {1,3},
    {4,6}, {5,7}, {0,4}, {1,5}, {2,6}, {3,7} };

/// which edge of an octant joins two of its z-order corners, -1 if they are not
/// the ends of one
inline int edgeOfCorners(int a, int b)
{
    for (int ed = 0; ed < 12; ++ed) {
        if ((brickEdgeCorners[ed][0] == a && brickEdgeCorners[ed][1] == b)
         || (brickEdgeCorners[ed][0] == b && brickEdgeCorners[ed][1] == a))
            return ed;
    }
    return -1;
}

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

/**
   Verifies that the simplices tile the domain without cracks: every face of
   every simplex must be shared with exactly one other simplex, or lie on the
   boundary. A quad face split one way by one element and the other way by its
   neighbour produces faces used once that are not on the boundary, which is
   exactly what this catches - and which no physics test would, since a globally
   linear function lies in both triangulations of a planar quad.

   Faces are keyed by their SORTED GLOBAL node ids, so a key means the same face
   on every rank. That is what lets the rank interfaces be checked too: a face
   there is used once locally and once on the neighbour, which looks exactly
   like a crack from either side alone. The ones that cannot be settled locally -
   used once, and not on the domain boundary - are routed by a hash of the key
   to whichever rank is to adjudicate, and each must arrive exactly twice.

   Only those are exchanged, not every face, so the traffic is the size of the
   interface rather than of the mesh. Nothing is lost by that: a face used twice
   on one rank and once on another is not sent by the rank that has it twice, so
   the one copy arrives alone and is reported - as a dangling face rather than
   as an over-used one, but reported.

   Collective, and the verdict is a reduction, so every rank throws or none does.
*/
template <int NF, int FN>
void checkConformity(const std::vector<index_t>& elementNodes, int nodesPerElem,
                     const int (*localFaces)[FN],
                     const std::vector<index_t>& boundaryNodes,
                     const std::vector<long>& gidOf,
                     const escript::JMPI& mpi, const char* what)
{
    typedef std::array<index_t, FN> Key;
    std::map<Key, int> count;
    const size_t ne = elementNodes.size() / nodesPerElem;
    for (size_t e = 0; e < ne; ++e) {
        for (int f = 0; f < NF; ++f) {
            Key k;
            for (int i = 0; i < FN; ++i)
                k[i] = (index_t) gidOf[elementNodes[e * nodesPerElem
                                                  + localFaces[f][i]]];
            std::sort(k.begin(), k.end());
            count[k]++;
        }
    }

    std::set<Key> boundary;
    const size_t nb = boundaryNodes.size() / FN;
    for (size_t b = 0; b < nb; ++b) {
        Key k;
        for (int i = 0; i < FN; ++i)
            k[i] = (index_t) gidOf[boundaryNodes[b * FN + i]];
        std::sort(k.begin(), k.end());
        boundary.insert(k);
    }

    long tooMany = 0;
    std::vector<Key> undecided;             // used once, and not on the boundary
    for (typename std::map<Key, int>::const_iterator it = count.begin();
         it != count.end(); ++it) {
        if (it->second > 2)
            ++tooMany;
        else if (it->second == 1 && boundary.find(it->first) == boundary.end())
            undecided.push_back(it->first);
    }

    long dangling = 0;
    if (mpi->size == 1) {
        dangling = (long) undecided.size();     // nobody else can be holding it
    }
#ifdef ESYS_MPI
    else {
        // route each key to the rank that will adjudicate it. Any function of
        // the key will do as long as every rank computes the same one, which is
        // why the key is global ids and sorted.
        const int size = mpi->size;
        std::vector<int> sendCount(size, 0), recvCount(size, 0);
        std::vector<int> owner(undecided.size());
        for (size_t i = 0; i < undecided.size(); ++i) {
            unsigned long h = 1469598103934665603UL;
            for (int j = 0; j < FN; ++j) {
                h ^= (unsigned long) undecided[i][j];
                h *= 1099511628211UL;
            }
            owner[i] = (int) (h % (unsigned long) size);
            sendCount[owner[i]] += FN;
        }
        MPI_Alltoall(&sendCount[0], 1, MPI_INT, &recvCount[0], 1, MPI_INT,
                     mpi->comm);

        std::vector<int> sendDispl(size, 0), recvDispl(size, 0);
        long sendTotal = 0, recvTotal = 0;
        for (int r = 0; r < size; ++r) {
            sendDispl[r] = (int) sendTotal;
            recvDispl[r] = (int) recvTotal;
            sendTotal += sendCount[r];
            recvTotal += recvCount[r];
        }
        std::vector<index_t> sendBuf(sendTotal), recvBuf(recvTotal);
        std::vector<int> fill(sendDispl);
        for (size_t i = 0; i < undecided.size(); ++i) {
            for (int j = 0; j < FN; ++j)
                sendBuf[fill[owner[i]]++] = undecided[i][j];
        }
        MPI_Alltoallv(sendBuf.empty() ? NULL : &sendBuf[0], &sendCount[0],
                      &sendDispl[0], MPI_DIM_T,
                      recvBuf.empty() ? NULL : &recvBuf[0], &recvCount[0],
                      &recvDispl[0], MPI_DIM_T, mpi->comm);

        std::map<Key, int> arrivals;
        for (long i = 0; i + FN <= recvTotal; i += FN) {
            Key k;
            for (int j = 0; j < FN; ++j)
                k[j] = recvBuf[i + j];
            arrivals[k]++;                  // already sorted by the sender
        }
        for (typename std::map<Key, int>::const_iterator it = arrivals.begin();
             it != arrivals.end(); ++it) {
            if (it->second != 2)
                ++dangling;
        }
    }

    if (mpi->size > 1) {
        long local[2] = { tooMany, dangling }, total[2] = { 0, 0 };
        MPI_Allreduce(local, total, 2, MPI_LONG, MPI_SUM, mpi->comm);
        tooMany = total[0];
        dangling = total[1];
    }
#endif

    if (tooMany || dangling) {
        std::stringstream ss;
        ss << "toFinley: the " << what << " mesh is not conforming - "
           << tooMany << " face(s) shared by more than two elements, "
           << dangling << " interior face(s) with only one element. "
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


// ---------------------------------------------------------------------------
// 3D transfers.
//
// 2D can rearrange: a boundary face corresponds to exactly one face of the
// export, and only Function has to evaluate anything. In 3D neither holds - one
// octant becomes six tetrahedra and one boundary quad becomes two triangles -
// so every space but ContinuousFunction has to EVALUATE, and the two sides'
// quadrature rules have nothing to do with each other.
//
// The rule used throughout, in both directions: a sample's values travel
// together with the COORDINATES of the points they were taken at, and the
// receiver rebuilds the polynomial those points are unisolvent for. Neither
// side then needs to know the other's quadrature rule or point order. The 2D
// code does know finley's, having measured it once (triangleQuadPoints); this
// costs a dozen doubles per sample and cannot go quietly wrong if a rule ever
// changes.
//
//   oxley -> finley: the octant's 2x2x2 points are a tensor grid, unisolvent
//       for a TRILINEAR function - bilinear on a boundary face, whose normal
//       axis simply has no spread - so two-point Lagrange per axis is exact.
//   finley -> oxley: a simplex's quadrature points are unisolvent for an
//       AFFINE function, so the affine weights through those points are exact.
//       Which simplex to read is decided by the octant point's position
//       against the simplex's VERTICES, since the quadrature points span only
//       part of it.
// ---------------------------------------------------------------------------

/**
   The split of one octant, and of one of its faces, recomputed rather than
   stored - the same spirit as splitOfElement in 2D. Both routines produce
   exactly what toFinley() emits, and in the same order, so simplex t here
   carries the id (global octant)*maxSimplices + t there.

   THE RULE, in one line: every face of an octant is fanned from its LOWEST-ID
   vertex, and a hanging face is not one polygon but the four faces of the finer
   octants across it, each fanned by that same rule.

   Why the rule has to be a property of the FACE and not of the octant looking at
   it: two octants share a face, and each has to cut it into the same triangles
   or the mesh is cracked along it. The lowest-id vertex is the same vertex from
   either side however the two octants number their own corners. For a plain quad
   this is the diagonal through the lowest-id corner, which is the rule the
   conforming split already used - a hanging face just has more vertices.

   And why the hanging face is done as four sub-quads rather than as one octagon:
   the fine side has already cut it, each of its four octants seeing a plain quad
   and fanning it. Fanning the octagon instead would be a perfectly good
   triangulation of the same eight vertices which simply does not MATCH, and a
   mismatched diagonal is a crack that no physics test can see - a linear field
   lies in both triangulations. Only the conformity check finds it.
*/

/// the node at the midpoint of the octant edge joining two of its corners, or -1
inline long edgeNode(const MeshAccess& m, long e, int a, int b)
{
    if (m.elementEdgeHangingNode.empty())
        return -1;
    const int ed = edgeOfCorners(a, b);
    return (ed < 0) ? -1 : m.elementEdgeHangingNode[(size_t) e * 12 + ed];
}

/**
   Fans a polygon, keeping the winding it came in with. Returns n-2 triangles.

   From the lowest-id vertex, except that a vertex sitting between two COLLINEAR
   neighbours - an edge midpoint, where the polygon runs corner, midpoint, corner
   along one straight edge - is preferred, and among those the lowest-id one.

   The exception is not a nicety. Fanning from either END of such an edge makes
   the triangle (corner, midpoint, corner), which is three points on a line: a
   tetrahedron on it has no volume and finley refuses the mesh. Dropping the
   sliver instead would drop the midpoint with it, and the fine side's node would
   then be a vertex of nothing on this face - which is the crack the midpoint was
   there to close. Fanning FROM the midpoint has neither problem: its two
   collinear neighbours end up at opposite ends of the fan.

   `isMid` marks those positions. Both octants sharing a face see the same
   vertices with the same ids and the same geometry, so both make the same
   choice - which is all conformity asks. (The test is structural here because
   the polygon was built by inserting them, but it is exactly the geometric
   one: a vertex whose two neighbours are collinear with it.)
*/
int fanPolygon(const std::vector<long>& gid, const long* poly, const bool* isMid,
               int n, std::array<long,3>* tris)
{
    int p = -1;
    for (int i = 0; i < n; ++i) {
        if (isMid != NULL && !isMid[i])
            continue;
        if (p < 0 || gid[poly[i]] < gid[poly[p]])
            p = i;
    }
    if (p < 0) {
        p = 0;
        for (int i = 1; i < n; ++i)
            if (gid[poly[i]] < gid[poly[p]])
                p = i;
    }
    for (int i = 0; i < n - 2; ++i) {
        tris[i][0] = poly[p];
        tris[i][1] = poly[(p + 1 + i) % n];
        tris[i][2] = poly[(p + 2 + i) % n];
    }
    return n - 2;
}

/// The triangles one face of an octant is cut into, wound counter-clockwise as
/// seen from OUTSIDE. At most eight.
int faceTriangles(const MeshAccess& m, const std::vector<long>& gid, long e,
                  int f, std::array<long,3>* tris)
{
    const int V = m.nodesPerElement;
    const long* en = &m.elementNodes[(size_t) e * V];
    const long centre = m.elementFaceHangingNode.empty() ? -1
                      : m.elementFaceHangingNode[(size_t) e * 6 + f];

    if (centre >= 0) {
        int n = 0;
        for (int k = 0; k < 4; ++k) {
            const int c = brickFaceCorners[f][k];
            const int prev = brickFaceCorners[f][(k + 3) & 3];
            const int next = brickFaceCorners[f][(k + 1) & 3];
            const long a = edgeNode(m, e, c, next);
            const long b = edgeNode(m, e, prev, c);
            if (a < 0 || b < 0)
                throw OxleyException("faceTriangles: a hanging face has an edge "
                        "with no midpoint. The finer octants across it subdivide "
                        "all four, so the seam list has missed one.");
            const long quad[4] = { en[c], a, centre, b };
            n += fanPolygon(gid, quad, NULL, 4, tris + n);
        }
        return n;
    }

    long poly[8];
    bool isMid[8];
    int np = 0;
    for (int k = 0; k < 4; ++k) {
        const int c = brickFaceCorners[f][k];
        const int next = brickFaceCorners[f][(k + 1) & 3];
        isMid[np] = false;
        poly[np++] = en[c];
        const long mid = edgeNode(m, e, c, next);
        if (mid >= 0) {
            isMid[np] = true;
            poly[np++] = mid;
        }
    }
    return fanPolygon(gid, poly, isMid, np, tris);
}

/// One tetrahedron of a cone: an apex over an outward-wound boundary triangle.
/// The orientation is fixed by MEASURING the signed volume, not by reasoning
/// about finley's winding convention.
inline void coneTet(const MeshAccess& m, long apex, const std::array<long,3>& tri,
                    std::array<long,4>& out)
{
    out[0] = apex; out[1] = tri[0]; out[2] = tri[1]; out[3] = tri[2];
    if (signedVolume(m.nodeCoords, out[0], out[1], out[2], out[3]) < 0.)
        std::swap(out[2], out[3]);
}

/**
   The tetrahedra one octant splits into. Returns how many, at most 48.

   A CONFORMING octant is coned from its lowest-id corner over the three faces
   that do not contain it: six tetrahedra, and the three faces that DO contain
   the apex are covered by the walls of the cone, which fan through the apex -
   which is the same triangulation they would get anyway, the apex being their
   lowest vertex too.

   A HANGING octant cannot be coned from any corner. A hanging face's
   triangulation is fixed by the finer octants across it and is NOT a fan
   through one of its vertices, so a corner apex lying on such a face would
   leave the part of the octant under the other sub-quads uncovered - and with
   all six faces hanging, every corner lies on one. The octant is convex, so its
   centre always serves; that is the whole reason elementCentreNode exists.
*/
int octantSplit(const MeshAccess& m, const std::vector<long>& gid, long e,
                std::array<long,4>* tets)
{
    const int V = m.nodesPerElement;
    const long* en = &m.elementNodes[(size_t) e * V];
    const long centre = m.elementCentreNode.empty() ? -1
                      : m.elementCentreNode[e];
    std::array<long,3> tris[8];
    int tet = 0;

    if (centre < 0) {
        int apex = 0;
        for (int c = 1; c < 8; ++c)
            if (gid[en[c]] < gid[en[apex]])
                apex = c;
        for (int f = 0; f < 6; ++f) {
            if (cornerOnFace(apex, f))
                continue;                    // near face, the cone covers it
            const int nt = faceTriangles(m, gid, e, f, tris);
            for (int t = 0; t < nt; ++t)
                coneTet(m, en[apex], tris[t], tets[tet++]);
        }
        return tet;
    }

    for (int f = 0; f < 6; ++f) {
        const int nt = faceTriangles(m, gid, e, f, tris);
        for (int t = 0; t < nt; ++t)
            coneTet(m, centre, tris[t], tets[tet++]);
    }
    return tet;
}

/// The triangles a boundary quad splits into, in the order the export gives
/// them the ids base + key*maxFaceSimplices + t. The octant's own split cuts
/// that face exactly the same way, so the boundary elements sit on the faces of
/// the tetrahedra rather than across them.
int boundaryTriangles(const MeshAccess& m, const std::vector<long>& gid, long f,
                      std::array<long,3>* tris)
{
    const long e = m.faceElements.empty() ? 0 : m.faceElements[f];
    const int dir = (int) (m.faceDirections.empty() ? 0 : m.faceDirections[f]);
    if (!m.elementFaceHangingNode.empty()
            && m.elementFaceHangingNode[(size_t) e * 6 + dir] >= 0)
        throw OxleyException("boundaryTriangles: a face on the domain boundary "
                "is hanging, so it has a finer neighbour - and then it is not "
                "on the boundary.");
    return faceTriangles(m, gid, e, dir, tris);
}

/**
   Appends the octant centres the hanging split cones from, one per octant that
   has any hanging node on it, and records them in elementCentreNode.

   These belong to the export alone. A centre lies strictly inside one octant,
   so no other octant and no other rank ever names it, and its id comes straight
   from the slot the numbering reserves for exactly this - no agreement needed,
   unlike a seam node.

   They are recorded as constrained on the octant's eight corners. That gives the
   tag inheritance and the nodal transfer the rule they already apply to a seam
   node, and it is the right rule here: the centre of a trilinear field is the
   mean of its corners, so a field the export can represent is exact there.
*/
void materialiseOctantCentres(MeshAccess& m, int rank)
{
    m.elementCentreNode.assign(m.numElements, -1);
    if (m.numDim != 3 || m.elementFaceHangingNode.empty())
        return;
    if ((long) m.nodeFinleyId.size() != m.numNodes)
        throw OxleyException("materialiseOctantCentres: the mesh view has no "
                "export numbering, so a centre node cannot be given an id.");

    // widen the constraint arrays: a seam node has at most four masters, an
    // octant centre has eight
    const int old = m.mastersPerConstrainedNode;
    const int mpc = 8;
    if (old < mpc) {
        const long n = (long) m.constrainedNodes.size();
        std::vector<long> masters((size_t) n * mpc, -1);
        std::vector<double> weights((size_t) n * mpc, 0.);
        for (long i = 0; i < n; ++i) {
            for (int k = 0; k < old; ++k) {
                masters[(size_t) i * mpc + k] =
                        m.constraintMasters[(size_t) i * old + k];
                weights[(size_t) i * mpc + k] =
                        m.constraintWeights[(size_t) i * old + k];
            }
        }
        m.constraintMasters.swap(masters);
        m.constraintWeights.swap(weights);
        m.mastersPerConstrainedNode = mpc;
    }

    const int V = m.nodesPerElement;
    const long slots = slotsPerElement(3);
    const long block = m.finleyDistribution[rank] + m.numOwnedNodes;
    for (long e = 0; e < m.numElements; ++e) {
        bool hanging = false;
        for (int f = 0; f < 6 && !hanging; ++f)
            hanging = m.elementFaceHangingNode[(size_t) e * 6 + f] >= 0;
        for (int ed = 0; ed < 12 && !hanging; ++ed)
            hanging = !m.elementEdgeHangingNode.empty()
                   && m.elementEdgeHangingNode[(size_t) e * 12 + ed] >= 0;
        if (!hanging)
            continue;

        const long* en = &m.elementNodes[(size_t) e * V];
        long masters[8];
        double mid[3] = {0., 0., 0.};
        for (int c = 0; c < 8; ++c) {
            masters[c] = en[c];
            for (int d = 0; d < 3; ++d)
                mid[d] += 0.125 * m.nodeCoords[(size_t) en[c] * 3 + d];
        }

        const long ni = m.numNodes++;
        for (int d = 0; d < 3; ++d)
            m.nodeCoords.push_back(mid[d]);
        m.nodeLnodesId.push_back(-1);
        if (!m.nodeTags.empty())
            m.nodeTags.push_back(inheritedTag(m, masters, 8));
        m.constrainedNodes.push_back(ni);
        for (int c = 0; c < 8; ++c) {
            m.constraintMasters.push_back(masters[c]);
            m.constraintWeights.push_back(0.125);
        }
        if (!m.hangingWriterRank.empty())
            m.hangingWriterRank.push_back(rank);
        m.nodeFinleyId.push_back(block + slots * e + OCTANT_CENTRE_SLOT_3D);
        m.elementCentreNode[e] = ni;
    }
}

/// coordinates of local node n
inline void nodePos(const MeshAccess& m, long n, double p[3])
{
    for (int d = 0; d < 3; ++d)
        p[d] = m.nodeCoords[(size_t) n * 3 + d];
}

/**
   Weights of the tensor-product interpolant through sample points that form a
   2x2x2 grid in the coordinate axes: w[k] is what sample k contributes at p.

   An axis whose points do not spread - the one normal to a boundary face -
   drops out, which is what lets the same routine do a face and an octant.
*/
void gridWeights(const double* pts, int numPts, const double* p,
                 std::vector<double>& w)
{
    w.assign(numPts, 1.);
    double lo[3], hi[3], widest = 0.;
    for (int d = 0; d < 3; ++d) {
        lo[d] = hi[d] = pts[d];
        for (int k = 1; k < numPts; ++k) {
            lo[d] = std::min(lo[d], pts[k*3 + d]);
            hi[d] = std::max(hi[d], pts[k*3 + d]);
        }
        widest = std::max(widest, hi[d] - lo[d]);
    }
    for (int d = 0; d < 3; ++d) {
        const double span = hi[d] - lo[d];
        if (span <= 1e-12 * widest)
            continue;                        // the axis normal to a face
        const double l0 = (hi[d] - p[d]) / span;
        const double l1 = (p[d] - lo[d]) / span;
        for (int k = 0; k < numPts; ++k) {
            const double t = pts[k*3 + d];
            w[k] *= (std::fabs(t - lo[d]) <= std::fabs(t - hi[d])) ? l0 : l1;
        }
    }
}

/**
   Weights of the affine interpolant through n points (4 in a tetrahedron, 3 in
   a triangle): the barycentric coordinates of p with respect to them.

   The three-point case solves the 2x2 Gram system instead of a 3x3, which is
   what makes it work on a triangle that lies in no coordinate plane - and it
   projects a point off the plane onto it, which is right here, since the values
   only ever describe a function ON that triangle.
*/
void affineWeights(const double* v, int n, const double* p, double* w)
{
    double e1[3], e2[3], e3[3], r[3];
    for (int d = 0; d < 3; ++d) {
        e1[d] = v[3 + d] - v[d];
        e2[d] = v[6 + d] - v[d];
        r[d]  = p[d] - v[d];
    }
    if (n == 4) {
        for (int d = 0; d < 3; ++d)
            e3[d] = v[9 + d] - v[d];
        const double det =
              e1[0]*(e2[1]*e3[2] - e2[2]*e3[1])
            - e2[0]*(e1[1]*e3[2] - e1[2]*e3[1])
            + e3[0]*(e1[1]*e2[2] - e1[2]*e2[1]);
        const double d1 =
              r[0]*(e2[1]*e3[2] - e2[2]*e3[1])
            - e2[0]*(r[1]*e3[2] - r[2]*e3[1])
            + e3[0]*(r[1]*e2[2] - r[2]*e2[1]);
        const double d2 =
              e1[0]*(r[1]*e3[2] - r[2]*e3[1])
            - r[0]*(e1[1]*e3[2] - e1[2]*e3[1])
            + e3[0]*(e1[1]*r[2] - e1[2]*r[1]);
        const double d3 =
              e1[0]*(e2[1]*r[2] - e2[2]*r[1])
            - e2[0]*(e1[1]*r[2] - e1[2]*r[1])
            + r[0]*(e1[1]*e2[2] - e1[2]*e2[1]);
        w[1] = d1 / det;
        w[2] = d2 / det;
        w[3] = d3 / det;
        w[0] = 1. - w[1] - w[2] - w[3];
        return;
    }
    const double a = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
    const double b = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
    const double c = e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2];
    const double u = e1[0]*r[0] + e1[1]*r[1] + e1[2]*r[2];
    const double t = e2[0]*r[0] + e2[1]*r[1] + e2[2]*r[2];
    const double det = a*c - b*b;
    w[1] = (c*u - b*t) / det;
    w[2] = (a*t - b*u) / det;
    w[0] = 1. - w[1] - w[2];
}

/// how deep p sits in a simplex given by its vertices: negative means outside,
/// and the largest value over the simplices of an octant names the one to read
inline double depthIn(const double* v, int n, const double* p)
{
    double w[4];
    affineWeights(v, n, p, w);
    double d = w[0];
    for (int k = 1; k < n; ++k)
        d = std::min(d, w[k]);
    return d;
}

/**
   For each local element, how many simplices the split emits and what fraction
   of the element each one covers.

   Recomputed rather than recorded: the split is a deterministic function of the
   hanging configuration and the element geometry, both of which MeshAccess
   already carries, so the weights can be worked out on the oxley side alone.
   That is what lets the inbound transfer weight by area without the finley side
   ever having to send areas along.

   The fractions sum to one per element, so a field that is constant over an
   element comes back unchanged whatever the split.
*/
void simplexWeights(const MeshAccess& m, const std::vector<long>& gid,
                    std::vector<int>& count,
                    std::vector<std::vector<double> >& weight)
{
    count.assign(m.numElements, 0);
    weight.assign(m.numElements, std::vector<double>());

    if (m.numDim == 3) {
        // Six tetrahedra of equal volume on a conforming octant, and anything
        // from twelve to forty-eight unequal ones on a hanging one - measured
        // rather than assumed, so the identity survives a change of split.
        std::array<long,4> tets[48];
        for (long e = 0; e < m.numElements; ++e) {
            const int nt = octantSplit(m, gid, e, tets);
            count[e] = nt;
            weight[e].resize(nt);
            double total = 0.;
            for (int t = 0; t < nt; ++t) {
                const double v = std::fabs(signedVolume(m.nodeCoords,
                        tets[t][0], tets[t][1], tets[t][2], tets[t][3]));
                weight[e][t] = v;
                total += v;
            }
            if (total <= 0.)
                throw OxleyException("simplexWeights: an octant has no volume.");
            for (int t = 0; t < nt; ++t)
                weight[e][t] /= total;
        }
        return;
    }

    const int V = m.nodesPerElement;
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

/// Where the face ids start: above every element id, which reserve
/// maxSimplices per octant. Derived from the octant count so that both the
/// converter and a later transfer arrive at the same number. Collective.
long faceIdBase(const MeshAccess& m, const escript::JMPI& mpi)
{
    long octants = m.numElements;
#ifdef ESYS_MPI
    if (mpi->size > 1) {
        long local = octants;
        MPI_Allreduce(&local, &octants, 1, MPI_LONG, MPI_SUM, mpi->comm);
    }
#endif
    return octants * maxSimplices(m.numDim);
}


/// The triangles of one element, as vertex coordinates. Same split the export
/// emitted, recomputed from the mesh view rather than stored.
void splitOfElement(const MeshAccess& m, long e,
                    std::vector<std::array<double,6> >& tris)
{
    const int V = m.nodesPerElement;
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
    tris.clear();
    for (int t = 0; t < pat->numTriangles; ++t) {
        std::array<double,6> v;
        for (int k = 0; k < 3; ++k) {
            const int sym = pat->tri[t][k];
            const long n = (sym < 4) ? V4[(sym + rot) & 3] : H4[(sym - 4 + rot) & 3];
            v[k*2]   = m.nodeCoords[(size_t) n * 2];
            v[k*2+1] = m.nodeCoords[(size_t) n * 2 + 1];
        }
        tris.push_back(v);
    }
}

/// the three points finley puts on a Tri3: its edge midpoints. Measured, on
/// conforming and graded meshes alike. Returned sorted by coordinate, which is
/// the order both sides agree on without exchanging anything.
void triangleQuadPoints(const std::array<double,6>& v,
                        std::vector<std::array<double,2> >& p)
{
    p.resize(3);
    for (int k = 0; k < 3; ++k) {
        const int a = k, b = (k + 1) % 3;
        p[k][0] = 0.5 * (v[a*2]   + v[b*2]);
        p[k][1] = 0.5 * (v[a*2+1] + v[b*2+1]);
    }
    std::sort(p.begin(), p.end(), [](const std::array<double,2>& a,
                                     const std::array<double,2>& b) {
        if (std::fabs(a[0] - b[0]) > 1e-12) return a[0] < b[0];
        return a[1] < b[1];
    });
}

/// barycentric coordinates of a point in a triangle
inline void barycentric(const std::array<double,6>& v, double x, double y,
                        double b[3])
{
    const double d = (v[2]-v[0])*(v[5]-v[1]) - (v[4]-v[0])*(v[3]-v[1]);
    b[1] = ((x-v[0])*(v[5]-v[1]) - (y-v[1])*(v[4]-v[0])) / d;
    b[2] = ((y-v[1])*(v[2]-v[0]) - (x-v[0])*(v[3]-v[1])) / d;
    b[0] = 1. - b[1] - b[2];
}

/// 1D Lagrange weights through the two Gauss abscissae of the element's rule,
/// so that evaluating at t reproduces any linear function of t exactly
inline void gaussLagrange(double t, double& l0, double& l1)
{
    const double g = 0.5 - 0.5 / std::sqrt(3.);       // 0.2113248654
    l0 = (1. - g - t) / (1. - 2. * g);
    l1 = (t - g) / (1. - 2. * g);
}

// ---------------------------------------------------------------------------
// Complex data needs no separate transfer.
//
// std::complex<double> is two doubles in memory, and every operation these
// transfers perform - copying a value, averaging masters, weighting by area,
// evaluating a shape function - is REAL-LINEAR, so it applies to the real and
// imaginary parts independently. A complex field is therefore just a real one
// with twice as many components, and the code below says so once rather than
// branching everywhere.
// ---------------------------------------------------------------------------

/// components counted as reals: twice the data point size when complex
inline int realComponents(const escript::Data& d)
{
    return d.getDataPointSize() * (d.isComplex() ? 2 : 1);
}

inline const double* readSample(const escript::Data& d, long i)
{
    return d.isComplex()
        ? reinterpret_cast<const double*>(
              d.getSampleDataRO(i, escript::DataTypes::cplx_t(0)))
        : d.getSampleDataRO(i, (double) 0);
}

inline double* writeSample(escript::Data& d, long i)
{
    return d.isComplex()
        ? reinterpret_cast<double*>(
              d.getSampleDataRW(i, escript::DataTypes::cplx_t(0)))
        : d.getSampleDataRW(i, (double) 0);
}

/// an empty Data of the same shape and complexity as src, on the given space
inline escript::Data makeLike(const escript::Data& src,
                              const escript::FunctionSpace& fs)
{
    if (src.isComplex())
        return escript::Data(escript::DataTypes::cplx_t(0),
                             src.getDataPointShape(), fs, true);
    return escript::Data(0., src.getDataPointShape(), fs, true);
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
    if (!conforming && !simplices)
        throw OxleyException("toFinley: the forest has hanging nodes. Only the "
                "simplex split handles them; the debug Rec4/Hex8 path needs a "
                "conforming forest, and always will - one element per octant "
                "cannot resolve a 2:1 seam.");

    // With hanging nodes present the mesh view must materialise them: a hanging
    // position is then a real node with a global id, so it can be a vertex of
    // the triangles on BOTH sides of a 2:1 seam. finley cannot represent a
    // hanging node (one ReferenceElementSet per ElementFile, and escript's q/r
    // is pointwise Dirichlet, not u = (u_a+u_b)/2), so the seam has to be
    // resolved by the triangulation instead - the node becomes an ordinary free
    // degree of freedom.
    MeshAccess m = dom.getMeshAccess(!conforming);
    // and in 3D the split needs one node more per hanging octant, the interior
    // apex it is coned from; see materialiseOctantCentres()
    materialiseOctantCentres(m, dom.getMPI()->rank);
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
    // how many face elements each boundary quad became. One, except in 3D,
    // where a quad is two triangles and up to six once its edges carry hanging
    // midpoints - and the face ids below have to know which.
    std::vector<int> facesPerBoundaryQuad(m.numFaces, 1);

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
                    (index_t)((m.globalElementOffset + e) * maxSimplices(m.numDim)));
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
                // an octant reserves maxSimplices ids and uses a few of them -
                // which finley does not mind: element ids are only stored and
                // handed back, never used to size anything.
                out.elementId.push_back(
                        (index_t)((m.globalElementOffset + e) * maxSimplices(m.numDim) + t));
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
        // ---- 3D: cone each octant, and cut every face from its lowest-id
        // vertex. A conforming octant is coned from its lowest-id CORNER over
        // the three faces that do not contain it - six tetrahedra - and a
        // hanging one from its centre over all six. See octantSplit().
        out.elementType = finley::Tet4;
        out.faceElementType = finley::Tri3;

        std::array<long,4> tets[48];
        std::array<long,3> tris[8];
        elems.reserve((size_t) m.numElements * 6 * 4);
        out.elementTag.reserve(m.numElements * 6);
        for (long e = 0; e < m.numElements; ++e) {
            const int nt = octantSplit(m, gidOf, e, tets);
            for (int t = 0; t < nt; ++t) {
                for (int k = 0; k < 4; ++k)
                    elems.push_back((index_t) tets[t][k]);
                out.elementTag.push_back((int) m.elementTags[e]);
                // an id that says which octant this tetrahedron came from, so a
                // transfer can find its way back without a stored map
                out.elementId.push_back((index_t)(
                        (m.globalElementOffset + e) * maxSimplices(m.numDim) + t));
            }
        }

        // Boundary quads are cut by the routine the octant's own split uses for
        // that face, so the triangles are the same triangles - the face elements
        // sit ON the tetrahedra rather than across them. getMeshAccess() wound
        // the corners counter-clockwise seen from outside and fanning preserves
        // that, so the normals stay outward.
        faces.reserve((size_t) m.numFaces * 2 * 3);
        out.faceTag.reserve(m.numFaces * 2);
        facesPerBoundaryQuad.resize(m.numFaces);
        for (long f = 0; f < m.numFaces; ++f) {
            const int nt = boundaryTriangles(m, gidOf, f, tris);
            facesPerBoundaryQuad[f] = nt;
            for (int t = 0; t < nt; ++t) {
                for (int k = 0; k < 3; ++k)
                    faces.push_back((index_t) tris[t][k]);
                out.faceTag.push_back((int) m.faceTags[f]);
            }
        }
    }

    // Collective, and exact on any number of ranks: the interface faces are
    // settled by an exchange rather than assumed. Keyed by global node ids, so
    // this has to come before the local indices are translated below - it does
    // its own translation.
    if (simplices) {
        if (m.numDim == 3)
            checkConformity<4, 3>(elems, 4, tetFaces, faces, gidOf,
                                  dom.getMPI(), "tetrahedral");
        else
            checkConformity<3, 2>(elems, 3, triFaces, faces, gidOf,
                                  dom.getMPI(), "triangular");
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
    // Clear of the element ids, which reserve maxSimplices per OCTANT. Based
    // on the octant count rather than the simplex count so that a transfer can
    // work the same offset out from the mesh view alone - see faceIdBase().
    // No per-rank scan any more: the ids below are derived from the global
    // element index, so they are already unique across ranks. Adding a scan -
    // which is what consecutive numbering needed - shifted every rank but the
    // first, and since the scan is zero in serial the mismatch only showed
    // under MPI.
    const index_t faceIdOffset = (index_t) faceIdBase(m, mpiInfo);
    // Ids that say which boundary face of which octant this came from, the
    // same trick the elements use: (global element, face direction) names a
    // boundary face identically on every rank, so a transfer needs no stored
    // map and survives finley redistributing the mesh. The trailing slot says
    // which triangle of that quad: two in 3D, and up to six once the quad's
    // edges carry hanging midpoints and it is an octagon.
    const int perQuad = maxFaceSimplices(m.numDim);
    out.faceId.resize(numFaces);
    for (long f = 0, k = 0; f < m.numFaces; ++f) {
        const long elem = m.faceElements.empty() ? 0 : m.faceElements[f];
        const long dir = m.faceDirections.empty() ? 0 : m.faceDirections[f];
        const long key = ((m.globalElementOffset + elem) * (2 * m.numDim) + dir)
                       * perQuad;
        for (int t = 0; t < facesPerBoundaryQuad[f]; ++t, ++k)
            out.faceId[k] = faceIdOffset + (index_t)(key + t);
    }

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
/// The mesh view the export was built from, rebuilt: hanging positions
/// materialised where the forest has them, and in 3D the octant centres the
/// split cones from. It must be the SAME view toFinley() used, node for node,
/// or the ids a transfer derives name something else.
MeshAccess viewMatchingExport(const OxleyDomain& dom)
{
    MeshAccess m = dom.getMeshAccess(!dom.isConforming());   // collective
    materialiseOctantCentres(m, dom.getMPI()->rank);
    return m;
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
    checkSameForest(*dom, *target, "toFinleyData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const std::vector<long>& gid = exportIds(m);
    const int numComp = realComponents(source);

    // What this rank can supply: its own nodes, plus the positions it
    // materialised - the seam nodes and the octant centres, which are no nodes
    // of the forest and so take the average of their masters.
    //
    // A master can itself be materialised. An octant that is the COARSE side of
    // one seam and the FINE side of another has a corner which is a seam node,
    // and its centre is the average of its eight corners including that one. So
    // the values are built up over ALL local nodes in creation order, which is
    // an order in which a master always comes before what depends on it: a seam
    // node's masters are corners of a coarse octant and always real, and the
    // centres are appended after every seam node.
    std::vector<double> nodeVal((size_t) m.numNodes * numComp, 0.);
    for (long i = 0; i < m.numRealNodes; ++i) {
        const double* in = readSample(source, i);
        for (int c = 0; c < numComp; ++c)
            nodeVal[(size_t) i * numComp + c] = in[c];
    }
    const int mpc = m.mastersPerConstrainedNode;
    for (size_t k = 0; k < m.constrainedNodes.size(); ++k) {
        const long node = m.constrainedNodes[k];
        double* v = &nodeVal[(size_t) node * numComp];
        for (int j = 0; j < mpc; ++j) {
            const long master = m.constraintMasters[k*mpc + j];
            const double w = m.constraintWeights[k*mpc + j];
            if (master < 0 || w == 0.)
                continue;
            if (master >= node)
                throw OxleyException("toFinleyData: a materialised node "
                        "depends on one created after it, so its masters have "
                        "no value yet.");
            for (int c = 0; c < numComp; ++c)
                v[c] += w * nodeVal[(size_t) master * numComp + c];
        }
    }

    std::vector<long> haveId;
    std::vector<double> haveVal;
    haveId.reserve(m.numNodes);
    haveVal.reserve((size_t) m.numNodes * numComp);
    for (long i = 0; i < m.numNodes; ++i) {
        haveId.push_back(gid[i]);
        for (int c = 0; c < numComp; ++c)
            haveVal.push_back(nodeVal[(size_t) i * numComp + c]);
    }

    escript::Data result = makeLike(source, escript::continuousFunction(*target));
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
        double* out = writeSample(result, j);
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
    checkSameForest(*dom, *(source.getFunctionSpace().getDomain()),
                    "fromFinleyData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const std::vector<long>& gid = exportIds(m);
    const int numComp = realComponents(source);

    const escript::FunctionSpace sourceFS = source.getFunctionSpace();
    const long ns = (long) source.getNumSamples();
    const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
            sourceFS.getTypeCode());
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * numComp);
    for (long j = 0; j < ns; ++j) {
        haveId[j] = (long) ids[j];
        const double* in = readSample(source, j);
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

    escript::Data result = makeLike(source, escript::continuousFunction(*target));
    result.requireWrite();
    for (long i = 0; i < m.numRealNodes; ++i) {
        double* out = writeSample(result, i);
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
    checkSameForest(*dom, *target, "toFinleyReducedData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = realComponents(source);
    std::vector<int> count;
    std::vector<std::vector<double> > weight;
    simplexWeights(m, exportIds(m), count, weight);

    // one octant value, repeated onto each simplex it was split into
    std::vector<long> haveId;
    std::vector<double> haveVal;
    for (long e = 0; e < m.numElements; ++e) {
        const double* in = readSample(source, e);
        for (int t = 0; t < count[e]; ++t) {
            haveId.push_back((m.globalElementOffset + e) * maxSimplices(m.numDim) + t);
            for (int c = 0; c < numComp; ++c)
                haveVal.push_back(in[c]);
        }
    }

    escript::Data result = makeLike(source, escript::reducedFunction(*target));
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
        double* out = writeSample(result, j);
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
    checkSameForest(*dom, *(source.getFunctionSpace().getDomain()),
                    "fromFinleyReducedData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = realComponents(source);
    std::vector<int> count;
    std::vector<std::vector<double> > weight;
    simplexWeights(m, exportIds(m), count, weight);

    const escript::FunctionSpace sourceFS = source.getFunctionSpace();
    const long ns = (long) source.getNumSamples();
    const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
            sourceFS.getTypeCode());
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * numComp);
    for (long j = 0; j < ns; ++j) {
        haveId[j] = (long) ids[j];
        const double* in = readSample(source, j);
        for (int c = 0; c < numComp; ++c)
            haveVal[(size_t) j*numComp + c] = in[c];
    }

    // ask for every simplex of every octant this rank owns
    std::vector<long> wantId;
    for (long e = 0; e < m.numElements; ++e)
        for (int t = 0; t < count[e]; ++t)
            wantId.push_back((m.globalElementOffset + e) * maxSimplices(m.numDim) + t);

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numComp, haveId, haveVal, wantId, wantVal);

    // and average them by the area each covers, so a field that is constant
    // over the octant comes back unchanged whatever the split
    escript::Data result = makeLike(source, escript::reducedFunction(*target));
    result.requireWrite();
    size_t pos = 0;
    for (long e = 0; e < m.numElements; ++e) {
        double* out = writeSample(result, e);
        for (int c = 0; c < numComp; ++c)
            out[c] = 0.;
        for (int t = 0; t < count[e]; ++t, ++pos)
            for (int c = 0; c < numComp; ++c)
                out[c] += weight[e][t] * wantVal[pos*numComp + c];
    }
    return result;
}


namespace {

/// the id the export gives the boundary faces of one oxley face
inline long faceKeyOf(const MeshAccess& m, long f, long base)
{
    const long elem = m.faceElements.empty() ? 0 : m.faceElements[f];
    const long dir = m.faceDirections.empty() ? 0 : m.faceDirections[f];
    return base + ((m.globalElementOffset + elem) * (2 * m.numDim) + dir) * 2;
}

/**
   The order in which a face's quadrature points should be read.

   The two meshes put the SAME physical points on a boundary face - measured,
   for both FunctionOnBoundary and its reduced form - but not necessarily in
   the same order, because the two describe the edge with their own winding.
   Sorting each face's points by coordinate gives both sides the same order
   without either having to send coordinates: the values can then be matched
   position by position.
*/
void pointOrder(const escript::Data& x, long sample, int numPoints, int dim,
                std::vector<int>& order)
{
    const double* p = x.getSampleDataRO(sample, (double) 0);
    order.resize(numPoints);
    for (int i = 0; i < numPoints; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        for (int d = 0; d < dim; ++d) {
            const double u = p[a*dim + d], v = p[b*dim + d];
            if (std::fabs(u - v) > 1e-12)
                return u < v;
        }
        return a < b;
    });
}

/// names one octant face across both meshes: the export builds the ids of its
/// triangles as base + key*maxFaceSimplices + t, so either side can form the
/// other's
inline long faceQuadKey(const MeshAccess& m, long f)
{
    const long elem = m.faceElements.empty() ? 0 : m.faceElements[f];
    const long dir = m.faceDirections.empty() ? 0 : m.faceDirections[f];
    return (m.globalElementOffset + elem) * (2 * m.numDim) + dir;
}

/**
   The 3D boundary transfer.

   Not a permutation, unlike 2D: a boundary quad becomes two triangles, and
   neither the point counts nor the rules match. Outbound, the quad's points are
   a tensor grid in the face's two axes and so determine a bilinear function,
   which the triangles' points read off. Inbound, each triangle's three points
   determine an affine function, and each of the quad's points is read off the
   triangle that contains it.

   The reduced spaces carry a single value per sample, which no evaluation can
   improve on: outbound it is replicated (gridWeights on one point does exactly
   that), inbound the two triangles are averaged by area. Replicate followed by
   average is the identity, as for ReducedFunction.
*/
escript::Data transferBoundary3D(const escript::Data& source,
                                 escript::Domain_ptr target, bool toFinley,
                                 const MeshAccess& m, long base,
                                 const OxleyDomain& dom, int fsCode,
                                 const char* what)
{
    const int numComp = realComponents(source);
    const int srcPts = source.getNumDataPointsPerSample();
    const std::vector<long>& gid = exportIds(m);
    const int stride = srcPts * (numComp + 3);

    const escript::FunctionSpace targetFS = (fsCode == FaceElements)
            ? escript::functionOnBoundary(*target)
            : escript::reducedFunctionOnBoundary(*target);
    escript::Data result = makeLike(source, targetFS);
    result.requireWrite();
    const int dstPts = result.getNumDataPointsPerSample();
    escript::Data sx = source.getFunctionSpace().getX();
    escript::Data rx = targetFS.getX();

    // what this side can supply, values and the points they were taken at
    const long ns = (long) source.getNumSamples();
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * stride);
    const dim_t* srcIds = toFinley ? NULL
            : source.getFunctionSpace().getDomain()->borrowSampleReferenceIDs(
                    source.getFunctionSpace().getTypeCode());
    for (long j = 0; j < ns; ++j) {
        haveId[j] = toFinley ? faceQuadKey(m, j) : (long) srcIds[j];
        const double* in = readSample(source, j);
        const double* xs = sx.getSampleDataRO(j, (double) 0);
        double* out = &haveVal[(size_t) j * stride];
        for (int q = 0; q < srcPts; ++q) {
            for (int c = 0; c < numComp; ++c)
                out[q*numComp + c] = in[q*numComp + c];
            for (int d = 0; d < 3; ++d)
                out[srcPts*numComp + q*3 + d] = xs[q*3 + d];
        }
    }

    if (toFinley) {
        // each triangle asks for the quad it was cut from
        const long n = (long) result.getNumSamples();
        const dim_t* ids = target->borrowSampleReferenceIDs(
                targetFS.getTypeCode());
        std::vector<long> wantId(n);
        for (long j = 0; j < n; ++j)
            wantId[j] = ((long) ids[j] - base) / maxFaceSimplices(m.numDim);

        std::vector<double> wantVal;
        exchangeByGlobalId(dom.getMPI(), stride, haveId, haveVal,
                           wantId, wantVal);

        std::vector<double> w;
        for (long j = 0; j < n; ++j) {
            const double* got = &wantVal[(size_t) j * stride];
            const double* pts = got + srcPts*numComp;
            const double* px = rx.getSampleDataRO(j, (double) 0);
            double* out = writeSample(result, j);
            for (int q = 0; q < dstPts; ++q) {
                gridWeights(pts, srcPts, &px[q*3], w);
                for (int c = 0; c < numComp; ++c) {
                    double v = 0.;
                    for (int k = 0; k < srcPts; ++k)
                        v += w[k] * got[k*numComp + c];
                    out[q*numComp + c] = v;
                }
            }
        }
        return result;
    }

    // inbound: every triangle of every boundary quad this rank holds. How many
    // that is varies once the quad's edges carry hanging midpoints, so the
    // triangles are addressed through a prefix sum rather than a fixed stride.
    const int perQuad = maxFaceSimplices(m.numDim);
    std::array<long,3> tris[8];
    std::vector<long> triOffset(m.numFaces + 1, 0);
    std::vector<long> wantId;
    wantId.reserve((size_t) m.numFaces * 2);
    for (long f = 0; f < m.numFaces; ++f) {
        const int nt = boundaryTriangles(m, gid, f, tris);
        triOffset[f+1] = triOffset[f] + nt;
        for (int t = 0; t < nt; ++t)
            wantId.push_back(base + faceQuadKey(m, f) * perQuad + t);
    }

    std::vector<double> wantVal;
    exchangeByGlobalId(dom.getMPI(), stride, haveId, haveVal, wantId, wantVal);

    for (long f = 0; f < m.numFaces; ++f) {
        const int nt = boundaryTriangles(m, gid, f, tris);
        double vert[8][9];
        for (int t = 0; t < nt; ++t)
            for (int k = 0; k < 3; ++k)
                nodePos(m, tris[t][k], &vert[t][k*3]);

        const double* px = rx.getSampleDataRO(f, (double) 0);
        double* out = writeSample(result, f);

        if (srcPts == 1) {
            // nothing to evaluate: average the triangles by the area each covers
            double area[8], total = 0.;
            for (int t = 0; t < nt; ++t) {
                double e1[3], e2[3], cr[3];
                for (int d = 0; d < 3; ++d) {
                    e1[d] = vert[t][3+d] - vert[t][d];
                    e2[d] = vert[t][6+d] - vert[t][d];
                }
                cr[0] = e1[1]*e2[2] - e1[2]*e2[1];
                cr[1] = e1[2]*e2[0] - e1[0]*e2[2];
                cr[2] = e1[0]*e2[1] - e1[1]*e2[0];
                area[t] = 0.5 * std::sqrt(cr[0]*cr[0] + cr[1]*cr[1]
                                        + cr[2]*cr[2]);
                total += area[t];
            }
            if (total <= 0.)
                throw OxleyException(std::string(what) + ": a boundary face "
                        "has no area.");
            for (int c = 0; c < numComp; ++c)
                out[c] = 0.;
            for (int t = 0; t < nt; ++t) {
                const double* got =
                        &wantVal[(size_t)(triOffset[f] + t) * stride];
                for (int c = 0; c < numComp; ++c)
                    out[c] += (area[t]/total) * got[c];
            }
            continue;
        }

        if (srcPts != 3)
            throw OxleyException(std::string(what) + ": a Tri3 must carry "
                    "three quadrature points for its values to determine an "
                    "affine function.");

        for (int q = 0; q < dstPts; ++q) {
            const double* p = &px[q*3];
            int best = 0;
            double bestDepth = -std::numeric_limits<double>::max();
            for (int t = 0; t < nt; ++t) {
                const double d = depthIn(vert[t], 3, p);
                if (d > bestDepth) { bestDepth = d; best = t; }
            }
            const double* got =
                    &wantVal[(size_t)(triOffset[f] + best) * stride];
            const double* pts = got + srcPts*numComp;
            double w[3];
            affineWeights(pts, 3, p, w);
            for (int c = 0; c < numComp; ++c) {
                double v = 0.;
                for (int k = 0; k < 3; ++k)
                    v += w[k] * got[k*numComp + c];
                out[q*numComp + c] = v;
            }
        }
    }
    return result;
}

/// shared by both directions: move boundary values keyed by (face, point)
escript::Data transferBoundary(const escript::Data& source,
                               escript::Domain_ptr target, bool toFinley,
                               const char* what)
{
    const int fsCode = source.getFunctionSpace().getTypeCode();
    if (fsCode != FaceElements && fsCode != ReducedFaceElements)
        throw OxleyException(std::string(what) + ": the data must live on "
                "FunctionOnBoundary or ReducedFunctionOnBoundary.");

    const OxleyDomain* dom = toFinley
            ? dynamic_cast<const OxleyDomain*>(
                    source.getFunctionSpace().getDomain().get())
            : dynamic_cast<const OxleyDomain*>(target.get());
    if (dom == NULL)
        throw OxleyException(std::string(what) + ": one side must be an oxley "
                "domain.");
    const escript::AbstractDomain& other = toFinley
            ? *target : *(source.getFunctionSpace().getDomain());
    checkSameForest(*dom, other, what);

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const long base = faceIdBase(m, dom->getMPI());         // collective
    const int dim = m.numDim;
    const int numComp = realComponents(source);

    if (dim == 3)
        return transferBoundary3D(source, target, toFinley, m, base, *dom,
                                  fsCode, what);

    // the same space on the other side
    escript::FunctionSpace targetFS = (fsCode == FaceElements)
            ? escript::functionOnBoundary(toFinley ? *target : *target)
            : escript::reducedFunctionOnBoundary(toFinley ? *target : *target);
    escript::Data result = makeLike(source, targetFS);
    result.requireWrite();

    const int numPoints = source.getNumDataPointsPerSample();
    if (result.getNumDataPointsPerSample() != numPoints)
        throw OxleyException(std::string(what) + ": the two meshes disagree "
                "about the number of quadrature points on a boundary face.");

    // keys: one per face element on each side, sorted point order within it
    escript::Data sx = source.getFunctionSpace().getX();
    escript::Data rx = targetFS.getX();

    std::vector<long> haveId, wantId;
    std::vector<double> haveVal;
    std::vector<int> order;

    const long ns = (long) source.getNumSamples();
    haveId.reserve(ns);
    haveVal.reserve((size_t) ns * numPoints * numComp);
    if (toFinley) {
        // supplying from the forest: one key per oxley boundary face
        for (long f = 0; f < ns; ++f) {
            haveId.push_back(faceKeyOf(m, f, base));
            pointOrder(sx, f, numPoints, dim, order);
            const double* in = readSample(source, f);
            for (int i = 0; i < numPoints; ++i)
                for (int c = 0; c < numComp; ++c)
                    haveVal.push_back(in[order[i]*numComp + c]);
        }
    } else {
        // supplying from the export: the face element ids say which face
        const dim_t* ids = source.getFunctionSpace().getDomain()
                ->borrowSampleReferenceIDs(fsCode);
        for (long f = 0; f < ns; ++f) {
            haveId.push_back((long) ids[f]);
            pointOrder(sx, f, numPoints, dim, order);
            const double* in = readSample(source, f);
            for (int i = 0; i < numPoints; ++i)
                for (int c = 0; c < numComp; ++c)
                    haveVal.push_back(in[order[i]*numComp + c]);
        }
    }

    const long nr = (long) result.getNumSamples();
    wantId.reserve(nr);
    if (toFinley) {
        const dim_t* ids = target->borrowSampleReferenceIDs(
                targetFS.getTypeCode());
        for (long f = 0; f < nr; ++f)
            wantId.push_back((long) ids[f]);
    } else {
        for (long f = 0; f < nr; ++f)
            wantId.push_back(faceKeyOf(m, f, base));
    }

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), numPoints * numComp,
                       haveId, haveVal, wantId, wantVal);

    for (long f = 0; f < nr; ++f) {
        pointOrder(rx, f, numPoints, dim, order);
        double* out = writeSample(result, f);
        for (int i = 0; i < numPoints; ++i)
            for (int c = 0; c < numComp; ++c)
                out[order[i]*numComp + c] =
                        wantVal[((size_t) f*numPoints + i)*numComp + c];
    }
    return result;
}

} // anonymous namespace

escript::Data toFinleyBoundaryData(const escript::Data& source,
                                   escript::Domain_ptr target)
{
    if (target.get() == NULL)
        throw OxleyException("toFinleyBoundaryData: no target domain given.");
    return transferBoundary(source, target, true, "toFinleyBoundaryData");
}

escript::Data fromFinleyBoundaryData(const escript::Data& source,
                                     escript::Domain_ptr target)
{
    if (target.get() == NULL)
        throw OxleyException("fromFinleyBoundaryData: no target domain given.");
    return transferBoundary(source, target, false, "fromFinleyBoundaryData");
}


escript::Data toFinleyFunctionData(const escript::Data& source,
                                   escript::Domain_ptr target)
{
    if (source.getFunctionSpace().getTypeCode() != Elements)
        throw OxleyException("toFinleyFunctionData: the data must live on "
                "Function.");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(
            source.getFunctionSpace().getDomain().get());
    if (dom == NULL)
        throw OxleyException("toFinleyFunctionData: the source must live on an "
                "oxley domain.");
    if (target.get() == NULL)
        throw OxleyException("toFinleyFunctionData: no target domain given.");
    checkSameForest(*dom, *target, "toFinleyFunctionData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = realComponents(source);
    const int srcPts = source.getNumDataPointsPerSample();

    if (m.numDim == 3) {
        // One message per OCTANT, not per tet: a tet's id already says which
        // octant it came from, so all six ask for the same key and the buffer
        // stays the size of the octant count.
        escript::Data sx = source.getFunctionSpace().getX();
        const int stride = srcPts * (numComp + 3);
        std::vector<long> haveId(m.numElements);
        std::vector<double> haveVal((size_t) m.numElements * stride);
        for (long e = 0; e < m.numElements; ++e) {
            haveId[e] = m.globalElementOffset + e;
            const double* in = readSample(source, e);
            const double* xs = sx.getSampleDataRO(e, (double) 0);
            double* out = &haveVal[(size_t) e * stride];
            for (int q = 0; q < srcPts; ++q) {
                for (int c = 0; c < numComp; ++c)
                    out[q*numComp + c] = in[q*numComp + c];
                for (int d = 0; d < 3; ++d)
                    out[srcPts*numComp + q*3 + d] = xs[q*3 + d];
            }
        }

        escript::Data result = makeLike(source, escript::function(*target));
        result.requireWrite();
        const escript::FunctionSpace targetFS = escript::function(*target);
        const long n = (long) result.getNumSamples();
        const int dstPts = result.getNumDataPointsPerSample();
        const dim_t* ids = target->borrowSampleReferenceIDs(
                targetFS.getTypeCode());
        std::vector<long> wantId(n);
        for (long j = 0; j < n; ++j)
            wantId[j] = (long) ids[j] / maxSimplices(m.numDim);  // the parent octant

        std::vector<double> wantVal;
        exchangeByGlobalId(dom->getMPI(), stride, haveId, haveVal,
                           wantId, wantVal);

        escript::Data rx = targetFS.getX();
        std::vector<double> w;
        for (long j = 0; j < n; ++j) {
            const double* got = &wantVal[(size_t) j * stride];
            const double* pts = got + srcPts*numComp;
            const double* px = rx.getSampleDataRO(j, (double) 0);
            double* out = writeSample(result, j);
            for (int q = 0; q < dstPts; ++q) {
                gridWeights(pts, srcPts, &px[q*3], w);
                for (int c = 0; c < numComp; ++c) {
                    double v = 0.;
                    for (int k = 0; k < srcPts; ++k)
                        v += w[k] * got[k*numComp + c];
                    out[q*numComp + c] = v;
                }
            }
        }
        return result;
    }

    if (srcPts != 4)
        throw OxleyException("toFinleyFunctionData: expected the 2x2 Gauss "
                "rule on the octants.");

    // The four values of an octant determine one bilinear function - the
    // Gauss abscissae are unisolvent for it - so evaluating that at the
    // triangles' points is exact for anything bilinear, linear included.
    std::vector<long> haveId;
    std::vector<double> haveVal;
    std::vector<std::array<double,6> > tris;
    std::vector<std::array<double,2> > pts;
    for (long e = 0; e < m.numElements; ++e) {
        const double* in = readSample(source, e);
        // the octant's own frame, from its corners
        const long* en = &m.elementNodes[(size_t) e * m.nodesPerElement];
        const double x0 = m.nodeCoords[(size_t) en[0] * 2];
        const double y0 = m.nodeCoords[(size_t) en[0] * 2 + 1];
        const double h = m.nodeCoords[(size_t) en[3] * 2] - x0;

        splitOfElement(m, e, tris);
        for (size_t t = 0; t < tris.size(); ++t) {
            triangleQuadPoints(tris[t], pts);
            haveId.push_back((m.globalElementOffset + e) * maxSimplices(m.numDim)
                             + (long) t);
            for (int q = 0; q < 3; ++q) {
                double lx0, lx1, ly0, ly1;
                gaussLagrange((pts[q][0] - x0) / h, lx0, lx1);
                gaussLagrange((pts[q][1] - y0) / h, ly0, ly1);
                const double w[4] = { lx0*ly0, lx1*ly0, lx0*ly1, lx1*ly1 };
                for (int c = 0; c < numComp; ++c) {
                    double v = 0.;
                    for (int k = 0; k < 4; ++k)
                        v += w[k] * in[k*numComp + c];
                    haveVal.push_back(v);
                }
            }
        }
    }

    escript::Data result = makeLike(source, escript::function(*target));
    result.requireWrite();
    const escript::FunctionSpace targetFS = escript::function(*target);
    if (result.getNumDataPointsPerSample() != 3)
        throw OxleyException("toFinleyFunctionData: expected three points on a "
                "Tri3.");
    const long n = (long) result.getNumSamples();
    const dim_t* ids = target->borrowSampleReferenceIDs(targetFS.getTypeCode());
    std::vector<long> wantId(n);
    for (long j = 0; j < n; ++j)
        wantId[j] = (long) ids[j];

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), 3 * numComp, haveId, haveVal,
                       wantId, wantVal);

    escript::Data rx = targetFS.getX();
    std::vector<int> order;
    for (long j = 0; j < n; ++j) {
        pointOrder(rx, j, 3, 2, order);
        double* out = writeSample(result, j);
        for (int q = 0; q < 3; ++q)
            for (int c = 0; c < numComp; ++c)
                out[order[q]*numComp + c] =
                        wantVal[((size_t) j*3 + q)*numComp + c];
    }
    return result;
}

escript::Data fromFinleyFunctionData(const escript::Data& source,
                                     escript::Domain_ptr target)
{
    if (source.getFunctionSpace().getTypeCode() != Elements)
        throw OxleyException("fromFinleyFunctionData: the data must live on "
                "Function.");
    const OxleyDomain* dom = dynamic_cast<const OxleyDomain*>(target.get());
    if (dom == NULL)
        throw OxleyException("fromFinleyFunctionData: the target must be an "
                "oxley domain.");
    checkSameForest(*dom, *(source.getFunctionSpace().getDomain()),
                    "fromFinleyFunctionData");

    const MeshAccess m = viewMatchingExport(*dom);          // collective
    const int numComp = realComponents(source);

    if (m.numDim == 3) {
        const int srcPts = source.getNumDataPointsPerSample();
        if (srcPts != 4)
            throw OxleyException("fromFinleyFunctionData: a Tet4 must carry "
                    "four quadrature points for its values to determine an "
                    "affine function.");

        const escript::FunctionSpace sourceFS = source.getFunctionSpace();
        escript::Data sx = sourceFS.getX();
        const long ns = (long) source.getNumSamples();
        const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
                sourceFS.getTypeCode());
        const int stride = srcPts * (numComp + 3);
        std::vector<long> haveId(ns);
        std::vector<double> haveVal((size_t) ns * stride);
        for (long j = 0; j < ns; ++j) {
            haveId[j] = (long) ids[j];
            const double* in = readSample(source, j);
            const double* xs = sx.getSampleDataRO(j, (double) 0);
            double* out = &haveVal[(size_t) j * stride];
            for (int q = 0; q < srcPts; ++q) {
                for (int c = 0; c < numComp; ++c)
                    out[q*numComp + c] = in[q*numComp + c];
                for (int d = 0; d < 3; ++d)
                    out[srcPts*numComp + q*3 + d] = xs[q*3 + d];
            }
        }

        // every tet of every octant this rank owns. How many that is varies
        // with the octant's hanging configuration, so they are addressed
        // through a prefix sum rather than a fixed stride.
        const std::vector<long>& gid = exportIds(m);
        std::array<long,4> tets[48];
        std::vector<long> tetOffset(m.numElements + 1, 0);
        std::vector<long> wantId;
        wantId.reserve((size_t) m.numElements * 6);
        for (long e = 0; e < m.numElements; ++e) {
            const int nt = octantSplit(m, gid, e, tets);
            tetOffset[e+1] = tetOffset[e] + nt;
            for (int t = 0; t < nt; ++t)
                wantId.push_back((m.globalElementOffset + e)
                                 * maxSimplices(m.numDim) + t);
        }

        std::vector<double> wantVal;
        exchangeByGlobalId(dom->getMPI(), stride, haveId, haveVal,
                           wantId, wantVal);

        escript::Data result = makeLike(source, escript::function(*target));
        result.requireWrite();
        const escript::FunctionSpace targetFS = escript::function(*target);
        const int dstPts = result.getNumDataPointsPerSample();
        escript::Data rx = targetFS.getX();
        for (long e = 0; e < m.numElements; ++e) {
            const int nt = octantSplit(m, gid, e, tets);
            std::vector<double> vert((size_t) nt * 12);
            for (int t = 0; t < nt; ++t)
                for (int k = 0; k < 4; ++k)
                    nodePos(m, tets[t][k], &vert[(size_t) t*12 + k*3]);

            const double* px = rx.getSampleDataRO(e, (double) 0);
            double* out = writeSample(result, e);
            for (int q = 0; q < dstPts; ++q) {
                const double* p = &px[q*3];
                // the tet this Gauss point sits in, decided against the tets'
                // VERTICES: their quadrature points span only part of them
                int best = 0;
                double bestDepth = -std::numeric_limits<double>::max();
                for (int t = 0; t < nt; ++t) {
                    const double d = depthIn(&vert[(size_t) t*12], 4, p);
                    if (d > bestDepth) { bestDepth = d; best = t; }
                }
                const double* got =
                        &wantVal[(size_t)(tetOffset[e] + best) * stride];
                const double* pts = got + srcPts*numComp;
                double w[4];
                affineWeights(pts, 4, p, w);
                for (int c = 0; c < numComp; ++c) {
                    double v = 0.;
                    for (int k = 0; k < 4; ++k)
                        v += w[k] * got[k*numComp + c];
                    out[q*numComp + c] = v;
                }
            }
        }
        return result;
    }

    if (source.getNumDataPointsPerSample() != 3)
        throw OxleyException("fromFinleyFunctionData: expected three points on "
                "a Tri3.");

    // each triangle supplies its three values, in its own sorted point order
    const escript::FunctionSpace sourceFS = source.getFunctionSpace();
    const long ns = (long) source.getNumSamples();
    const dim_t* ids = sourceFS.getDomain()->borrowSampleReferenceIDs(
            sourceFS.getTypeCode());
    escript::Data sx = sourceFS.getX();
    std::vector<long> haveId(ns);
    std::vector<double> haveVal((size_t) ns * 3 * numComp);
    std::vector<int> order;
    for (long j = 0; j < ns; ++j) {
        haveId[j] = (long) ids[j];
        pointOrder(sx, j, 3, 2, order);
        const double* in = readSample(source, j);
        for (int q = 0; q < 3; ++q)
            for (int c = 0; c < numComp; ++c)
                haveVal[((size_t) j*3 + q)*numComp + c] =
                        in[order[q]*numComp + c];
    }

    std::vector<std::array<double,6> > tris;
    std::vector<std::array<double,2> > pts;
    std::vector<long> wantId;
    std::vector<int> perElement(m.numElements, 0);
    for (long e = 0; e < m.numElements; ++e) {
        splitOfElement(m, e, tris);
        perElement[e] = (int) tris.size();
        for (size_t t = 0; t < tris.size(); ++t)
            wantId.push_back((m.globalElementOffset + e) * maxSimplices(m.numDim)
                             + (long) t);
    }

    std::vector<double> wantVal;
    exchangeByGlobalId(dom->getMPI(), 3 * numComp, haveId, haveVal,
                       wantId, wantVal);

    // Rebuild the linear function on each triangle - three edge midpoints are
    // unisolvent for it - and read the octant's own Gauss points off whichever
    // triangle contains them. Exact for a field that is linear on the split.
    escript::Data result = makeLike(source, escript::function(*target));
    result.requireWrite();
    const double g = 0.5 - 0.5 / std::sqrt(3.);
    size_t pos = 0;
    for (long e = 0; e < m.numElements; ++e) {
        splitOfElement(m, e, tris);
        const long* en = &m.elementNodes[(size_t) e * m.nodesPerElement];
        const double x0 = m.nodeCoords[(size_t) en[0] * 2];
        const double y0 = m.nodeCoords[(size_t) en[0] * 2 + 1];
        const double h = m.nodeCoords[(size_t) en[3] * 2] - x0;

        // vertex values of each triangle, from its midpoint values
        std::vector<std::vector<double> > vertexVal(tris.size());
        for (size_t t = 0; t < tris.size(); ++t, ++pos) {
            triangleQuadPoints(tris[t], pts);
            vertexVal[t].assign(3 * numComp, 0.);
            for (int k = 0; k < 3; ++k) {
                // which sorted midpoint is which edge
                const int a = k, b = (k + 1) % 3;
                const double mx = 0.5*(tris[t][a*2] + tris[t][b*2]);
                const double my = 0.5*(tris[t][a*2+1] + tris[t][b*2+1]);
                int which = 0;
                for (int q = 0; q < 3; ++q)
                    if (std::fabs(pts[q][0]-mx) < 1e-12
                            && std::fabs(pts[q][1]-my) < 1e-12)
                        which = q;
                for (int c = 0; c < numComp; ++c)
                    vertexVal[t][k*numComp + c] =
                            wantVal[(pos*3 + which)*numComp + c];
            }
            // midpoint values m_ab, m_bc, m_ca -> vertex values
            for (int c = 0; c < numComp; ++c) {
                const double mab = vertexVal[t][0*numComp + c];
                const double mbc = vertexVal[t][1*numComp + c];
                const double mca = vertexVal[t][2*numComp + c];
                vertexVal[t][0*numComp + c] = mab + mca - mbc;   // a
                vertexVal[t][1*numComp + c] = mab + mbc - mca;   // b
                vertexVal[t][2*numComp + c] = mbc + mca - mab;   // c
            }
        }

        double* out = writeSample(result, e);
        const double gx[4] = { g, 1.-g, g, 1.-g };
        const double gy[4] = { g, g, 1.-g, 1.-g };
        for (int q = 0; q < 4; ++q) {
            const double px = x0 + gx[q]*h, py = y0 + gy[q]*h;
            int best = 0;
            double bestScore = -1e30;
            double bc[3];
            for (size_t t = 0; t < tris.size(); ++t) {
                barycentric(tris[t], px, py, bc);
                const double score = std::min(bc[0], std::min(bc[1], bc[2]));
                if (score > bestScore) { bestScore = score; best = (int) t; }
            }
            barycentric(tris[best], px, py, bc);
            for (int c = 0; c < numComp; ++c)
                out[q*numComp + c] = bc[0]*vertexVal[best][0*numComp + c]
                                   + bc[1]*vertexVal[best][1*numComp + c]
                                   + bc[2]*vertexVal[best][2*numComp + c];
        }
    }
    return result;
}

} // namespace oxley
