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

#include <finley/FinleyDomain.h>

#include <algorithm>
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
    // nodeGlobalId cannot do that, since it numbers materialised nodes in each
    // rank's own creation order. Brick does not build it yet, hence the fallback.
    const std::vector<long>& gidOf =
            m.nodeExportId.empty() ? m.nodeGlobalId : m.nodeExportId;

    finley::MeshArrays out;
    out.numDim = m.numDim;

    // nodes: hand over every local node, owned or ghost. finley tolerates a
    // node being supplied by more than one rank and fetches any it still needs.
    out.nodeId.resize(m.numNodes);
    for (long i = 0; i < m.numNodes; ++i)
        out.nodeId[i] = (index_t) gidOf[i];
    out.nodeCoords = m.nodeCoords;

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
    index_t faceIdOffset = globalNumElements;
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

    out.tagMap["left"] = 1;
    out.tagMap["right"] = 2;
    out.tagMap["bottom"] = 10;
    out.tagMap["top"] = 20;
    if (m.numDim == 3) {
        out.tagMap["front"] = 100;
        out.tagMap["back"] = 200;
    }

    std::stringstream name;
    name << "finley mesh from oxley " << (m.numDim == 2 ? "Rectangle" : "Brick");

    return finley::FinleyDomain::createFromArrays(out, name.str(), order,
                                                  reducedOrder, optimize,
                                                  mpiInfo);
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

} // namespace oxley
