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
#include <oxley/OxleyException.h>

#include <finley/FinleyDomain.h>

#include <algorithm>
#include <array>
#include <map>
#include <sstream>

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
    if (!dom.isConforming())
        throw OxleyException("toFinley: the forest has hanging nodes. Only "
                "conforming forests can be converted so far - the simplex split "
                "does not yet number the hanging positions.");

    const MeshAccess m = dom.getMeshAccess();
    if (m.numDim != 2 && m.numDim != 3) {
        std::stringstream ss;
        ss << "toFinley: unsupported dimension " << m.numDim;
        throw OxleyException(ss.str());
    }

    const int V = m.nodesPerElement;
    const int FV = m.nodesPerFace;
    // global id of a local node, the key the split rule is built on
    const std::vector<long>& gidOf = m.nodeGlobalId;

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
        // ---- 2D: each quad fans into triangles from its lowest-id corner ---
        out.elementType = finley::Tri3;
        out.faceElementType = finley::Line2;

        elems.reserve((size_t) m.numElements * 2 * 3);
        out.elementTag.reserve(m.numElements * 2);
        for (long e = 0; e < m.numElements; ++e) {
            const long* en = &m.elementNodes[(size_t) e * V];
            const int mpos = lowestIdPos(rectPolygon, 4,
                    [&](int c) { return gidOf[en[c]]; });
            for (int t = 0; t < 2; ++t) {
                long a = en[rectPolygon[mpos]];
                long b = en[rectPolygon[(mpos + 1 + t) % 4]];
                long c = en[rectPolygon[(mpos + 2 + t) % 4]];
                if (signedArea(m.nodeCoords, a, b, c) < 0.)
                    std::swap(b, c);
                elems.push_back((index_t) a);
                elems.push_back((index_t) b);
                elems.push_back((index_t) c);
                out.elementTag.push_back((int) m.elementTags[e]);
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

} // namespace oxley
