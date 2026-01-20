
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "FinleyDomain.h"

#include <escript/index.h>

using escript::DataTypes::real_t;

namespace finley {

escript::Domain_ptr FinleyDomain::createRec8(dim_t NE0, dim_t NE1,
                                             real_t l0, real_t l1,
                                             bool periodic0, bool periodic1,
                                             int order, int reducedOrder,
                                             bool useElementsOnFace,
                                             bool useFullElementOrder,
                                             bool useMacroElements,
                                             bool optimize,
                                             escript::JMPI mpiInfo)
{
    const int N_PER_E = 2;
    const int DIM = 2;
    const int LEFTTAG = 1;    // boundary x1=0
    const int RIGHTTAG = 2;   // boundary x1=1
    const int BOTTOMTAG = 10; // boundary x2=0
    const int TOPTAG = 20;    // boundary x2=1
    dim_t Nstride0 = 0, Nstride1 = 0, local_NE0, local_NE1;
    index_t e_offset0 = 0, e_offset1 = 0;
    const bool generateAllNodes = useFullElementOrder || useMacroElements;

    const int myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    NE0 = std::max((dim_t)1, NE0);
    NE1 = std::max((dim_t)1, NE1);
    const dim_t N0 = N_PER_E*NE0 + 1;
    const dim_t N1 = N_PER_E*NE1 + 1;

    // allocate mesh
    std::stringstream name;
    name << "Rectangular " << N0 << " x " << N1 << " mesh";
    FinleyDomain* out = new FinleyDomain(name.str(), DIM, mpiInfo);

    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    if (generateAllNodes) {
        if (useMacroElements) {
            refElements.reset(new ReferenceElementSet(Rec9Macro, order, reducedOrder));
        } else {
            refElements.reset(new ReferenceElementSet(Rec9, order, reducedOrder));
        }
        if (useElementsOnFace) {
            throw escript::NotImplementedError("rich elements for Rec9 elements are not supported yet.");
        } else {
            if (useMacroElements) {
                refFaceElements.reset(new ReferenceElementSet(Line3Macro, order, reducedOrder));
            } else {
                refFaceElements.reset(new ReferenceElementSet(Line3, order, reducedOrder));
            }
            refContactElements.reset(new ReferenceElementSet(Line3_Contact, order, reducedOrder));
        }
    } else { // !generateAllNodes
        refElements.reset(new ReferenceElementSet(Rec8, order, reducedOrder));
        if (useElementsOnFace) {
            refFaceElements.reset(new ReferenceElementSet(Rec8Face, order, reducedOrder));
            refContactElements.reset(new ReferenceElementSet(Rec8Face_Contact, order, reducedOrder));
        } else {
            refFaceElements.reset(new ReferenceElementSet(Line3, order, reducedOrder));
            refContactElements.reset(new ReferenceElementSet(Line3_Contact, order, reducedOrder));
        }
    }
    refPoints.reset(new ReferenceElementSet(Point1, order, reducedOrder));
    if (!refPoints->referenceElement)
        throw FinleyException("ERRRRORRRRR!!");

    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    out->setElements(elements);
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    out->setFaceElements(faces);
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));
    out->setPoints(new ElementFile(refPoints, mpiInfo));

    // work out the largest dimension
    if (N1 == std::max(N0, N1)) {
        Nstride0 = 1;
        Nstride1 = N0;
        local_NE0 = NE0;
        e_offset0 = 0;
        mpiInfo->split(NE1, &local_NE1, &e_offset1);
    } else {
        Nstride0 = N1;
        Nstride1 = 1;
        mpiInfo->split(NE0, &local_NE0, &e_offset0);
        local_NE1 = NE1;
        e_offset1 = 0;
    }
    const index_t offset0 = e_offset0 * N_PER_E;
    const index_t offset1 = e_offset1 * N_PER_E;
    const dim_t local_N0 = local_NE0 > 0 ? local_NE0*N_PER_E+1 : 0;
    const dim_t local_N1 = local_NE1 > 0 ? local_NE1*N_PER_E+1 : 0;
    dim_t NDOF0 = 0, NDOF1 = 0;

    // get the number of surface elements
    dim_t NFaceElements = 0;
    if (!periodic0 && local_NE0 > 0) {
        NDOF0 = N0;
        if (e_offset0 == 0)
            NFaceElements += local_NE1;
        if (local_NE0 + e_offset0 == NE0)
            NFaceElements += local_NE1;
    } else {
        NDOF0 = N0 - 1;
    }

    if (!periodic1 && local_NE1 > 0) {
        NDOF1 = N1;
        if (e_offset1 == 0)
            NFaceElements += local_NE0;
        if (local_NE1 + e_offset1 == NE1)
            NFaceElements += local_NE0;
    } else {
        NDOF1 = N1 - 1;
    }

    // allocate tables
    NodeFile* nodes = out->getNodes();
    nodes->allocTable(local_N0 * local_N1);
    elements->allocTable(local_NE0 * local_NE1);
    faces->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i1 = 0; i1 < local_N1; i1++) {
        for (index_t i0 = 0; i0 < local_N0; i0++) {
            const dim_t k = i0 + local_N0 * i1;
            const index_t global_i0 = i0 + offset0;
            const index_t global_i1 = i1 + offset1;
            nodes->Coordinates[INDEX2(0, k, DIM)] = (real_t)global_i0 / (real_t)(N0 - 1) * l0;
            nodes->Coordinates[INDEX2(1, k, DIM)] = (real_t)global_i1 / (real_t)(N1 - 1) * l1;
            nodes->Id[k] = Nstride0 * global_i0 + Nstride1 * global_i1;
            nodes->Tag[k] = 0;
            nodes->globalDegreesOfFreedom[k] = Nstride0 * (global_i0 % NDOF0)
                                             + Nstride1 * (global_i1 % NDOF1);
        }
    }

    // set the elements
    dim_t NN = elements->numNodes;
#pragma omp parallel for
    for (index_t i1 = 0; i1 < local_NE1; i1++) {
        for (index_t i0 = 0; i0 < local_NE0; i0++) {
            const dim_t k = i0 + local_NE0 * i1;
            const index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0)
                                + Nstride1 * N_PER_E * (i1 + e_offset1);

            elements->Id[k] = (i0 + e_offset0) + NE0*(i1 + e_offset1);
            elements->Tag[k] = 0;
            elements->Owner[k] = myRank;

            elements->Nodes[INDEX2(0, k, NN)] = node0;
            elements->Nodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride0;
            elements->Nodes[INDEX2(2, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
            elements->Nodes[INDEX2(3, k, NN)] = node0 + 2 * Nstride1;
            elements->Nodes[INDEX2(4, k, NN)] = node0 + 1 * Nstride0;
            elements->Nodes[INDEX2(5, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
            elements->Nodes[INDEX2(6, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
            elements->Nodes[INDEX2(7, k, NN)] = node0 + Nstride1;
            if (generateAllNodes) {
                elements->Nodes[INDEX2(8, k, NN)] = node0 + Nstride1 + Nstride0;
            }
        }
    }

    // face elements
    NN = faces->numNodes;
    dim_t totalNECount = NE0 * NE1;
    dim_t faceNECount = 0;
    index_t* eNodes = faces->Nodes;

    if (!periodic0 && local_NE0 > 0) {
        // ** elements on boundary 001 (x1=0)
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                const dim_t k = i1 + faceNECount;
                const index_t node0 = Nstride1 * N_PER_E * (i1 + e_offset1);

                faces->Id[k] = i1 + e_offset1 + totalNECount;
                faces->Tag[k] = LEFTTAG;
                faces->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(1, k, NN)] = node0;
                    eNodes[INDEX2(2, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(3, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(4, k, NN)] = node0 + Nstride1;
                    eNodes[INDEX2(5, k, NN)] = node0 + Nstride0;
                    eNodes[INDEX2(6, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(7, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
                } else {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(1, k, NN)] = node0;
                    eNodes[INDEX2(2, k, NN)] = node0 + Nstride1;
                }
            }
            faceNECount += local_NE1;
        }
        totalNECount += NE1;
        // ** elements on boundary 002 (x1=1)
        if (local_NE0 + e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                const dim_t k = i1 + faceNECount;
                const index_t node0 = Nstride0 * N_PER_E * (NE0 - 1)
                                    + Nstride1 * N_PER_E * (i1 + e_offset1);

                faces->Id[k] = (i1 + e_offset1) + totalNECount;
                faces->Tag[k] = RIGHTTAG;
                faces->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(2, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(3, k, NN)] = node0;
                    eNodes[INDEX2(4, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(5, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
                    eNodes[INDEX2(6, k, NN)] = node0 + Nstride1;
                    eNodes[INDEX2(7, k, NN)] = node0 + Nstride0;
                } else {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(2, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
                }
            }
            faceNECount += local_NE1;
        }
        totalNECount += NE1;
    }

    if (!periodic1 && local_NE1 > 0) {
        // ** elements on boundary 010 (x2=0)
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const dim_t k = i0 + faceNECount;
                const index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0);
                faces->Id[k] = e_offset0 + i0 + totalNECount;
                faces->Tag[k] = BOTTOMTAG;
                faces->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0, k, NN)] = node0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(2, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(3, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(4, k, NN)] = node0 + Nstride0;
                    eNodes[INDEX2(5, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(6, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
                    eNodes[INDEX2(7, k, NN)] = node0 + Nstride1;
                } else {
                    eNodes[INDEX2(0, k, NN)] = node0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(2, k, NN)] = node0 + Nstride0;
                }
            }
            faceNECount += local_NE0;
        }
        totalNECount += NE0;
        // ** elements on boundary 020 (x2=1)
        if (local_NE1 + e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const dim_t k = i0 + faceNECount;
                const index_t node0 = Nstride0 * N_PER_E * (i0 + e_offset0)
                                    + Nstride1 * N_PER_E * (NE1 - 1);

                faces->Id[k] = i0 + e_offset0 + totalNECount;
                faces->Tag[k] = TOPTAG;
                faces->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(2, k, NN)] = node0;
                    eNodes[INDEX2(3, k, NN)] = node0 + 2 * Nstride0;
                    eNodes[INDEX2(4, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
                    eNodes[INDEX2(5, k, NN)] = node0 + Nstride1;
                    eNodes[INDEX2(6, k, NN)] = node0 + Nstride0;
                    eNodes[INDEX2(7, k, NN)] = node0 + Nstride1 + 2 * Nstride0;
                } else {
                    eNodes[INDEX2(0, k, NN)] = node0 + 2 * Nstride1 + 2 * Nstride0;
                    eNodes[INDEX2(1, k, NN)] = node0 + 2 * Nstride1;
                    eNodes[INDEX2(2, k, NN)] = node0 + 2 * Nstride1 + Nstride0;
                }
            }
            faceNECount += local_NE0;
        }
        totalNECount += NE0;
    }

    // add tag names
    out->setTagMap("top", TOPTAG);
    out->setTagMap("bottom", BOTTOMTAG);
    out->setTagMap("left", LEFTTAG);
    out->setTagMap("right", RIGHTTAG);

    // prepare mesh for further calculations
    out->resolveNodeIds();
    out->prepare(optimize);
    return out->getPtr();
}

} // namespace finley

