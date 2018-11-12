
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "DudleyDomain.h"

#include <escript/index.h>

using escript::DataTypes::real_t;

namespace dudley {

escript::Domain_ptr DudleyDomain::create2D(dim_t NE0, dim_t NE1,
                                           real_t l0, real_t l1,
                                           bool optimize,
                                           escript::JMPI mpiInfo)
{
    const int DIM = 2;
    const int LEFTTAG = 1;    // boundary x1=0
    const int RIGHTTAG = 2;   // boundary x1=1
    const int BOTTOMTAG = 10; // boundary x2=0
    const int TOPTAG = 20;    // boundary x2=1

#ifdef Dudley_TRACE
    double time0 = Dudley_timer();
#endif

    const int myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    NE0 = std::max((dim_t)1, NE0);
    NE1 = std::max((dim_t)1, NE1);
    const dim_t N0 = NE0 + 1;
    const dim_t N1 = NE1 + 1;

    // This code was originally copied from Finley's Rec4 constructor.
    // NE? refers to the number of rectangular elements in each direction.
    // The number of nodes produced is the same but the number of non-face
    // elements will double since each "rectangle" is split into two triangles.

    // allocate mesh
    std::stringstream name;
    name << "Triangular " << N0 << " x " << N1 << " (x 2) mesh";
    DudleyDomain* out = new DudleyDomain(name.str(), DIM, mpiInfo);

    ElementFile* elements = new ElementFile(Dudley_Tri3, mpiInfo);
    out->setElements(elements);
    ElementFile* faces = new ElementFile(Dudley_Line2, mpiInfo);
    out->setFaceElements(faces);
    ElementFile* points = new ElementFile(Dudley_Point1, mpiInfo);
    out->setPoints(points);

    const dim_t Nstride0 = 1;
    const dim_t Nstride1 = N0;
    dim_t local_NE0, local_NE1;
    index_t e_offset0 = 0, e_offset1 = 0;
    if (N1 == std::max(N0, N1)) {
        local_NE0 = NE0;
        e_offset0 = 0;
        mpiInfo->split(NE1, &local_NE1, &e_offset1);
    } else {
        mpiInfo->split(NE0, &local_NE0, &e_offset0);
        local_NE1 = NE1;
        e_offset1 = 0;
    }
    const index_t offset0 = e_offset0;
    const index_t offset1 = e_offset1;
    const dim_t local_N0 = local_NE0 > 0 ? local_NE0 + 1 : 0;
    const dim_t local_N1 = local_NE1 > 0 ? local_NE1 + 1 : 0;

    // get the number of surface elements
    dim_t NFaceElements = 0;
    dim_t NDOF0, NDOF1;
    if (local_NE0 > 0) {
        NDOF0 = N0;
        if (e_offset0 == 0)
            NFaceElements += local_NE1;
        if (local_NE0 + e_offset0 == NE0)
            NFaceElements += local_NE1;
    } else {
        NDOF0 = N0 - 1;
    }
    if (local_NE1 > 0) {
        NDOF1 = N1;
        if (e_offset1 == 0)
            NFaceElements += local_NE0;
        if (local_NE1 + e_offset1 == NE1)
            NFaceElements += local_NE0;
    } else {
        NDOF1 = N1 - 1;
    }

    NodeFile* nodes = out->getNodes();
    nodes->allocTable(local_N0 * local_N1);
    elements->allocTable(local_NE0 * local_NE1 * 2);
    faces->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i1 = 0; i1 < local_N1; i1++) {
        for (index_t i0 = 0; i0 < local_N0; i0++) {
            const dim_t k = i0 + local_N0 * i1;
            const dim_t global_i0 = i0 + offset0;
            const dim_t global_i1 = i1 + offset1;
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
    const index_t global_adjustment = (offset0 + offset1) % 2;

#pragma omp parallel for
    for (index_t i1 = 0; i1 < local_NE1; i1++) {
        for (index_t i0 = 0; i0 < local_NE0; i0++) {
            // we will split this "rectangle" into two triangles
            const dim_t k = 2 * (i0 + local_NE0 * i1);
            const index_t node0 = Nstride0 * (i0 + e_offset0)
                                + Nstride1 * (i1 + e_offset1);

            elements->Id[k] = 2 * (i0 + e_offset0) + NE0*(i1 + e_offset1);
            elements->Tag[k] = 0;
            elements->Owner[k] = myRank;
            elements->Id[k + 1] = elements->Id[k] + 1;
            elements->Tag[k + 1] = 0;
            elements->Owner[k + 1] = myRank;

            // a,b,c,d gives the nodes of the rectangle in clockwise order
            const index_t a = node0;
            const index_t b = node0 + Nstride0;
            const index_t c = node0 + Nstride1 + Nstride0;
            const index_t d = node0 + Nstride1;
            // For a little bit of variety
            if ((global_adjustment + node0) % 2) {
                elements->Nodes[INDEX2(0, k, NN)] = a;
                elements->Nodes[INDEX2(1, k, NN)] = b;
                elements->Nodes[INDEX2(2, k, NN)] = d;
                elements->Nodes[INDEX2(0, k + 1, NN)] = b;
                elements->Nodes[INDEX2(1, k + 1, NN)] = c;
                elements->Nodes[INDEX2(2, k + 1, NN)] = d;
            } else {
                elements->Nodes[INDEX2(0, k, NN)] = a;
                elements->Nodes[INDEX2(1, k, NN)] = b;
                elements->Nodes[INDEX2(2, k, NN)] = c;
                elements->Nodes[INDEX2(0, k + 1, NN)] = a;
                elements->Nodes[INDEX2(1, k + 1, NN)] = c;
                elements->Nodes[INDEX2(2, k + 1, NN)] = d;
            }
        }
    }

    // face elements
    NN = faces->numNodes;
    dim_t totalNECount = 2 * NE0 * NE1; // because we have split the rectangles
    dim_t faceNECount = 0;
    if (local_NE0 > 0) {
        // ** elements on boundary 001 (x1=0)
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                const dim_t k = i1 + faceNECount;
                const index_t node0 = Nstride1 * (i1 + e_offset1);
                faces->Id[k] = i1 + e_offset1 + totalNECount;
                faces->Tag[k] = LEFTTAG;
                faces->Owner[k] = myRank;
                faces->Nodes[INDEX2(0, k, NN)] = node0 + Nstride1;
                faces->Nodes[INDEX2(1, k, NN)] = node0;
            }
            faceNECount += local_NE1;
        }
        totalNECount += NE1;
        // ** elements on boundary 002 (x1=1)
        if (local_NE0 + e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                const dim_t k = i1 + faceNECount;
                const index_t node0 = Nstride0 * (NE0 - 1)
                                    + Nstride1 * (i1 + e_offset1);

                faces->Id[k] = (i1 + e_offset1) + totalNECount;
                faces->Tag[k] = RIGHTTAG;
                faces->Owner[k] = myRank;
                faces->Nodes[INDEX2(0, k, NN)] = node0 + Nstride0;
                faces->Nodes[INDEX2(1, k, NN)] = node0 + Nstride1 + Nstride0;
            }
            faceNECount += local_NE1;
        }
        totalNECount += NE1;
    }
    if (local_NE1 > 0) {
        // ** elements on boundary 010 (x2=0)
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const dim_t k = i0 + faceNECount;
                const index_t node0 = Nstride0 * (i0 + e_offset0);
                faces->Id[k] = e_offset0 + i0 + totalNECount;
                faces->Tag[k] = BOTTOMTAG;
                faces->Owner[k] = myRank;
                faces->Nodes[INDEX2(0, k, NN)] = node0;
                faces->Nodes[INDEX2(1, k, NN)] = node0 + Nstride0;
            }
            faceNECount += local_NE0;
        }
        totalNECount += NE0;
        // ** elements on boundary 020 (x2=1)
        if (local_NE1 + e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const dim_t k = i0 + faceNECount;
                const index_t node0 = Nstride0 * (i0 + e_offset0)
                                    + Nstride1 * (NE1 - 1);

                faces->Id[k] = i0 + e_offset0 + totalNECount;
                faces->Tag[k] = TOPTAG;
                faces->Owner[k] = myRank;
                faces->Nodes[INDEX2(0, k, NN)] = node0 + Nstride1 + Nstride0;
                faces->Nodes[INDEX2(1, k, NN)] = node0 + Nstride1;
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
    return escript::Domain_ptr(out);
}

} // namespace dudley

