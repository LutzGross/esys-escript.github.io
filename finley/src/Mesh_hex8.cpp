
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


/****************************************************************************

  Finley: generates rectangular meshes

  Generates a numElements[0] x numElements[1] x numElements[2] mesh with
  first order elements (Hex8) in the brick
  [0,Length[0]] x [0,Length[1]] x [0,Length[2]].
  order is the desired accuracy of the integration scheme.

*****************************************************************************/

#include "FinleyDomain.h"

#include <escript/index.h>

#define MAX3(_arg1_,_arg2_,_arg3_) std::max(_arg1_,std::max(_arg2_,_arg3_))

using escript::DataTypes::real_t;

namespace finley {

escript::Domain_ptr FinleyDomain::createHex8(dim_t NE0, dim_t NE1, dim_t NE2,
                            double l0, double l1, double l2,
                            bool periodic0, bool periodic1, bool periodic2,
                            int order, int reduced_order,
                            bool useElementsOnFace,
                            bool optimize, escript::JMPI mpiInfo)
{
    const int N_PER_E = 1;
    const int DIM = 3;
    dim_t Nstride0=0, Nstride1=0, Nstride2=0, local_NE0, local_NE1, local_NE2;
    index_t e_offset0, e_offset1, e_offset2;

    const int myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    NE0 = std::max(dim_t(1), NE0);
    NE1 = std::max(dim_t(1), NE1);
    NE2 = std::max(dim_t(1), NE2);
    const dim_t N0 = N_PER_E*NE0+1;
    const dim_t N1 = N_PER_E*NE1+1;
    const dim_t N2 = N_PER_E*NE2+1;

    // allocate mesh
    std::stringstream name;
    name << "Rectangular " << N0 << " x " << N1 << " x " << N2 << " mesh";
    FinleyDomain* out = new FinleyDomain(name.str(), DIM, mpiInfo);

    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    if (useElementsOnFace) {
        refFaceElements.reset(new ReferenceElementSet(Hex8Face, order, reduced_order));
        refContactElements.reset(new ReferenceElementSet(Hex8Face_Contact, order, reduced_order));
    } else {
        refFaceElements.reset(new ReferenceElementSet(Rec4, order, reduced_order));
        refContactElements.reset(new ReferenceElementSet(Rec4_Contact, order, reduced_order));
    }
    refElements.reset(new ReferenceElementSet(Hex8, order, reduced_order));
    refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));

    out->setPoints(new ElementFile(refPoints, mpiInfo));
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    out->setFaceElements(faces);
    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    out->setElements(elements);

    // work out the largest dimension
    if (N2==MAX3(N0,N1,N2)) {
        Nstride0 = 1;
        Nstride1 = N0;
        Nstride2 = N0*N1;
        local_NE0 = NE0;
        e_offset0 = 0;
        local_NE1 = NE1;
        e_offset1 = 0;
        mpiInfo->split(NE2, &local_NE2, &e_offset2);
    } else if (N1==MAX3(N0,N1,N2)) {
        Nstride0 = N2;
        Nstride1 = N0*N2;
        Nstride2 = 1;
        local_NE0 = NE0;
        e_offset0 = 0;
        mpiInfo->split(NE1, &local_NE1, &e_offset1);
        local_NE2 = NE2;
        e_offset2 = 0;
    } else {
        Nstride0 = N1*N2;
        Nstride1 = 1;
        Nstride2 = N1;
        mpiInfo->split(NE0, &local_NE0, &e_offset0);
        local_NE1 = NE1;
        e_offset1 = 0;
        local_NE2 = NE2;
        e_offset2 = 0;
    }
    const index_t offset0 = e_offset0*N_PER_E;
    const index_t offset1 = e_offset1*N_PER_E;
    const index_t offset2 = e_offset2*N_PER_E;
    const dim_t local_N0 = local_NE0>0 ? local_NE0*N_PER_E+1 : 0;
    const dim_t local_N1 = local_NE1>0 ? local_NE1*N_PER_E+1 : 0;
    const dim_t local_N2 = local_NE2>0 ? local_NE2*N_PER_E+1 : 0;
    dim_t NDOF0=0, NDOF1=0, NDOF2=0;

    // get the number of surface elements
    dim_t NFaceElements = 0;
    if (!periodic2 && local_NE2 > 0) {
        NDOF2=N2;
        if (offset2==0)
            NFaceElements+=local_NE1*local_NE0;
        if (local_NE2+e_offset2 == NE2)
            NFaceElements+=local_NE1*local_NE0;
    } else {
        NDOF2=N2-1;
    }

    if (!periodic0 && local_NE0 > 0) {
        NDOF0=N0;
        if (e_offset0 == 0)
            NFaceElements+=local_NE1*local_NE2;
        if (local_NE0+e_offset0 == NE0)
            NFaceElements+=local_NE1*local_NE2;
    } else {
        NDOF0=N0-1;
    }
    if (!periodic1 && local_NE1 > 0) {
        NDOF1=N1;
        if (e_offset1 == 0)
            NFaceElements+=local_NE0*local_NE2;
        if (local_NE1+e_offset1 == NE1)
            NFaceElements+=local_NE0*local_NE2;
    } else {
        NDOF1=N1-1;
    }

    // allocate tables
    NodeFile* nodes = out->getNodes();
    nodes->allocTable(local_N0*local_N1*local_N2);
    elements->allocTable(local_NE0*local_NE1*local_NE2);
    faces->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i2=0; i2<local_N2; i2++) {
        for (index_t i1=0; i1<local_N1; i1++) {
            for (index_t i0=0; i0<local_N0; i0++) {
                const dim_t k = i0+local_N0*i1+local_N0*local_N1*i2;
                const index_t global_i0 = i0+offset0;
                const index_t global_i1 = i1+offset1;
                const index_t global_i2 = i2+offset2;
                nodes->Coordinates[INDEX2(0,k,DIM)]=(real_t)global_i0/(real_t)(N0-1)*l0;
                nodes->Coordinates[INDEX2(1,k,DIM)]=(real_t)global_i1/(real_t)(N1-1)*l1;
                nodes->Coordinates[INDEX2(2,k,DIM)]=(real_t)global_i2/(real_t)(N2-1)*l2;
                nodes->Id[k] = Nstride0*global_i0+Nstride1*global_i1+Nstride2*global_i2;
                nodes->Tag[k] = 0;
                nodes->globalDegreesOfFreedom[k] = Nstride0*(global_i0%NDOF0)
                                                 + Nstride1*(global_i1%NDOF1)
                                                 + Nstride2*(global_i2%NDOF2);
            }
        }
    }

    // set the elements
    dim_t NN = elements->numNodes;
    index_t* eNodes = elements->Nodes;
#pragma omp parallel for
    for (index_t i2 = 0; i2 < local_NE2; i2++) {
        for (index_t i1 = 0; i1 < local_NE1; i1++) {
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const dim_t k = i0+local_NE0*i1+local_NE0*local_NE1*i2;
                const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                    + Nstride1*N_PER_E*(i1+e_offset1)
                                    + Nstride2*N_PER_E*(i2+e_offset2);

                elements->Id[k] = (i0+e_offset0) + NE0*(i1+e_offset1)
                                                 + NE0*NE1*(i2+e_offset2);
                elements->Tag[k] = 0;
                elements->Owner[k] = myRank;

                eNodes[INDEX2(0,k,NN)] = node0;
                eNodes[INDEX2(1,k,NN)] = node0+Nstride0;
                eNodes[INDEX2(2,k,NN)] = node0+Nstride1+Nstride0;
                eNodes[INDEX2(3,k,NN)] = node0+Nstride1;
                eNodes[INDEX2(4,k,NN)] = node0+Nstride2;
                eNodes[INDEX2(5,k,NN)] = node0+Nstride2+Nstride0;
                eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                eNodes[INDEX2(7,k,NN)] = node0+Nstride2+Nstride1;
            }
        }
    }

    // face elements
    NN = faces->numNodes;
    dim_t totalNECount = NE0*NE1*NE2;
    dim_t faceNECount = 0;
    eNodes = faces->Nodes;

    // these are the quadrilateral elements on boundary 1 (x3=0):
    if (!periodic2 && local_NE2 > 0) {
        // **  elements on boundary 100 (x3=0):
        if (e_offset2 == 0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                for (index_t i0=0; i0 < local_NE0; i0++) {
                    const dim_t k = i0 + local_NE0*i1 + faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(i1+e_offset1);

                    faces->Id[k] = (i0+e_offset0) + NE0*(i1+e_offset1)
                                                  + totalNECount;
                    faces->Tag[k] = 100;
                    faces->Owner[k] = myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(4,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(7,k,NN)] = node0+Nstride2+Nstride0;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride0;
                    }
                }
            }
            faceNECount += local_NE1*local_NE0;
        }
        totalNECount += NE1*NE0;

        // **  elements on boundary 200 (x3=1):
        if (local_NE2+e_offset2 == NE2) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const dim_t k = i0 + local_NE0*i1 + faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(NE2-1);

                    faces->Id[k] = (i0+e_offset0) + NE0*(i1+e_offset1)
                                                  + totalNECount;
                    faces->Tag[k] = 200;
                    faces->Owner[k] = myRank;
                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2+         Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2+Nstride1;

                        eNodes[INDEX2(4,k,NN)] = node0;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(6,k,NN)] = node0+         Nstride1+Nstride0;
                        eNodes[INDEX2(7,k,NN)] = node0+         Nstride1;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2         +Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2+Nstride1;
                    }
                }
            }
            faceNECount += local_NE1*local_NE0;
        }
        totalNECount += NE1*NE0;
    } // !periodic2 && local_NE2 > 0

    if (!periodic0 && local_NE0 > 0) {
        // **  elements on boundary 001 (x1=0):
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i1 = 0; i1 < local_NE1; i1++) {
                    const dim_t k = i1 + local_NE1*i2 + faceNECount;
                    const index_t node0 = Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);
                    faces->Id[k] = (i1+e_offset1) + NE1*(i2+e_offset2)
                                                  + totalNECount;
                    faces->Tag[k] = 1;
                    faces->Owner[k] = myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(4,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride2+Nstride0;
                        eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(7,k,NN)] = node0+Nstride1+Nstride0;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride1;
                    }
                }
            }
            faceNECount += local_NE1*local_NE2;
        }
        totalNECount += NE1*NE2;

        // **  elements on boundary 002 (x1=1):
        if (local_NE0+e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i1 = 0; i1 < local_NE1; i1++) {
                    const dim_t k = i1 + local_NE1*i2 + faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(NE0-1)
                                        + Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);
                    faces->Id[k] = (i1+e_offset1) + NE1*(i2+e_offset2)
                                                  + totalNECount;
                    faces->Tag[k] = 2;
                    faces->Owner[k] = myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride1+Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2+Nstride0;

                        eNodes[INDEX2(4,k,NN)] = node0;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(7,k,NN)] = node0+Nstride2;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride1+Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2+Nstride0;
                    }
                }
            }
            faceNECount += local_NE1*local_NE2;
        }
        totalNECount += NE1*NE2;
    } // !periodic0 && local_NE0 > 0

    if (!periodic1 && local_NE1 > 0) {
        // **  elements on boundary 010 (x2=0):
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const dim_t k = i0 + local_NE0*i2 + faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride2*N_PER_E*(i2+e_offset2);

                    faces->Id[k] = (i2+e_offset2) + NE2*(e_offset0+i0)
                                                  + totalNECount;
                    faces->Tag[k] = 10;
                    faces->Owner[k] = myRank;
                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2;

                        eNodes[INDEX2(4,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride1+Nstride0;
                        eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(7,k,NN)] = node0+Nstride2+Nstride1;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride0;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride2;
                    }
                }
            }
            faceNECount += local_NE0*local_NE2;
        }
        totalNECount += NE0*NE2;

        // **  elements on boundary 020 (x2=1):
        if (local_NE1+e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const dim_t k = i0 + local_NE0*i2 + faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(NE1-1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);

                    faces->Id[k] = (i2+e_offset2) + NE2*(i0+e_offset0)
                                                  + totalNECount;
                    faces->Tag[k] = 20;
                    faces->Owner[k] = myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride1+Nstride0;

                        eNodes[INDEX2(4,k,NN)] = node0;
                        eNodes[INDEX2(5,k,NN)] = node0+Nstride2;
                        eNodes[INDEX2(6,k,NN)] = node0+Nstride2+Nstride0;
                        eNodes[INDEX2(7,k,NN)] = node0+Nstride0;
                    } else {
                        eNodes[INDEX2(0,k,NN)] = node0+Nstride1;
                        eNodes[INDEX2(1,k,NN)] = node0+Nstride2+Nstride1;
                        eNodes[INDEX2(2,k,NN)] = node0+Nstride2+Nstride1+Nstride0;
                        eNodes[INDEX2(3,k,NN)] = node0+Nstride1+Nstride0;
                    }
                }
            }
            faceNECount += local_NE0*local_NE2;
        }
        totalNECount += NE0*NE2;
    }

    // add tag names
    out->setTagMap("top", 200);
    out->setTagMap("bottom", 100);
    out->setTagMap("left", 1);
    out->setTagMap("right", 2);
    out->setTagMap("front", 10);
    out->setTagMap("back", 20);

    // prepare mesh for further calculations
    out->resolveNodeIds();
    out->prepare(optimize);
    return out->getPtr();
}

} // namespace finley

