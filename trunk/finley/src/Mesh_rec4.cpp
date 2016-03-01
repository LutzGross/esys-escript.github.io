
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************

  Finley: generates rectangular meshes

  Generates a numElements[0] x numElements[1] mesh with first order elements
  (Rec4) in the rectangle [0,Length[0]] x [0,Length[1]].
  order is the desired accuracy of the integration scheme.

*****************************************************************************/

#include "RectangularMesh.h"

using escript::DataTypes::real_t;

namespace finley {

Mesh* RectangularMesh_Rec4(const dim_t* numElements, const double* Length,
                           const bool* periodic, int order, int reduced_order,
                           bool useElementsOnFace, bool useFullElementOrder,
                           bool optimize, escript::JMPI& mpiInfo)
{
    const int N_PER_E = 1;
    const int DIM = 2;
    dim_t Nstride0=0, Nstride1=0, local_NE0, local_NE1;
    index_t e_offset0=0, e_offset1=0;

    const int myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    const dim_t NE0 = std::max(dim_t(1),numElements[0]);
    const dim_t NE1 = std::max(dim_t(1),numElements[1]);
    const dim_t N0 = N_PER_E*NE0+1;
    const dim_t N1 = N_PER_E*NE1+1;

    // allocate mesh
    std::stringstream name;
    name << "Rectangular " << N0 << " x " << N1 << " mesh";
    Mesh* out = new Mesh(name.str(), DIM, mpiInfo);

    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    if (useElementsOnFace) {
        refFaceElements.reset(new ReferenceElementSet(Rec4Face, order, reduced_order));
        refContactElements.reset(new ReferenceElementSet(Rec4Face_Contact, order, reduced_order));
    } else {
        refFaceElements.reset(new ReferenceElementSet(Line2, order, reduced_order));
        refContactElements.reset(new ReferenceElementSet(Line2_Contact, order, reduced_order));
    }
    refElements.reset(new ReferenceElementSet(Rec4, order, reduced_order));
    refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));

    out->setPoints(new ElementFile(refPoints, mpiInfo));
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));
    out->setFaceElements(new ElementFile(refFaceElements, mpiInfo));
    out->setElements(new ElementFile(refElements, mpiInfo));

    // work out the largest dimension
    if (N1 == std::max(N0,N1)) {
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
    const index_t offset0 = e_offset0*N_PER_E;
    const index_t offset1 = e_offset1*N_PER_E;
    const dim_t local_N0 = local_NE0>0 ? local_NE0*N_PER_E+1 : 0;
    const dim_t local_N1 = local_NE1>0 ? local_NE1*N_PER_E+1 : 0;
    dim_t NDOF0=0, NDOF1=0;

    // get the number of surface elements
    dim_t NFaceElements = 0;
    if (!periodic[0] && local_NE0>0) {
        NDOF0=N0;
        if (e_offset0 == 0)
            NFaceElements+=local_NE1;
        if (local_NE0+e_offset0 == NE0)
            NFaceElements+=local_NE1;
    } else {
        NDOF0=N0-1;
    }

    if (!periodic[1] && local_NE1>0) {
        NDOF1=N1;
        if (e_offset1 == 0)
            NFaceElements+=local_NE0;
        if (local_NE1+e_offset1 == NE1)
            NFaceElements+=local_NE0;
    } else {
        NDOF1=N1-1;
    }

    // allocate tables
    out->Nodes->allocTable(local_N0*local_N1);
    out->Elements->allocTable(local_NE0*local_NE1);
    out->FaceElements->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i1=0; i1<local_N1; i1++) {
        for (index_t i0=0; i0<local_N0; i0++) {
            const dim_t k = i0+local_N0*i1;
            const index_t global_i0 = i0+offset0;
            const index_t global_i1 = i1+offset1;
            out->Nodes->Coordinates[INDEX2(0,k,DIM)]=(real_t)global_i0/(real_t)(N0-1)*Length[0];
            out->Nodes->Coordinates[INDEX2(1,k,DIM)]=(real_t)global_i1/(real_t)(N1-1)*Length[1];
            out->Nodes->Id[k] = Nstride0*global_i0 + Nstride1*global_i1;
            out->Nodes->Tag[k]=0;
            out->Nodes->globalDegreesOfFreedom[k] = Nstride0*(global_i0%NDOF0)
                                                  + Nstride1*(global_i1%NDOF1);
        }
    }

    // set the elements
    dim_t NN = out->Elements->numNodes;
#pragma omp parallel for
    for (index_t i1=0; i1<local_NE1; i1++) {
        for (index_t i0=0; i0<local_NE0; i0++) {
            const dim_t k = i0+local_NE0*i1;
            const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                + Nstride1*N_PER_E*(i1+e_offset1);

            out->Elements->Id[k] = (i0+e_offset0) + NE0*(i1+e_offset1);
            out->Elements->Tag[k] = 0;
            out->Elements->Owner[k] = myRank;

            out->Elements->Nodes[INDEX2(0,k,NN)]=node0;
            out->Elements->Nodes[INDEX2(1,k,NN)]=node0+Nstride0;
            out->Elements->Nodes[INDEX2(2,k,NN)]=node0+Nstride1+Nstride0;
            out->Elements->Nodes[INDEX2(3,k,NN)]=node0+Nstride1;
        }
    }

    // face elements
    NN=out->FaceElements->numNodes;
    dim_t totalNECount=NE0*NE1;
    dim_t faceNECount = 0;
    index_t* eNodes = out->FaceElements->Nodes;

    if (!periodic[0] && local_NE0>0) {
        // **  elements on boundary 001 (x1=0):
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i1=0; i1<local_NE1; i1++) {
                const dim_t k = i1+faceNECount;
                const index_t node0 = Nstride1*N_PER_E*(i1+e_offset1);

                out->FaceElements->Id[k] = i1+e_offset1+totalNECount;
                out->FaceElements->Tag[k] = 1;
                out->FaceElements->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride1;
                    eNodes[INDEX2(1,k,NN)]=node0;
                    eNodes[INDEX2(2,k,NN)]=node0+Nstride0;
                    eNodes[INDEX2(3,k,NN)]=node0+Nstride1+Nstride0;
                } else {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride1;
                    eNodes[INDEX2(1,k,NN)]=node0;
                }
            }
            faceNECount+=local_NE1;
        }
        totalNECount+=NE1;

        // **  elements on boundary 002 (x1=1):
        if (local_NE0+e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i1=0; i1<local_NE1; i1++) {
                const dim_t k = i1+faceNECount;
                const index_t node0 = Nstride0*N_PER_E*(NE0-1)
                                    + Nstride1*N_PER_E*(i1+e_offset1);

                out->FaceElements->Id[k] = (i1+e_offset1)+totalNECount;
                out->FaceElements->Tag[k] = 2;
                out->FaceElements->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride1+Nstride0;
                    eNodes[INDEX2(2,k,NN)]=node0+Nstride1;
                    eNodes[INDEX2(3,k,NN)]=node0;
                } else {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride1+Nstride0;
                }
            }
            faceNECount+=local_NE1;
        }
        totalNECount+=NE1;
    }

    if (!periodic[1] && local_NE1>0) {
        // **  elements on boundary 010 (x2=0):
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i0=0; i0<local_NE0; i0++) {
                const dim_t k = i0+faceNECount;
                const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0);
                out->FaceElements->Id[k] = e_offset0+i0+totalNECount;
                out->FaceElements->Tag[k] = 10;
                out->FaceElements->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0,k,NN)]=node0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride0;
                    eNodes[INDEX2(2,k,NN)]=node0+Nstride1+Nstride0;
                    eNodes[INDEX2(3,k,NN)]=node0+Nstride1;
                } else {
                    eNodes[INDEX2(0,k,NN)]=node0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride0;
                }
            }
            faceNECount+=local_NE0;
        }
        totalNECount+=NE0;

        // **  elements on boundary 020 (x2=1):
        if (local_NE1+e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i0=0; i0<local_NE0; i0++) {
                const dim_t k = i0+faceNECount;
                const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                    + Nstride1*N_PER_E*(NE1-1);

                out->FaceElements->Id[k] = i0+e_offset0+totalNECount;
                out->FaceElements->Tag[k] = 20;
                out->FaceElements->Owner[k] = myRank;
                if (useElementsOnFace) {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride1+Nstride0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride1;
                    eNodes[INDEX2(2,k,NN)]=node0;
                    eNodes[INDEX2(3,k,NN)]=node0+Nstride0;
                } else {
                    eNodes[INDEX2(0,k,NN)]=node0+Nstride1+Nstride0;
                    eNodes[INDEX2(1,k,NN)]=node0+Nstride1;
                }
            }
            faceNECount+=local_NE0;
        }
        totalNECount+=NE0;
    }

    // add tag names
    out->addTagMap("top", 20);
    out->addTagMap("bottom", 10);
    out->addTagMap("left", 1);
    out->addTagMap("right", 2);

    // prepare mesh for further calculations
    out->resolveNodeIds();
    out->prepare(optimize);
    return out;
}

} // namespace finley

