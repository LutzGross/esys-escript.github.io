
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

  Generates a numElements[0] x numElements[1] x numElements[2] mesh with
  second order elements (Hex20) in the brick
  [0,Length[0]] x [0,Length[1]] x [0,Length[2]].
  order is the desired accuracy of the integration scheme.

*****************************************************************************/

#define ESNEEDPYTHON
#include <esysUtils/first.h>

#include "RectangularMesh.h"

namespace finley {

Mesh* RectangularMesh_Hex20(const dim_t* numElements, const double* Length,
                            const bool* periodic, int order, int reduced_order,
                            bool useElementsOnFace, bool useFullElementOrder,
                            bool useMacroElements, bool optimize,
                            esysUtils::JMPI& mpiInfo)
{
    const int N_PER_E = 2;
    const int DIM = 3;
    dim_t Nstride0=0, Nstride1=0, Nstride2=0, local_NE0, local_NE1, local_NE2;
    index_t e_offset0, e_offset1, e_offset2;

    const Esys_MPI_rank myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    const dim_t NE0 = std::max(dim_t(1),numElements[0]);
    const dim_t NE1 = std::max(dim_t(1),numElements[1]);
    const dim_t NE2 = std::max(dim_t(1),numElements[2]);
    const dim_t N0 = N_PER_E*NE0+1;
    const dim_t N1 = N_PER_E*NE1+1;
    const dim_t N2 = N_PER_E*NE2+1;

    // allocate mesh
    std::stringstream name;
    name << "Brick " << N0 << " x " << N1 << " x " << N2;
    Mesh* out = new Mesh(name.str(), DIM, mpiInfo);

    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    bool generateAllNodes=(useFullElementOrder || useMacroElements);

    if (generateAllNodes) {
        if (useMacroElements) {
            refElements.reset(new ReferenceElementSet(Hex27Macro, order, reduced_order));
        } else {
            refElements.reset(new ReferenceElementSet(Hex27, order, reduced_order));
        }
        if (useElementsOnFace) {
            setError(SYSTEM_ERROR, "rich elements for Hex27 elements are not supported.");
            delete out;
            return NULL;
        } else {
            if (useMacroElements) {
                refFaceElements.reset(new ReferenceElementSet(Rec9Macro, order, reduced_order));
            } else {
                refFaceElements.reset(new ReferenceElementSet(Rec9, order, reduced_order));
            }
            refContactElements.reset(new ReferenceElementSet(Rec9_Contact, order, reduced_order));
        }
    } else { // !generateAllNodes
        refElements.reset(new ReferenceElementSet(Hex20, order, reduced_order));
        if (useElementsOnFace) {
            refFaceElements.reset(new ReferenceElementSet(Hex20Face, order, reduced_order));
            refContactElements.reset(new ReferenceElementSet(Hex20Face_Contact, order, reduced_order));
        } else {
            refFaceElements.reset(new ReferenceElementSet(Rec8, order, reduced_order));
            refContactElements.reset(new ReferenceElementSet(Rec8_Contact, order, reduced_order));
        }
    }
    refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));

    out->setPoints(new ElementFile(refPoints, mpiInfo));
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));
    out->setFaceElements(new ElementFile(refFaceElements, mpiInfo));
    out->setElements(new ElementFile(refElements, mpiInfo));

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
    if (!periodic[2] && local_NE2>0) {
        NDOF2=N2;
        if (offset2==0) NFaceElements+=local_NE1*local_NE0;
        if (local_NE2+e_offset2 == NE2) NFaceElements+=local_NE1*local_NE0;
    } else {
        NDOF2=N2-1;
    }

    if (!periodic[0] && local_NE0>0) {
        NDOF0=N0;
        if (e_offset0 == 0) NFaceElements+=local_NE1*local_NE2;
        if (local_NE0+e_offset0 == NE0) NFaceElements+=local_NE1*local_NE2;
    } else {
        NDOF0=N0-1;
    }
    if (!periodic[1] && local_NE1>0) {
        NDOF1=N1;
        if (e_offset1 == 0) NFaceElements+=local_NE0*local_NE2;
        if (local_NE1+e_offset1 == NE1) NFaceElements+=local_NE0*local_NE2;
    } else {
        NDOF1=N1-1;
    }

    // allocate tables
    out->Nodes->allocTable(local_N0*local_N1*local_N2);
    out->Elements->allocTable(local_NE0*local_NE1*local_NE2);
    out->FaceElements->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i2=0; i2<local_N2; i2++) {
        for (index_t i1=0; i1<local_N1; i1++) {
            for (index_t i0=0; i0<local_N0; i0++) {
                const dim_t k = i0+local_N0*i1+local_N0*local_N1*i2;
                const index_t global_i0 = i0+offset0;
                const index_t global_i1 = i1+offset1;
                const index_t global_i2 = i2+offset2;
                out->Nodes->Coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*Length[0];
                out->Nodes->Coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*Length[1];
                out->Nodes->Coordinates[INDEX2(2,k,DIM)]=DBLE(global_i2)/DBLE(N2-1)*Length[2];
                out->Nodes->Id[k]=Nstride0*global_i0+Nstride1*global_i1+Nstride2*global_i2;
                out->Nodes->Tag[k]=0;
                out->Nodes->globalDegreesOfFreedom[k]=Nstride0*(global_i0%NDOF0)
                                                +Nstride1*(global_i1%NDOF1)
                                                +Nstride2*(global_i2%NDOF2);
            }
        }
    }

    // set the elements
    dim_t NN = out->Elements->numNodes;
    index_t* eNodes = out->Elements->Nodes;
#pragma omp parallel for
    for (index_t i2=0; i2<local_NE2; i2++) {
        for (index_t i1=0; i1<local_NE1; i1++) {
            for (index_t i0=0; i0<local_NE0; i0++) {
                const dim_t k = i0+local_NE0*i1+local_NE0*local_NE1*i2;
                const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                    + Nstride1*N_PER_E*(i1+e_offset1)
                                    + Nstride2*N_PER_E*(i2+e_offset2);

                out->Elements->Id[k] = (i0+e_offset0)
                                     + NE0*(i1+e_offset1)
                                     + NE0*NE1*(i2+e_offset2);
                out->Elements->Tag[k]=0;
                out->Elements->Owner[k]=myRank;

                eNodes[INDEX2(0,k,NN)] =node0;
                eNodes[INDEX2(1,k,NN)] =node0+                      2*Nstride0;
                eNodes[INDEX2(2,k,NN)] =node0+           2*Nstride1+2*Nstride0;
                eNodes[INDEX2(3,k,NN)] =node0+           2*Nstride1;
                eNodes[INDEX2(4,k,NN)] =node0+2*Nstride2;
                eNodes[INDEX2(5,k,NN)] =node0+2*Nstride2           +2*Nstride0;
                eNodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                eNodes[INDEX2(7,k,NN)] =node0+2*Nstride2+2*Nstride1;
                eNodes[INDEX2(8,k,NN)] =node0+                      1*Nstride0;
                eNodes[INDEX2(9,k,NN)] =node0+           1*Nstride1+2*Nstride0;
                eNodes[INDEX2(10,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                eNodes[INDEX2(11,k,NN)]=node0+           1*Nstride1;
                eNodes[INDEX2(12,k,NN)]=node0+1*Nstride2;
                eNodes[INDEX2(13,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                eNodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                eNodes[INDEX2(15,k,NN)]=node0+1*Nstride2+2*Nstride1;
                eNodes[INDEX2(16,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                eNodes[INDEX2(17,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                eNodes[INDEX2(18,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                eNodes[INDEX2(19,k,NN)]=node0+2*Nstride2+1*Nstride1;
                if (generateAllNodes) {
                    eNodes[INDEX2(20,k,NN)]=node0+           1*Nstride1+1*Nstride0;
                    eNodes[INDEX2(21,k,NN)]=node0+1*Nstride2           +1*Nstride0;
                    eNodes[INDEX2(22,k,NN)]=node0+1*Nstride2+1*Nstride1+2*Nstride0;
                    eNodes[INDEX2(23,k,NN)]=node0+1*Nstride2+2*Nstride1+1*Nstride0;
                    eNodes[INDEX2(24,k,NN)]=node0+1*Nstride2+1*Nstride1;
                    eNodes[INDEX2(25,k,NN)]=node0+2*Nstride2+1*Nstride1+1*Nstride0;
                    eNodes[INDEX2(26,k,NN)]=node0+1*Nstride2+1*Nstride1+1*Nstride0;
                }
            }
        }
    }

    // face elements
    NN=out->FaceElements->numNodes;
    dim_t totalNECount=NE0*NE1*NE2;
    dim_t faceNECount = 0;
    eNodes = out->FaceElements->Nodes;

    // these are the quadrilateral elements on boundary 1 (x3=0):
    if (!periodic[2] && local_NE2>0) {
        // **  elements on boundary 100 (x3=0):
        if (offset2==0) {
#pragma omp parallel for
            for (index_t i1=0; i1<local_NE1; i1++) {
                for (index_t i0=0; i0<local_NE0; i0++) {
                    const dim_t k = i0+local_NE0*i1+faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(i1+e_offset1);

                    out->FaceElements->Id[k] = (i0+e_offset0)
                                             + NE0*(i1+e_offset1)+totalNECount;
                    out->FaceElements->Tag[k]=100;
                    out->FaceElements->Owner[k]=myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] =node0;
                        eNodes[INDEX2(1,k,NN)] =node0           +2*Nstride1;
                        eNodes[INDEX2(2,k,NN)] =node0           +2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)] =node0+           2*Nstride0;
                        eNodes[INDEX2(4,k,NN)] =node0+2*Nstride2;
                        eNodes[INDEX2(5,k,NN)] =node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)] =node0+2*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(8,k,NN)] =node0+           1*Nstride1;
                        eNodes[INDEX2(9,k,NN)] =node0+           2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(10,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(11,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(12,k,NN)]=node0+1*Nstride2;
                        eNodes[INDEX2(13,k,NN)]=node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(16,k,NN)]=node0+2*Nstride2+1*Nstride1;
                        eNodes[INDEX2(17,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(18,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(19,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)] =node0;
                        eNodes[INDEX2(1,k,NN)] =node0+           2*Nstride1;
                        eNodes[INDEX2(2,k,NN)] =node0+           2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)] =node0+                      2*Nstride0;
                        eNodes[INDEX2(4,k,NN)] =node0+           1*Nstride1;
                        eNodes[INDEX2(5,k,NN)] =node0+           2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(6,k,NN)] =node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)] =node0+                      1*Nstride0;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)]=node0+        1*Nstride1+1*Nstride0;
                        }
                    }
                }
            }
            faceNECount+=local_NE1*local_NE0;
        }
        totalNECount+=NE1*NE0;

        // **  elements on boundary 200 (x3=1):
        if (local_NE2+e_offset2 == NE2) {
#pragma omp parallel for
            for (index_t i1=0; i1<local_NE1; i1++) {
                for (index_t i0=0; i0<local_NE0; i0++) {
                    const dim_t k = i0+local_NE0*i1+faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(NE2-1);

                    out->FaceElements->Id[k] = (i0+e_offset0)
                                             + NE0*(i1+e_offset1)+totalNECount;
                    out->FaceElements->Tag[k]=200;
                    out->FaceElements->Owner[k]=myRank;
                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] =node0+2*Nstride2;
                        eNodes[INDEX2(1,k,NN)] =node0+2*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)] =node0+2*Nstride2+2*Nstride1;

                        eNodes[INDEX2(4,k,NN)] =node0;
                        eNodes[INDEX2(5,k,NN)] =node0+2*Nstride0;
                        eNodes[INDEX2(6,k,NN)] =node0+           2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)] =node0+           2*Nstride1;

                        eNodes[INDEX2(8,k,NN)] =node0+2*Nstride2+           1*Nstride0;
                        eNodes[INDEX2(9,k,NN)] =node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(10,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(11,k,NN)]=node0+2*Nstride2+1*Nstride1;

                        eNodes[INDEX2(12,k,NN)]=node0+1*Nstride2;
                        eNodes[INDEX2(13,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+1*Nstride2+2*Nstride1;

                        eNodes[INDEX2(16,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(17,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(18,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(19,k,NN)]=node0+           1*Nstride1;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)] =node0+2*Nstride2;
                        eNodes[INDEX2(1,k,NN)] =node0+2*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)] =node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(4,k,NN)] =node0+2*Nstride2           +1*Nstride0;
                        eNodes[INDEX2(5,k,NN)] =node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(7,k,NN)] =node0+2*Nstride2+1*Nstride1;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)]=node0+2*Nstride2+1*Nstride1+1*Nstride0;
                        }
                    }
                }
            }
            faceNECount+=local_NE1*local_NE0;
        }
        totalNECount+=NE1*NE0;
    } // !periodic[2] && local_NE2>0

    if (!periodic[0] && local_NE0>0) {
        // **  elements on boundary 001 (x1=0):
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i2=0; i2<local_NE2; i2++) {
                for (index_t i1=0; i1<local_NE1; i1++) {
                    const dim_t k = i1+local_NE1*i2+faceNECount;
                    const index_t node0 = Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);
                    out->FaceElements->Id[k] = (i1+e_offset1)
                                             + NE1*(i2+e_offset2)+totalNECount;
                    out->FaceElements->Tag[k]=1;
                    out->FaceElements->Owner[k]=myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)] =node0;
                        eNodes[INDEX2(1,k,NN)] =node0+2*Nstride2;
                        eNodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(3,k,NN)] =node0+2*Nstride1;

                        eNodes[INDEX2(4,k,NN)] =node0+2*Nstride0;
                        eNodes[INDEX2(5,k,NN)] =node0+2*Nstride2+2*Nstride0;
                        eNodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)] =node0+2*Nstride1+2*Nstride0;

                        eNodes[INDEX2(8,k,NN)] =node0+1*Nstride2;
                        eNodes[INDEX2(9,k,NN)] =node0+2*Nstride2+1*Nstride1;
                        eNodes[INDEX2(10,k,NN)]=node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(11,k,NN)]=node0+           1*Nstride1;

                        eNodes[INDEX2(12,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(13,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                        eNodes[INDEX2(14,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+2*Nstride1+           1*Nstride0;

                        eNodes[INDEX2(16,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(17,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(18,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(19,k,NN)]=node0+1*Nstride1+           2*Nstride0;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)] =node0;
                        eNodes[INDEX2(1,k,NN)] =node0+2*Nstride2;
                        eNodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(3,k,NN)] =node0+           2*Nstride1;
                        eNodes[INDEX2(4,k,NN)] =node0+1*Nstride2;
                        eNodes[INDEX2(5,k,NN)] =node0+2*Nstride2+1*Nstride1;
                        eNodes[INDEX2(6,k,NN)] =node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(7,k,NN)] =node0+           1*Nstride1;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)] =node0+1*Nstride2+1*Nstride1;
                        }
                    }
                }
            }
            faceNECount+=local_NE1*local_NE2;
        } // e_offset0 == 0
        totalNECount+=NE1*NE2;

        // **  elements on boundary 002 (x1=1):
        if (local_NE0+e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i2=0; i2<local_NE2; i2++) {
                for (index_t i1=0; i1<local_NE1; i1++) {
                    const dim_t k = i1+local_NE1*i2+faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(NE0-1)
                                        + Nstride1*N_PER_E*(i1+e_offset1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);
                    out->FaceElements->Id[k] = (i1+e_offset1)
                                             + NE1*(i2+e_offset2)+totalNECount;
                    out->FaceElements->Tag[k]=2;
                    out->FaceElements->Owner[k]=myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)]=node0+                      2*Nstride0;
                        eNodes[INDEX2(1,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+2*Nstride2+           2*Nstride0;

                        eNodes[INDEX2(4,k,NN)]=node0;
                        eNodes[INDEX2(5,k,NN)]=node0+           2*Nstride1;
                        eNodes[INDEX2(6,k,NN)]=node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(7,k,NN)]=node0+2*Nstride2;

                        eNodes[INDEX2(8,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(9,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(10,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(11,k,NN)]=node0+1*Nstride2+           2*Nstride0;

                        eNodes[INDEX2(12,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(13,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(14,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+2*Nstride2+           1*Nstride0;

                        eNodes[INDEX2(16,k,NN)]=node0+           1*Nstride1;
                        eNodes[INDEX2(17,k,NN)]=node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(18,k,NN)]=node0+2*Nstride2+1*Nstride1;
                        eNodes[INDEX2(19,k,NN)]=node0+1*Nstride2;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)]=node0                      +2*Nstride0;
                        eNodes[INDEX2(1,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(4,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(5,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(6,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)]=node0+1*Nstride2+1*Nstride1+2*Nstride0;
                        }
                    }
                }
            }
            faceNECount+=local_NE1*local_NE2;
        }
        totalNECount+=NE1*NE2;
    } // !periodic[0] && local_NE0>0

    if (!periodic[1] && local_NE1>0) {
        // **  elements on boundary 010 (x2=0):
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i2=0; i2<local_NE2; i2++) {
                for (index_t i0=0; i0<local_NE0; i0++) {
                    const dim_t k = i0+local_NE0*i2+faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride2*N_PER_E*(i2+e_offset2);

                    out->FaceElements->Id[k] = (i2+e_offset2)
                                             + NE2*(e_offset0+i0)+totalNECount;
                    out->FaceElements->Tag[k]=10;
                    out->FaceElements->Owner[k]=myRank;
                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)]=node0;
                        eNodes[INDEX2(1,k,NN)]=node0+                      2*Nstride0;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+2*Nstride2;

                        eNodes[INDEX2(4,k,NN)]=node0+           2*Nstride1;
                        eNodes[INDEX2(5,k,NN)]=node0+2*Nstride1+           2*Nstride0;
                        eNodes[INDEX2(6,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)]=node0+2*Nstride2+2*Nstride1;

                        eNodes[INDEX2(8,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(9,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(10,k,NN)]=node0+2*Nstride2+           1*Nstride0;
                        eNodes[INDEX2(11,k,NN)]=node0+1*Nstride2;

                        eNodes[INDEX2(12,k,NN)]=node0+           1*Nstride1;
                        eNodes[INDEX2(13,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(14,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+2*Nstride2+1*Nstride1;

                        eNodes[INDEX2(16,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(17,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(18,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(19,k,NN)]=node0+1*Nstride2+2*Nstride1;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)]=node0;
                        eNodes[INDEX2(1,k,NN)]=node0+                      2*Nstride0;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+2*Nstride2;
                        eNodes[INDEX2(4,k,NN)]=node0+                      1*Nstride0;
                        eNodes[INDEX2(5,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(6,k,NN)]=node0+2*Nstride2+           1*Nstride0;
                        eNodes[INDEX2(7,k,NN)]=node0+1*Nstride2;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)]=node0+1*Nstride2+         1*Nstride0;
                        }
                    }
                }
            }
            faceNECount+=local_NE0*local_NE2;
        } // e_offset1==0
        totalNECount+=NE0*NE2;

        // **  elements on boundary 020 (x2=1):
        if (local_NE1+e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i2=0; i2<local_NE2; i2++) {
                for (index_t i0=0; i0<local_NE0; i0++) {
                    const dim_t k = i0+local_NE0*i2+faceNECount;
                    const index_t node0 = Nstride0*N_PER_E*(i0+e_offset0)
                                        + Nstride1*N_PER_E*(NE1-1)
                                        + Nstride2*N_PER_E*(i2+e_offset2);

                    out->FaceElements->Id[k] = (i2+e_offset2)
                                             + NE2*(i0+e_offset0)+totalNECount;
                    out->FaceElements->Tag[k]=20;
                    out->FaceElements->Owner[k]=myRank;

                    if (useElementsOnFace) {
                        eNodes[INDEX2(0,k,NN)]=node0+           2*Nstride1;
                        eNodes[INDEX2(1,k,NN)]=node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+2*Nstride1+2*Nstride0;

                        eNodes[INDEX2(4,k,NN)]=node0;
                        eNodes[INDEX2(5,k,NN)]=node0+2*Nstride2;
                        eNodes[INDEX2(6,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                        eNodes[INDEX2(7,k,NN)]=node0+                      2*Nstride0;

                        eNodes[INDEX2(8,k,NN)]=node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(9,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(10,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(11,k,NN)]=node0+           2*Nstride1+1*Nstride0;

                        eNodes[INDEX2(12,k,NN)]=node0+           1*Nstride1;
                        eNodes[INDEX2(13,k,NN)]=node0+2*Nstride2+1*Nstride1;
                        eNodes[INDEX2(14,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                        eNodes[INDEX2(15,k,NN)]=node0+           1*Nstride1+2*Nstride0;

                        eNodes[INDEX2(16,k,NN)]=node0+1*Nstride2;
                        eNodes[INDEX2(17,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                        eNodes[INDEX2(18,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                        eNodes[INDEX2(19,k,NN)]=node0+                      1*Nstride0;
                    } else { // !useElementsOnFace
                        eNodes[INDEX2(0,k,NN)]=node0+           2*Nstride1;
                        eNodes[INDEX2(1,k,NN)]=node0+2*Nstride2+2*Nstride1;
                        eNodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(3,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(4,k,NN)]=node0+1*Nstride2+2*Nstride1;
                        eNodes[INDEX2(5,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                        eNodes[INDEX2(6,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                        eNodes[INDEX2(7,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                        if (generateAllNodes) {
                            eNodes[INDEX2(8,k,NN)]=node0+1*Nstride2+2*Nstride1+1*Nstride0;
                        }
                    }
                }
            }
            faceNECount+=local_NE0*local_NE2;
        }
        totalNECount+=NE0*NE2;
    }

    // add tag names
    out->addTagMap("top", 200);
    out->addTagMap("bottom", 100);
    out->addTagMap("left", 1);
    out->addTagMap("right", 2);
    out->addTagMap("front", 10);
    out->addTagMap("back", 20);

    // prepare mesh for further calculations
    out->resolveNodeIds();
    if (noError()) {
        out->prepare(optimize);
    }

    if (!noError()) {
        delete out;
        out=NULL;
    }
    return out;
}

} // namespace finley

