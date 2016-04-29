
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "TriangularMesh.h"

#include <escript/index.h>

using escript::DataTypes::real_t;

#define MAX3(_arg1_,_arg2_,_arg3_) std::max(_arg1_,std::max(_arg2_,_arg3_))

namespace dudley {

// Be careful reading this function. The X? and NStride? are 1,2,3
// but the loop vars are 0,1,2
Mesh* TriangularMesh_Tet4(const dim_t* numElements, const double* length,
                          bool optimize, escript::JMPI mpiInfo)
{
    const int DIM = 3;
#ifdef Dudley_TRACE
    double time0 = Dudley_timer();
#endif

    const int LEFTTAG = 1;      /* boundary x1=0 */
    const int RIGHTTAG = 2;     /* boundary x1=1 */
    const int BOTTOMTAG = 100;  /* boundary x3=1 */
    const int TOPTAG = 200;     /* boundary x3=0 */
    const int FRONTTAG = 10;    /* boundary x2=0 */
    const int BACKTAG = 20;     /* boundary x2=1 */

    const int myRank = mpiInfo->rank;

    // set up the global dimensions of the mesh
    const dim_t NE0 = std::max(dim_t(1), numElements[0]);
    const dim_t NE1 = std::max(dim_t(1), numElements[1]);
    const dim_t NE2 = std::max(dim_t(1), numElements[2]);
    const dim_t N0 = NE0 + 1;
    const dim_t N1 = NE1 + 1;
    const dim_t N2 = NE2 + 1;

    // allocate mesh
    std::stringstream name;
    name << "Triangular " << N0 << " x " << N1 << " x " << N2 << " (x 5) mesh";
    Mesh* out = new Mesh(name.str(), DIM, mpiInfo);

    out->setPoints(new ElementFile(Dudley_Point1, mpiInfo));
    out->setFaceElements(new ElementFile(Dudley_Tri3, mpiInfo));
    out->setElements(new ElementFile(Dudley_Tet4, mpiInfo));

    dim_t Nstride0, Nstride1, Nstride2;
    dim_t local_NE0, local_NE1, local_NE2;
    index_t e_offset0, e_offset1, e_offset2;
    // work out the largest dimension
    if (N2 == MAX3(N0, N1, N2)) {
        Nstride0 = 1;
        Nstride1 = N0;
        Nstride2 = N0 * N1;
        local_NE0 = NE0;
        e_offset0 = 0;
        local_NE1 = NE1;
        e_offset1 = 0;
        mpiInfo->split(NE2, &local_NE2, &e_offset2);
    } else if (N1 == MAX3(N0, N1, N2)) {
        Nstride0 = N2;
        Nstride1 = N0 * N2;
        Nstride2 = 1;
        local_NE0 = NE0;
        e_offset0 = 0;
        mpiInfo->split(NE1, &local_NE1, &e_offset1);
        local_NE2 = NE2;
        e_offset2 = 0;
    } else {
        Nstride0 = N1 * N2;
        Nstride1 = 1;
        Nstride2 = N1;
        mpiInfo->split(NE0, &local_NE0, &e_offset0);
        local_NE1 = NE1;
        e_offset1 = 0;
        local_NE2 = NE2;
        e_offset2 = 0;
    }
    const index_t offset0 = e_offset0;
    const index_t offset1 = e_offset1;
    const index_t offset2 = e_offset2;
    const dim_t local_N0 = local_NE0 > 0 ? local_NE0 + 1 : 0;
    const dim_t local_N1 = local_NE1 > 0 ? local_NE1 + 1 : 0;
    const dim_t local_N2 = local_NE2 > 0 ? local_NE2 + 1 : 0;

    // get the number of surface elements
    dim_t NFaceElements = 0;
    dim_t NDOF0, NDOF1, NDOF2;
    if (local_NE2 > 0) {
        NDOF2 = N2;
        if (offset2 == 0)
            NFaceElements += 2 * local_NE1 * local_NE0; // each face is split
        if (local_NE2 + e_offset2 == NE2)
            NFaceElements += 2 * local_NE1 * local_NE0;
    } else {
        NDOF2 = N2 - 1;
    }

    if (local_NE0 > 0) {
        NDOF0 = N0;
        if (e_offset0 == 0)
            NFaceElements += 2 * local_NE1 * local_NE2;
        if (local_NE0 + e_offset0 == NE0)
            NFaceElements += 2 * local_NE1 * local_NE2;
    } else {
        NDOF0 = N0 - 1;
    }

    if (local_NE1 > 0) {
        NDOF1 = N1;
        if (e_offset1 == 0)
            NFaceElements += 2 * local_NE0 * local_NE2;
        if (local_NE1 + e_offset1 == NE1)
            NFaceElements += 2 * local_NE0 * local_NE2;
    } else {
        NDOF1 = N1 - 1;
    }

    // allocate tables
    out->Nodes->allocTable(local_N0 * local_N1 * local_N2);
    // we split the rectangular prism this code used to produce into 5
    // tetrahedra
    out->Elements->allocTable(local_NE0 * local_NE1 * local_NE2 * 5);
    // each border face will be split in half
    out->FaceElements->allocTable(NFaceElements);

    // create nodes
#pragma omp parallel for
    for (index_t i2 = 0; i2 < local_N2; i2++) {
        for (index_t i1 = 0; i1 < local_N1; i1++) {
            for (index_t i0 = 0; i0 < local_N0; i0++) {
                const index_t k = i0 + local_N0 * i1 + local_N0 * local_N1 * i2;
                const index_t global_i0 = i0 + offset0;
                const index_t global_i1 = i1 + offset1;
                const index_t global_i2 = i2 + offset2;
                out->Nodes->Coordinates[INDEX2(0, k, DIM)] = (real_t)global_i0 / (real_t)(N0 - 1) * length[0];
                out->Nodes->Coordinates[INDEX2(1, k, DIM)] = (real_t)global_i1 / (real_t)(N1 - 1) * length[1];
                out->Nodes->Coordinates[INDEX2(2, k, DIM)] = (real_t)global_i2 / (real_t)(N2 - 1) * length[2];
                out->Nodes->Id[k] = Nstride0 * global_i0 + Nstride1 * global_i1 + Nstride2 * global_i2;
                out->Nodes->Tag[k] = 0;
                out->Nodes->globalDegreesOfFreedom[k] =
                                    Nstride0 * (global_i0 % NDOF0)
                                    + Nstride1 * (global_i1 % NDOF1)
                                    + Nstride2 * (global_i2 % NDOF2);
            }
        }
    }

    // set the elements
    // If we are not the only rank we may need to shift our pattern to match
    // neighbours
    int global_adjustment = (offset0 + offset1 + offset2) % 2;

    int NN = out->Elements->numNodes;
#pragma omp parallel for
    for (index_t i2 = 0; i2 < local_NE2; i2++) {
        for (index_t i1 = 0; i1 < local_NE1; i1++) {
            for (index_t i0 = 0; i0 < local_NE0; i0++) {
                const index_t k = 5 * (i0 + local_NE0 * i1 + local_NE0 * local_NE1 * i2);
                const index_t node0 = Nstride0 * (i0 + e_offset0)
                                    + Nstride1 * (i1 + e_offset1)
                                    + Nstride2 * (i2 + e_offset2);

                const index_t res = 5 * ((i0 + e_offset0)
                                  + NE0 * (i1 + e_offset1)
                                  + NE0 * NE1 * (i2 + e_offset2));
                for (int j = 0; j < 5; ++j) {
                    out->Elements->Id[k + j] = res + j;
                    out->Elements->Tag[k + j] = 0;
                    out->Elements->Owner[k + j] = myRank;
                }

                // in non-rotated orientation the points are numbered as
                // follows:
                // The bottom face (anticlockwise = 0,1,3,2),
                // top face (anticlockwise 4,5,7,6)
                index_t v[8];
                if ((global_adjustment + i0 + i1 + i2) % 2 == 0) {
                    v[0] = node0;
                    v[1] = node0 + Nstride0;
                    v[2] = node0 + Nstride1;
                    v[3] = node0 + Nstride1 + Nstride0;
                    v[4] = node0 + Nstride2;
                    v[5] = node0 + Nstride0 + Nstride2;
                    v[6] = node0 + Nstride1 + Nstride2;
                    v[7] = node0 + Nstride2 + Nstride1 + Nstride0;
                } else {
                    // this form is rotated around the 0,2,4,6 face clockwise
                    // 90 degrees
                    v[0] = node0 + Nstride1; // node 0 ends up in position 2
                    v[2] = node0 + Nstride1 + Nstride2; // node 2 ends up in position 6
                    v[6] = node0 + Nstride2; // node 6 ends up in position 4
                    v[4] = node0; // node 4 ends up in position 0
                    v[1] = node0 + Nstride0 + Nstride1; // node 1 -> pos 3
                    v[3] = node0 + Nstride2 + Nstride1 + Nstride0; // node 3 -> pos 7
                    v[7] = node0 + Nstride0 + Nstride2; // node 7 -> pos 5
                    v[5] = node0 + Nstride0; // node 5 -> pos 1
                }

                // elements nodes are numbered: centre, x, y, z
                out->Elements->Nodes[INDEX2(0, k, NN)] = v[4];
                out->Elements->Nodes[INDEX2(1, k, NN)] = v[5];
                out->Elements->Nodes[INDEX2(2, k, NN)] = v[6];
                out->Elements->Nodes[INDEX2(3, k, NN)] = v[0];

                out->Elements->Nodes[INDEX2(0, k + 1, NN)] = v[7];
                out->Elements->Nodes[INDEX2(1, k + 1, NN)] = v[6];
                out->Elements->Nodes[INDEX2(2, k + 1, NN)] = v[5];
                out->Elements->Nodes[INDEX2(3, k + 1, NN)] = v[3];

                out->Elements->Nodes[INDEX2(0, k + 2, NN)] = v[2];
                out->Elements->Nodes[INDEX2(1, k + 2, NN)] = v[3];
                out->Elements->Nodes[INDEX2(2, k + 2, NN)] = v[0];
                out->Elements->Nodes[INDEX2(3, k + 2, NN)] = v[6];

                out->Elements->Nodes[INDEX2(0, k + 3, NN)] = v[1];
                out->Elements->Nodes[INDEX2(1, k + 3, NN)] = v[0];
                out->Elements->Nodes[INDEX2(2, k + 3, NN)] = v[3];
                out->Elements->Nodes[INDEX2(3, k + 3, NN)] = v[5];

                // I can't work out where the center is for this one
                out->Elements->Nodes[INDEX2(0, k + 4, NN)] = v[5];
                out->Elements->Nodes[INDEX2(1, k + 4, NN)] = v[0];
                out->Elements->Nodes[INDEX2(2, k + 4, NN)] = v[6];
                out->Elements->Nodes[INDEX2(3, k + 4, NN)] = v[3];
            }
        }
    } // for all elements

    // face elements
    NN = out->FaceElements->numNodes;
    dim_t totalNECount = 5 * NE0 * NE1 * NE2;
    dim_t faceNECount = 0;

    // these are the quadrilateral elements on boundary 1 (x3=0)
    if (local_NE2 > 0) {
        // ** elements on boundary 100 (x3=0)
        if (e_offset2 == 0) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const index_t k = 2 * (i0 + local_NE0 * i1) + faceNECount;
                    const index_t node0 = Nstride0 * (i0 + e_offset0)
                                        + Nstride1 * (i1 + e_offset1);
                    const index_t res = 2 * (i0 + e_offset0)
                                     + NE0 * (i1 + e_offset1) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = BOTTOMTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = BOTTOMTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n0 = node0;
                    const index_t n1 = node0 + Nstride0;
                    const index_t n2 = node0 + Nstride1;
                    const index_t n3 = node0 + Nstride0 + Nstride1;

                    if ((global_adjustment + i0 + i1) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n3;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n1;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n3;

                    } else {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n1;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n3;

                    }
                }
            }
            faceNECount += 2 * local_NE1 * local_NE0;
        }
        totalNECount += 2 * NE1 * NE0;
        // ** elements on boundary 200 (x3=1) - Top
        if (local_NE2 + e_offset2 == NE2) {
#pragma omp parallel for
            for (index_t i1 = 0; i1 < local_NE1; i1++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const index_t k = 2 * (i0 + local_NE0 * i1) + faceNECount;
                    const index_t node0 = Nstride0 * (i0 + e_offset0)
                                        + Nstride1 * (i1 + e_offset1)
                                        + Nstride2 * (NE2 - 1);

                    const index_t res = 2 * (i0 + e_offset0)
                                      + NE0 * (i1 + e_offset1) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = TOPTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = TOPTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n4 = node0 + Nstride2;
                    const index_t n5 = node0 + Nstride0 + Nstride2;
                    const index_t n6 = node0 + Nstride1 + Nstride2;
                    const index_t n7 = node0 + Nstride1 + Nstride0 + Nstride2;

                    if ((global_adjustment + i0 + i1 + local_NE2 - 1) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n5;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n6;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n5;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n6;
                    } else {

                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n5;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n7;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n6;
                    }
                }
            }
            faceNECount += 2 * local_NE1 * local_NE0;
        }
        totalNECount += 2 * NE1 * NE0;
    }

    if (local_NE0 > 0) {
        // ** elements on boundary 001 (x1=0) - Left
        if (e_offset0 == 0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i1 = 0; i1 < local_NE1; i1++) {
                    const index_t k = 2 * (i1 + local_NE1 * i2) + faceNECount;
                    const index_t node0 = Nstride1 * (i1 + e_offset1)
                                        + Nstride2 * (i2 + e_offset2);
                    const index_t res = 2 * (i1 + e_offset1)
                                      + NE1 * (i2 + e_offset2) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = LEFTTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = LEFTTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n0 = node0;
                    const index_t n2 = node0 + Nstride1;
                    const index_t n4 = node0 + Nstride2;
                    const index_t n6 = node0 + Nstride1 + Nstride2;

                    if ((global_adjustment + 0 + i1 + i2) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n6;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n6;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n2;
                    } else {
                        // this form is rotated around the 0,2,4,6 face
                        // clockwise 90 degrees
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n2;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n4;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n6;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n2;
                    }
                }
            }
            faceNECount += 2 * local_NE1 * local_NE2;
        }
        totalNECount += 2 * NE1 * NE2;
        // ** elements on boundary 002 (x1=1) - Right
        if (local_NE0 + e_offset0 == NE0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i1 = 0; i1 < local_NE1; i1++) {
                    const index_t k = 2 * (i1 + local_NE1 * i2) + faceNECount;
                    const index_t node0 = Nstride0 * (NE0 - 1)
                                        + Nstride1 * (i1 + e_offset1)
                                        + Nstride2 * (i2 + e_offset2);
                    const index_t res = 2 * (i1 + e_offset1)
                                      + NE1 * (i2 + e_offset2) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = RIGHTTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = RIGHTTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n1 = node0 + Nstride0;
                    const index_t n3 = node0 + Nstride0 + Nstride1;
                    const index_t n5 = node0 + Nstride0 + Nstride2;
                    const index_t n7 = node0 + Nstride0 + Nstride1 + Nstride2;

                    if ((global_adjustment + local_NE0 - 1 + i1 + i2) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n3;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n5;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n3;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n5;
                    } else {
                        // this form is rotated around the 0,2,4,6 face
                        // clockwise 90 degrees
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n5;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n3;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n7;
                    }
                }
            }
            faceNECount += 2 * local_NE1 * local_NE2;
        }
        totalNECount += 2 * NE1 * NE2;
    }
    if (local_NE1 > 0) {
        // ** elements on boundary 010 (x2=0) - Front
        if (e_offset1 == 0) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const index_t k = 2 * (i0 + local_NE0 * i2) + faceNECount;
                    const index_t node0 = Nstride0 * (i0 + e_offset0)
                                        + Nstride2 * (i2 + e_offset2);
                    const index_t res = 2 * (i2 + e_offset2)
                                      + NE2 * (e_offset0 + i0) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = FRONTTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = FRONTTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n0 = node0;
                    const index_t n1 = node0 + Nstride0;
                    const index_t n4 = node0 + Nstride2;
                    const index_t n5 = node0 + Nstride0 + Nstride2;

                    if ((global_adjustment + i0 + 0 + i2) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n5;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n5;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n4;

                    } else {
                        // this form is rotated around the 0,2,4,6 face
                        // clockwise 90 degrees
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n0;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n4;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n1;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n5;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n4;

                    }
                }
            }
            faceNECount += 2 * local_NE0 * local_NE2;
        }
        totalNECount += 2 * NE0 * NE2;
        // ** elements on boundary 020 (x2=1) - Back
        if (local_NE1 + e_offset1 == NE1) {
#pragma omp parallel for
            for (index_t i2 = 0; i2 < local_NE2; i2++) {
                for (index_t i0 = 0; i0 < local_NE0; i0++) {
                    const index_t k = 2 * (i0 + local_NE0 * i2) + faceNECount;
                    const index_t node0 = Nstride0 * (i0 + e_offset0)
                                        + Nstride1 * (NE1 - 1)
                                        + Nstride2 * (i2 + e_offset2);
                    const index_t res = 2 * (i2 + e_offset2)
                                      + NE2 * (i0 + e_offset0) + totalNECount;
                    out->FaceElements->Id[k] = res;
                    out->FaceElements->Tag[k] = BACKTAG;
                    out->FaceElements->Owner[k] = myRank;
                    out->FaceElements->Id[k + 1] = res + 1;
                    out->FaceElements->Tag[k + 1] = BACKTAG;
                    out->FaceElements->Owner[k + 1] = myRank;

                    const index_t n2 = node0 + Nstride1;
                    const index_t n6 = node0 + Nstride1 + Nstride2;
                    const index_t n7 = node0 + Nstride0 + Nstride1 + Nstride2;
                    const index_t n3 = node0 + Nstride0 + Nstride1;

                    if ((global_adjustment + i0 + local_NE1 - 1 + i2) % 2 == 0) {
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n6;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n3;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n6;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n3;

                    } else {
                        // this form is rotated around the 0,2,4,6 face
                        // clockwise 90 degrees
                        out->FaceElements->Nodes[INDEX2(0, k, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(1, k, NN)] = n6;
                        out->FaceElements->Nodes[INDEX2(2, k, NN)] = n7;

                        out->FaceElements->Nodes[INDEX2(0, k + 1, NN)] = n2;
                        out->FaceElements->Nodes[INDEX2(1, k + 1, NN)] = n7;
                        out->FaceElements->Nodes[INDEX2(2, k + 1, NN)] = n3;
                    }
                }
            }
            faceNECount += 2 * local_NE0 * local_NE2;
        }
        totalNECount += 2 * NE0 * NE2;
    }

    // add tag names
    out->addTagMap("top", TOPTAG);
    out->addTagMap("bottom", BOTTOMTAG);
    out->addTagMap("left", LEFTTAG);
    out->addTagMap("right", RIGHTTAG);
    out->addTagMap("front", FRONTTAG);
    out->addTagMap("back", BACKTAG);

    // prepare mesh for further calculations
    out->resolveNodeIds();
    out->prepare(optimize);
    return out;
}

} // namespace dudley

