
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

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

#include <escript/index.h>

// Unless the loops in here get complicated again, this file should be
// compiled with loop unrolling

/* input: 

double* coordinates[DIM*(*)]
dim_t numQuad
double* QuadWeights[numQuad]
dim_t numShape
dim_t numElements
dim_t numNodes
index_t* nodes[numNodes*numElements]  where NUMSIDES*numShape<=numNodes
double* DSDv[numShape*DIM*numQuad]
dim_t numTest
double* DTDv[LOCDIM*numTest*numQuad] 
index_t* elementId[numElements]

output:

double* dTdX[DIM*numTest*NUMSIDES*numQuad*numElements]
double* volume[numQuad*numElements]

*/

#define SCALING(_nsub_,_dim_) pow(1./(double)(_nsub_),1./(double)(_dim_))

namespace dudley {

/****************************************************************************/
//
//  Jacobian 2D with area element
//
void Assemble_jacobians_2D(const double* coordinates, int numQuad,
                       dim_t numElements, int numNodes, const index_t* nodes,
                       double* dTdX, double* absD, double* quadweight,
                       const index_t* elementId)
{
    const int DIM = 2;
    const int numTest = 3; // hoping this is used in constant folding
    *quadweight = (numQuad == 1) ? 1. / 2 : 1. / 6; // numQuad is 1 or 3
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
#define COMPDXDV0(P) coordinates[INDEX2(P,nodes[INDEX2(0,e,numNodes)],DIM)]*(-1)+\
coordinates[INDEX2(P,nodes[INDEX2(1,e,numNodes)],DIM)]*1+\
coordinates[INDEX2(P,nodes[INDEX2(2,e,numNodes)],DIM)]*(0)

#define COMPDXDV1(P)  coordinates[INDEX2(P,nodes[INDEX2(0,e,numNodes)],DIM)]*(-1)+\
coordinates[INDEX2(P,nodes[INDEX2(1,e,numNodes)],DIM)]*(0)+\
coordinates[INDEX2(P,nodes[INDEX2(2,e,numNodes)],DIM)]*(1)

        double dXdv00 = COMPDXDV0(0);
        double dXdv10 = COMPDXDV0(1);
        double dXdv01 = COMPDXDV1(0);
        double dXdv11 = COMPDXDV1(1);
        const double D = dXdv00 * dXdv11 - dXdv01 * dXdv10;
        absD[e] = std::abs(D);
        if (D == 0.) {
            std::stringstream ss;
            ss << "Assemble_jacobians_2D: element " << e
                << " (id " << elementId[e] << ") has area zero.";
            throw DudleyException(ss.str());
        } else {
            const double invD = 1. / D;
            const double dvdX00 = dXdv11 * invD;
            const double dvdX10 = -dXdv10 * invD;
            const double dvdX01 = -dXdv01 * invD;
            const double dvdX11 = dXdv00 * invD;
            if (numQuad == 1) {
                dTdX[INDEX4(0, 0, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[0][0] * dvdX00 + DTDV_2D[1][1] * dvdX10;
                dTdX[INDEX4(1, 0, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[0][1] * dvdX00 + DTDV_2D[1][0] * dvdX10;
                dTdX[INDEX4(2, 0, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[2][0] * dvdX00 + DTDV_2D[2][1] * dvdX10;

                dTdX[INDEX4(0, 1, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[0][0] * dvdX01 + DTDV_2D[1][1] * dvdX11;
                dTdX[INDEX4(1, 1, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[0][1] * dvdX01 + DTDV_2D[1][0] * dvdX11;
                dTdX[INDEX4(2, 1, 0, e, numTest, DIM, numQuad)] =
                    DTDV_2D[2][0] * dvdX01 + DTDV_2D[2][1] * dvdX11;

            } else { // numQuad == 3
                // relying on unroll loops to optimise this
                for (int q = 0; q < numTest; ++q) {
                    dTdX[INDEX4(0, 0, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[0][0] * dvdX00 + DTDV_2D[1][1] * dvdX10;
                    dTdX[INDEX4(1, 0, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[0][1] * dvdX00 + DTDV_2D[1][0] * dvdX10;
                    dTdX[INDEX4(2, 0, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[2][0] * dvdX00 + DTDV_2D[2][1] * dvdX10;

                    dTdX[INDEX4(0, 1, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[0][0] * dvdX01 + DTDV_2D[1][1] * dvdX11;
                    dTdX[INDEX4(1, 1, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[0][1] * dvdX01 + DTDV_2D[1][0] * dvdX11;
                    dTdX[INDEX4(2, 1, q, e, numTest, DIM, numQuad)] =
                        DTDV_2D[2][0] * dvdX01 + DTDV_2D[2][1] * dvdX11;

                }
            }
        }
    } // end parallel for
#undef COMPDXDV0
#undef COMPDXDV1
}

//
// Jacobian 1D manifold in 2D and 1D elements
//
void Assemble_jacobians_2D_M1D_E1D(const double* coordinates, int numQuad,
                                   dim_t numElements, int numNodes,
                                   const index_t* nodes, double* dTdX,
                                   double* absD, double* quadweight,
                                   const index_t* elementId)
{
    const int DIM = 2;
    const int numTest = 2;
    *quadweight = (numQuad == 1) ? 1.0 : 0.5; // numQuad is 1 or 2
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        double dXdv00 =
            coordinates[INDEX2(0, nodes[INDEX2(0, e, numNodes)], DIM)] * (-1.) +
            coordinates[INDEX2(0, nodes[INDEX2(1, e, numNodes)], DIM)];
        double dXdv10 =
            coordinates[INDEX2(1, nodes[INDEX2(0, e, numNodes)], DIM)] * (-1.) +
            coordinates[INDEX2(1, nodes[INDEX2(1, e, numNodes)], DIM)];
        const double D = dXdv00 * dXdv00 + dXdv10 * dXdv10;
        if (D == 0.) {
            std::stringstream ss;
            ss << "Assemble_jacobians_2D_M1D_E1D: element " << e
                << " (id " << elementId[e] << ") has length zero.";
            throw DudleyException(ss.str());
        } else {
            const double invD = 1. / D;
            const double dvdX00 = dXdv00 * invD;
            const double dvdX01 = dXdv10 * invD;
            // The number of quad points is 1 or 2
            dTdX[INDEX4(0, 0, 0, e, numTest, DIM, numQuad)] = -1 * dvdX00;
            dTdX[INDEX4(0, 1, 0, e, numTest, DIM, numQuad)] = -1 * dvdX01;
            dTdX[INDEX4(1, 0, 0, e, numTest, DIM, numQuad)] = -1 * dvdX00;
            dTdX[INDEX4(1, 1, 0, e, numTest, DIM, numQuad)] = -1 * dvdX01;
            absD[e] = sqrt(D);
            if (numQuad == 2) {
                dTdX[INDEX4(0, 0, 1, e, numTest, DIM, numQuad)] = dvdX00;
                dTdX[INDEX4(0, 1, 1, e, numTest, DIM, numQuad)] = dvdX01;
                dTdX[INDEX4(1, 0, 1, e, numTest, DIM, numQuad)] = dvdX00;
                dTdX[INDEX4(1, 1, 1, e, numTest, DIM, numQuad)] = dvdX01;
            }
        }
    } // end parallel for
}

//
// Jacobian 3D
//
void Assemble_jacobians_3D(const double* coordinates, int numQuad,
                           dim_t numElements, int numNodes,
                           const index_t* nodes, double* dTdX, double* absD,
                           double* quadweight, const index_t* elementId)
{
    const int DIM = 3;
    const int numShape = 4;
    const int numTest = 4;
    *quadweight = (numQuad == 1) ? 1. / 6 : 1. / 24; // numQuad is 1 or 4

#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        double dXdv00 = 0;
        double dXdv10 = 0;
        double dXdv20 = 0;
        double dXdv01 = 0;
        double dXdv11 = 0;
        double dXdv21 = 0;
        double dXdv02 = 0;
        double dXdv12 = 0;
        double dXdv22 = 0;
        for (int s = 0; s < numShape; s++) {
            const double X0_loc = coordinates[INDEX2(0, nodes[INDEX2(s, e, numNodes)], DIM)];
            const double X1_loc = coordinates[INDEX2(1, nodes[INDEX2(s, e, numNodes)], DIM)];
            const double X2_loc = coordinates[INDEX2(2, nodes[INDEX2(s, e, numNodes)], DIM)];
            dXdv00 += X0_loc * DTDV_3D[s][0];
            dXdv10 += X1_loc * DTDV_3D[s][0];
            dXdv20 += X2_loc * DTDV_3D[s][0];
            dXdv01 += X0_loc * DTDV_3D[s][1];
            dXdv11 += X1_loc * DTDV_3D[s][1];
            dXdv21 += X2_loc * DTDV_3D[s][1];
            dXdv02 += X0_loc * DTDV_3D[s][2];
            dXdv12 += X1_loc * DTDV_3D[s][2];
            dXdv22 += X2_loc * DTDV_3D[s][2];
        }
        const double D = dXdv00 * (dXdv11 * dXdv22 - dXdv12 * dXdv21)
                       + dXdv01 * (dXdv20 * dXdv12 - dXdv10 * dXdv22)
                       + dXdv02 * (dXdv10 * dXdv21 - dXdv20 * dXdv11);
        absD[e] = std::abs(D);
        if (D == 0.) {
            std::stringstream ss;
            ss << "Assemble_jacobians_3D: element " << e
                << " (id " << elementId[e] << ") has volume zero.";
            throw DudleyException(ss.str());
        } else {
            const double invD = 1. / D;
            const double dvdX00 = (dXdv11 * dXdv22 - dXdv12 * dXdv21) * invD;
            const double dvdX10 = (dXdv20 * dXdv12 - dXdv10 * dXdv22) * invD;
            const double dvdX20 = (dXdv10 * dXdv21 - dXdv20 * dXdv11) * invD;
            const double dvdX01 = (dXdv02 * dXdv21 - dXdv01 * dXdv22) * invD;
            const double dvdX11 = (dXdv00 * dXdv22 - dXdv20 * dXdv02) * invD;
            const double dvdX21 = (dXdv01 * dXdv20 - dXdv00 * dXdv21) * invD;
            const double dvdX02 = (dXdv01 * dXdv12 - dXdv02 * dXdv11) * invD;
            const double dvdX12 = (dXdv02 * dXdv10 - dXdv00 * dXdv12) * invD;
            const double dvdX22 = (dXdv00 * dXdv11 - dXdv01 * dXdv10) * invD;
            for (int q = 0; q < numQuad; q++) {
                for (int s = 0; s < numTest; s++) {
                    dTdX[INDEX4(s, 0, q, e, numTest, DIM, numQuad)] =
                        DTDV_3D[s][0] * dvdX00 + DTDV_3D[s][1] * dvdX10
                        + DTDV_3D[s][2] * dvdX20;
                    dTdX[INDEX4(s, 1, q, e, numTest, DIM, numQuad)] =
                        DTDV_3D[s][0] * dvdX01 + DTDV_3D[s][1] * dvdX11
                        + DTDV_3D[s][2] * dvdX21;
                    dTdX[INDEX4(s, 2, q, e, numTest, DIM, numQuad)] =
                        DTDV_3D[s][0] * dvdX02 + DTDV_3D[s][1] * dvdX12
                        + DTDV_3D[s][2] * dvdX22;
                }
            }
        }
    } // end parallel for
}

//
// Jacobian 2D manifold in 3D with 2D elements
//
void Assemble_jacobians_3D_M2D_E2D(const double* coordinates, int numQuad,
                                   dim_t numElements, int numNodes,
                                   const index_t* nodes, double* dTdX,
                                   double* absD, double* quadweight,
                                   const index_t* elementId)
{
    const int DIM = 3;
    const double DTDV[3][2] = { {-1., -1.}, {1., 0.}, {0., 1.} };
    const int numShape = 3;
    const int numTest = 3;
    *quadweight = (numQuad == 1) ? 1. / 2 : 1. / 6; // numQuad is 1 or 3
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        double dXdv00 = 0;
        double dXdv10 = 0;
        double dXdv20 = 0;
        double dXdv01 = 0;
        double dXdv11 = 0;
        double dXdv21 = 0;
        for (int s = 0; s < numShape; s++) {
            const double X0_loc = coordinates[INDEX2(0, nodes[INDEX2(s, e, numNodes)], DIM)];
            const double X1_loc = coordinates[INDEX2(1, nodes[INDEX2(s, e, numNodes)], DIM)];
            const double X2_loc = coordinates[INDEX2(2, nodes[INDEX2(s, e, numNodes)], DIM)];
            dXdv00 += X0_loc * DTDV[s][0];
            dXdv10 += X1_loc * DTDV[s][0];
            dXdv20 += X2_loc * DTDV[s][0];
            dXdv01 += X0_loc * DTDV[s][1];
            dXdv11 += X1_loc * DTDV[s][1];
            dXdv21 += X2_loc * DTDV[s][1];
        }
        const double m00 = dXdv00 * dXdv00 + dXdv10 * dXdv10 + dXdv20 * dXdv20;
        const double m01 = dXdv00 * dXdv01 + dXdv10 * dXdv11 + dXdv20 * dXdv21;
        const double m11 = dXdv01 * dXdv01 + dXdv11 * dXdv11 + dXdv21 * dXdv21;
        const double D = m00 * m11 - m01 * m01;
        absD[e] = sqrt(D);
        if (D == 0.) {
            std::stringstream ss;
            ss << "Assemble_jacobians_3D_M2D: element " << e
                << " (id " << elementId[e] << ") has area zero.";
            throw DudleyException(ss.str());
        } else {
            const double invD = 1. / D;
            const double dvdX00 = (m00 * dXdv00 - m01 * dXdv01) * invD;
            const double dvdX01 = (m00 * dXdv10 - m01 * dXdv11) * invD;
            const double dvdX02 = (m00 * dXdv20 - m01 * dXdv21) * invD;
            const double dvdX10 = (-m01 * dXdv00 + m11 * dXdv01) * invD;
            const double dvdX11 = (-m01 * dXdv10 + m11 * dXdv11) * invD;
            const double dvdX12 = (-m01 * dXdv20 + m11 * dXdv21) * invD;
            for (int q = 0; q < numQuad; q++) {
                for (int s = 0; s < numTest; s++) {
                    dTdX[INDEX4(s, 0, q, e, numTest, DIM, numQuad)] =
                        DTDV[s][0] * dvdX00 + DTDV[s][1] * dvdX10;
                    dTdX[INDEX4(s, 1, q, e, numTest, DIM, numQuad)] =
                        DTDV[s][0] * dvdX01 + DTDV[s][1] * dvdX11;
                    dTdX[INDEX4(s, 2, q, e, numTest, DIM, numQuad)] =
                        DTDV[s][0] * dvdX02 + DTDV[s][1] * dvdX12;
                }
            }
        }
    } // end parallel for
}

} // namespace dudley

