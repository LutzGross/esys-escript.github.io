
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
**
*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

#include <sstream>

/*
  input: 
    const double* coordinates[DIM*(*)]
    int numQuad
    const double* QuadWeights[numQuad]
    int numShape
    dim_t numElements
    int numNodes
    const index_t* nodes[numNodes*numElements]  where NUMSIDES*numShape<=numNodes
    const double* DSDv[numShape*DIM*numQuad]
    int numTest
    double* DTDv[LOCDIM*numTest*numQuad] 
    const index_t* elementId[numElements]

  output:
    double* dTdX[DIM*numTest*NUMSIDES*numQuad*numElements]
    double* volume[numQuad*numElements]
*/

namespace finley {

/****************************************************************************/
//
//  Jacobian 1D
//
void Assemble_jacobians_1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=1;
    const int LOCDIM=1;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double D=0.;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                D += X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
            }
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_1D: element " << e
                    << " (id " << elementId[e] << ") has length zero.";
                std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                for (int s=0; s<numTest; s++)
                     dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                         DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*invD;
            }
            volume[INDEX2(q,e,numQuad)]=std::abs(D)*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 2D with area element
//
void Assemble_jacobians_2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM = 2;
    const int LOCDIM = 2;
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        for (int q = 0; q < numQuad; q++) {
            double dXdv00 = 0.;
            double dXdv10 = 0.;
            double dXdv01 = 0.;
            double dXdv11 = 0.;
            for (int s = 0; s < numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
            }
            const double D = dXdv00*dXdv11 - dXdv01*dXdv10;
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_2D: element " << e
                    << " (id " << elementId[e] << ") has length zero.";
                throw FinleyException(ss.str());
            } else {
                const double invD = 1./D;
                const double dvdX00 = dXdv11*invD;
                const double dvdX10 =-dXdv10*invD;
                const double dvdX01 =-dXdv01*invD;
                const double dvdX11 = dXdv00*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                }
            }
            volume[INDEX2(q,e,numQuad)] = std::abs(D)*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 1D manifold in 2D and 1D elements
//
void Assemble_jacobians_2D_M1D_E1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=2;
    const int LOCDIM=1;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00=0.;
            double dXdv10=0.;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
            }
            const double D = dXdv00*dXdv00 + dXdv10*dXdv10;
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_2D_M1D_E1D: element " << e
                   << " (id " << elementId[e] << ") has length zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                const double dvdX00 = dXdv00*invD;
                const double dvdX01 = dXdv10*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01;
                }
                volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
            }
        }
    }
}

/****************************************************************************/
//
//  Jacobian 1D manifold in 2D and 1D elements with contact
//
void Assemble_jacobians_2D_M1D_E1D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=2;
    const int LOCDIM=1;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00_0=0;
            double dXdv10_0=0;
            double dXdv00_1=0;
            double dXdv10_1=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
            }
            const double D_0 = dXdv00_0*dXdv00_0 + dXdv10_0*dXdv10_0;
            const double D_1 = dXdv00_1*dXdv00_1 + dXdv10_1*dXdv10_1;
            if (D_0 == 0. || D_1 == 0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_2D_M1D_E1D_C: element " << e
                    << " (id " << elementId[e] << ") has length zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD_0 = 1./D_0;
                const double dvdX00_0=dXdv00_0*invD_0;
                const double dvdX01_0=dXdv10_0*invD_0;
                const double invD_1 = 1./D_1;
                const double dvdX00_1=dXdv00_1*invD_1;
                const double dvdX01_1=dXdv10_1*invD_1;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0;
                    dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0;
                    dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1;
                    dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1;
                }
                volume[INDEX2(q,e,numQuad)]=(sqrt(D_0)+sqrt(D_1))/2.*QuadWeights[q];
            }
        }
    }
}

/*****************************************************************************/
//
//  Jacobian 1D manifold in 2D and 2D elements
//
void Assemble_jacobians_2D_M1D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=2;
    const int LOCDIM=2;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00=0;
            double dXdv10=0;
            double dXdv01=0;
            double dXdv11=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
            }
            const double D = dXdv00*dXdv11 - dXdv01*dXdv10;
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_2D_M1D_E2D: element " << e
                    << " (id " << elementId[e] << ") has area zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                const double dvdX00 = dXdv11*invD;
                const double dvdX10 =-dXdv10*invD;
                const double dvdX01 =-dXdv01*invD;
                const double dvdX11 = dXdv00*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                }
            }
            volume[INDEX2(q,e,numQuad)]=sqrt(dXdv00*dXdv00+dXdv10*dXdv10)*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 1D manifold in 2D and 2D elements with contact
//
void Assemble_jacobians_2D_M1D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=2;
    const int LOCDIM=2;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00_0=0;
            double dXdv10_0=0;
            double dXdv01_0=0;
            double dXdv11_0=0;
            double dXdv00_1=0;
            double dXdv10_1=0;
            double dXdv01_1=0;
            double dXdv11_1=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
            }
            const double D_0 = dXdv00_0*dXdv11_0 - dXdv01_0*dXdv10_0;
            const double D_1 = dXdv00_1*dXdv11_1 - dXdv01_1*dXdv10_1;
            if (D_0 == 0. || D_1 == 0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_2D_M1D_E2D_C: element " << e
                    << " (id " << elementId[e] << ") has area zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD_0=1./D_0;
                const double dvdX00_0= dXdv11_0*invD_0;
                const double dvdX10_0=-dXdv10_0*invD_0;
                const double dvdX01_0=-dXdv01_0*invD_0;
                const double dvdX11_0= dXdv00_0*invD_0;
                const double invD_1=1./D_1;
                const double dvdX00_1= dXdv11_1*invD_1;
                const double dvdX10_1=-dXdv10_1*invD_1;
                const double dvdX01_1=-dXdv01_1*invD_1;
                const double dvdX11_1= dXdv00_1*invD_1;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0;
                    dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0;
                    dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1;
                    dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1;
                }
            }
            volume[INDEX2(q,e,numQuad)]=(sqrt(dXdv00_0*dXdv00_0+dXdv10_0*dXdv10_0)+sqrt(dXdv00_1*dXdv00_1+dXdv10_1*dXdv10_1))/2.*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 3D
//
void Assemble_jacobians_3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=3;
    const int LOCDIM=3;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00=0;
            double dXdv10=0;
            double dXdv20=0;
            double dXdv01=0;
            double dXdv11=0;
            double dXdv21=0;
            double dXdv02=0;
            double dXdv12=0;
            double dXdv22=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv02+=X0_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv12+=X1_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv22+=X2_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
            }
            const double D = dXdv00*(dXdv11*dXdv22-dXdv12*dXdv21)+ dXdv01*(dXdv20*dXdv12-dXdv10*dXdv22)+dXdv02*(dXdv10*dXdv21-dXdv20*dXdv11);
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_3D: element " << e
                    << " (id " << elementId[e] << ") has volume zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                const double dvdX00=(dXdv11*dXdv22-dXdv12*dXdv21)*invD;
                const double dvdX10=(dXdv20*dXdv12-dXdv10*dXdv22)*invD;
                const double dvdX20=(dXdv10*dXdv21-dXdv20*dXdv11)*invD;
                const double dvdX01=(dXdv02*dXdv21-dXdv01*dXdv22)*invD;
                const double dvdX11=(dXdv00*dXdv22-dXdv20*dXdv02)*invD;
                const double dvdX21=(dXdv01*dXdv20-dXdv00*dXdv21)*invD;
                const double dvdX02=(dXdv01*dXdv12-dXdv02*dXdv11)*invD;
                const double dvdX12=(dXdv02*dXdv10-dXdv00*dXdv12)*invD;
                const double dvdX22=(dXdv00*dXdv11-dXdv01*dXdv10)*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21;
                    dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22;
                }
                volume[INDEX2(q,e,numQuad)]=std::abs(D)*QuadWeights[q];
            }
        }
    }
}

/****************************************************************************/
//
//  Jacobian 2D manifold in 3D with 3D elements
//
void Assemble_jacobians_3D_M2D_E3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=3;
    const int LOCDIM=3;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00=0;
            double dXdv10=0;
            double dXdv20=0;
            double dXdv01=0;
            double dXdv11=0;
            double dXdv21=0;
            double dXdv02=0;
            double dXdv12=0;
            double dXdv22=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv02+=X0_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv12+=X1_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv22+=X2_loc*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
            }
            const double D = dXdv00*(dXdv11*dXdv22-dXdv12*dXdv21) +
                             dXdv01*(dXdv20*dXdv12-dXdv10*dXdv22) +
                             dXdv02*(dXdv10*dXdv21-dXdv20*dXdv11);
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_M2D_E3D: element " << e
                    << " (id " << elementId[e] << ") has volume zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                const double dvdX00=(dXdv11*dXdv22-dXdv12*dXdv21)*invD;
                const double dvdX10=(dXdv20*dXdv12-dXdv10*dXdv22)*invD;
                const double dvdX20=(dXdv10*dXdv21-dXdv20*dXdv11)*invD;
                const double dvdX01=(dXdv02*dXdv21-dXdv01*dXdv22)*invD;
                const double dvdX11=(dXdv00*dXdv22-dXdv20*dXdv02)*invD;
                const double dvdX21=(dXdv01*dXdv20-dXdv00*dXdv21)*invD;
                const double dvdX02=(dXdv01*dXdv12-dXdv02*dXdv11)*invD;
                const double dvdX12=(dXdv02*dXdv10-dXdv00*dXdv12)*invD;
                const double dvdX22=(dXdv00*dXdv11-dXdv01*dXdv10)*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21;
                    dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22;
                }
            }
            const double m0=dXdv10*dXdv21-dXdv20*dXdv11;
            const double m1=dXdv20*dXdv01-dXdv00*dXdv21;
            const double m2=dXdv00*dXdv11-dXdv10*dXdv01;
            volume[INDEX2(q,e,numQuad)]=sqrt(m0*m0+m1*m1+m2*m2)*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 2D manifold in 3D with 3D elements on contact
//
void Assemble_jacobians_3D_M2D_E3D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=3;
    const int LOCDIM=3;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00_0=0;
            double dXdv10_0=0;
            double dXdv20_0=0;
            double dXdv01_0=0;
            double dXdv11_0=0;
            double dXdv21_0=0;
            double dXdv02_0=0;
            double dXdv12_0=0;
            double dXdv22_0=0;
            double dXdv00_1=0;
            double dXdv10_1=0;
            double dXdv20_1=0;
            double dXdv01_1=0;
            double dXdv11_1=0;
            double dXdv21_1=0;
            double dXdv02_1=0;
            double dXdv12_1=0;
            double dXdv22_1=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X2_loc_0=coordinates[INDEX2(2,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X2_loc_1=coordinates[INDEX2(2,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20_0+=X2_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21_0+=X2_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv02_0+=X0_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv12_0+=X1_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv22_0+=X2_loc_0*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20_1+=X2_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21_1+=X2_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv02_1+=X0_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv12_1+=X1_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
                dXdv22_1+=X2_loc_1*DSDv[INDEX3(s,2,q,numShape,LOCDIM)];
            }

            const double D_0=dXdv00_0*(dXdv11_0*dXdv22_0-dXdv12_0*dXdv21_0) +
                             dXdv01_0*(dXdv20_0*dXdv12_0-dXdv10_0*dXdv22_0) +
                             dXdv02_0*(dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0);
            const double D_1=dXdv00_1*(dXdv11_1*dXdv22_1-dXdv12_1*dXdv21_1) +
                             dXdv01_1*(dXdv20_1*dXdv12_1-dXdv10_1*dXdv22_1) +
                             dXdv02_1*(dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1);
            if (D_0 == 0. || D_1 == 0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_M2D_E3D_C: element " << e
                    << " (id " << elementId[e] << ") has volume zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD_0=1./D_0;
                const double dvdX00_0=(dXdv11_0*dXdv22_0-dXdv12_0*dXdv21_0)*invD_0;
                const double dvdX10_0=(dXdv20_0*dXdv12_0-dXdv10_0*dXdv22_0)*invD_0;
                const double dvdX20_0=(dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0)*invD_0;
                const double dvdX01_0=(dXdv02_0*dXdv21_0-dXdv01_0*dXdv22_0)*invD_0;
                const double dvdX11_0=(dXdv00_0*dXdv22_0-dXdv20_0*dXdv02_0)*invD_0;
                const double dvdX21_0=(dXdv01_0*dXdv20_0-dXdv00_0*dXdv21_0)*invD_0;
                const double dvdX02_0=(dXdv01_0*dXdv12_0-dXdv02_0*dXdv11_0)*invD_0;
                const double dvdX12_0=(dXdv02_0*dXdv10_0-dXdv00_0*dXdv12_0)*invD_0;
                const double dvdX22_0=(dXdv00_0*dXdv11_0-dXdv01_0*dXdv10_0)*invD_0;
                const double invD_1=1./D_1;
                const double dvdX00_1=(dXdv11_1*dXdv22_1-dXdv12_1*dXdv21_1)*invD_1;
                const double dvdX10_1=(dXdv20_1*dXdv12_1-dXdv10_1*dXdv22_1)*invD_1;
                const double dvdX20_1=(dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1)*invD_1;
                const double dvdX01_1=(dXdv02_1*dXdv21_1-dXdv01_1*dXdv22_1)*invD_1;
                const double dvdX11_1=(dXdv00_1*dXdv22_1-dXdv20_1*dXdv02_1)*invD_1;
                const double dvdX21_1=(dXdv01_1*dXdv20_1-dXdv00_1*dXdv21_1)*invD_1;
                const double dvdX02_1=(dXdv01_1*dXdv12_1-dXdv02_1*dXdv11_1)*invD_1;
                const double dvdX12_1=(dXdv02_1*dXdv10_1-dXdv00_1*dXdv12_1)*invD_1;
                const double dvdX22_1=(dXdv00_1*dXdv11_1-dXdv01_1*dXdv10_1)*invD_1;

                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20_0;
                    dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21_0;
                    dTdX[INDEX4(        s,2,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_0 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22_0;
                    dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX20_1;
                    dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX21_1;
                    dTdX[INDEX4(numTest+s,2,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_1 +
                        DTDv[INDEX3(s,2,q,numTest,LOCDIM)]*dvdX22_1;
                }
            }
            const double m0_0=dXdv10_0*dXdv21_0-dXdv20_0*dXdv11_0;
            const double m1_0=dXdv20_0*dXdv01_0-dXdv00_0*dXdv21_0;
            const double m2_0=dXdv00_0*dXdv11_0-dXdv10_0*dXdv01_0;
            const double m0_1=dXdv10_1*dXdv21_1-dXdv20_1*dXdv11_1;
            const double m1_1=dXdv20_1*dXdv01_1-dXdv00_1*dXdv21_1;
            const double m2_1=dXdv00_1*dXdv11_1-dXdv10_1*dXdv01_1;
            volume[INDEX2(q,e,numQuad)]=(sqrt(m0_0*m0_0+m1_0*m1_0+m2_0*m2_0) +
                    sqrt(m0_1*m0_1+m1_1*m1_1+m2_1*m2_1))/2.*QuadWeights[q];
        }
    }
}

/****************************************************************************/
//
//  Jacobian 2D manifold in 3D with 2D elements
//
void Assemble_jacobians_3D_M2D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=3;
    const int LOCDIM=2;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00=0;
            double dXdv10=0;
            double dXdv20=0;
            double dXdv01=0;
            double dXdv11=0;
            double dXdv21=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc=coordinates[INDEX2(0,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X1_loc=coordinates[INDEX2(1,nodes[INDEX2(s,e,numNodes)],DIM)];
                const double X2_loc=coordinates[INDEX2(2,nodes[INDEX2(s,e,numNodes)],DIM)];
                dXdv00+=X0_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10+=X1_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20+=X2_loc*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01+=X0_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11+=X1_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21+=X2_loc*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
            }
            const double m00=dXdv00*dXdv00+dXdv10*dXdv10+dXdv20*dXdv20;
            const double m01=dXdv00*dXdv01+dXdv10*dXdv11+dXdv20*dXdv21;
            const double m11=dXdv01*dXdv01+dXdv11*dXdv11+dXdv21*dXdv21;
            const double D=m00*m11-m01*m01;
            if (D==0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_3D_M2D_E2D: element " << e
                    << " (id " << elementId[e] << ") has area zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD = 1./D;
                const double dvdX00=( m00*dXdv00-m01*dXdv01)*invD;
                const double dvdX01=( m00*dXdv10-m01*dXdv11)*invD;
                const double dvdX02=( m00*dXdv20-m01*dXdv21)*invD;
                const double dvdX10=(-m01*dXdv00+m11*dXdv01)*invD;
                const double dvdX11=(-m01*dXdv10+m11*dXdv11)*invD;
                const double dvdX12=(-m01*dXdv20+m11*dXdv21)*invD;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(s,0,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10;
                    dTdX[INDEX4(s,1,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11;
                    dTdX[INDEX4(s,2,q,e,numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12;
                }
                volume[INDEX2(q,e,numQuad)]=sqrt(D)*QuadWeights[q];
            }
        }
    }
}

/****************************************************************************/
//
//  Jacobian 2D manifold in 3D with 2D elements with contact
//
void Assemble_jacobians_3D_M2D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId)
{
    const int DIM=3;
    const int LOCDIM=2;
#pragma omp parallel for
    for (index_t e=0; e<numElements; e++) {
        for (int q=0; q<numQuad; q++) {
            double dXdv00_0=0;
            double dXdv10_0=0;
            double dXdv20_0=0;
            double dXdv01_0=0;
            double dXdv11_0=0;
            double dXdv21_0=0;
            double dXdv00_1=0;
            double dXdv10_1=0;
            double dXdv20_1=0;
            double dXdv01_1=0;
            double dXdv11_1=0;
            double dXdv21_1=0;
            for (int s=0; s<numShape; s++) {
                const double X0_loc_0=coordinates[INDEX2(0,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X1_loc_0=coordinates[INDEX2(1,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X2_loc_0=coordinates[INDEX2(2,nodes[INDEX2(s         ,e,numNodes)],DIM)];
                const double X0_loc_1=coordinates[INDEX2(0,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X1_loc_1=coordinates[INDEX2(1,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                const double X2_loc_1=coordinates[INDEX2(2,nodes[INDEX2(s+numShape,e,numNodes)],DIM)];
                dXdv00_0+=X0_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_0+=X1_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20_0+=X2_loc_0*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_0+=X0_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_0+=X1_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21_0+=X2_loc_0*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv00_1+=X0_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv10_1+=X1_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv20_1+=X2_loc_1*DSDv[INDEX3(s,0,q,numShape,LOCDIM)];
                dXdv01_1+=X0_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv11_1+=X1_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
                dXdv21_1+=X2_loc_1*DSDv[INDEX3(s,1,q,numShape,LOCDIM)];
            }
            const double m00_0=dXdv00_0*dXdv00_0+dXdv10_0*dXdv10_0+dXdv20_0*dXdv20_0;
            const double m01_0=dXdv00_0*dXdv01_0+dXdv10_0*dXdv11_0+dXdv20_0*dXdv21_0;
            const double m11_0=dXdv01_0*dXdv01_0+dXdv11_0*dXdv11_0+dXdv21_0*dXdv21_0;
            const double D_0=m00_0*m11_0-m01_0*m01_0;
            const double m00_1=dXdv00_1*dXdv00_1+dXdv10_1*dXdv10_1+dXdv20_1*dXdv20_1;
            const double m01_1=dXdv00_1*dXdv01_1+dXdv10_1*dXdv11_1+dXdv20_1*dXdv21_1;
            const double m11_1=dXdv01_1*dXdv01_1+dXdv11_1*dXdv11_1+dXdv21_1*dXdv21_1;
            const double D_1=m00_1*m11_1-m01_1*m01_1;
            if (D_0 == 0. || D_1 == 0.) {
                std::stringstream ss;
                ss << "Assemble_jacobians_3D_M2D_E2D_C: element " << e
                    << " (id " << elementId[e] << ") has area zero.";
                const std::string errorMsg = ss.str();
                throw FinleyException(errorMsg);
            } else {
                const double invD_0=1./D_0;
                const double dvdX00_0=( m00_0*dXdv00_0-m01_0*dXdv01_0)*invD_0;
                const double dvdX01_0=( m00_0*dXdv10_0-m01_0*dXdv11_0)*invD_0;
                const double dvdX02_0=( m00_0*dXdv20_0-m01_0*dXdv21_0)*invD_0;
                const double dvdX10_0=(-m01_0*dXdv00_0+m11_0*dXdv01_0)*invD_0;
                const double dvdX11_0=(-m01_0*dXdv10_0+m11_0*dXdv11_0)*invD_0;
                const double dvdX12_0=(-m01_0*dXdv20_0+m11_0*dXdv21_0)*invD_0;
                const double invD_1=1./D_1;
                const double dvdX00_1=( m00_1*dXdv00_1-m01_1*dXdv01_1)*invD_1;
                const double dvdX01_1=( m00_1*dXdv10_1-m01_1*dXdv11_1)*invD_1;
                const double dvdX02_1=( m00_1*dXdv20_1-m01_1*dXdv21_1)*invD_1;
                const double dvdX10_1=(-m01_1*dXdv00_1+m11_1*dXdv01_1)*invD_1;
                const double dvdX11_1=(-m01_1*dXdv10_1+m11_1*dXdv11_1)*invD_1;
                const double dvdX12_1=(-m01_1*dXdv20_1+m11_1*dXdv21_1)*invD_1;
                for (int s=0; s<numTest; s++) {
                    dTdX[INDEX4(        s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_0;
                    dTdX[INDEX4(        s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_0;
                    dTdX[INDEX4(        s,2,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_0 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_0;
                    dTdX[INDEX4(numTest+s,0,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX00_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX10_1;
                    dTdX[INDEX4(numTest+s,1,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX01_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX11_1;
                    dTdX[INDEX4(numTest+s,2,q,e,2*numTest,DIM,numQuad)] =
                        DTDv[INDEX3(s,0,q,numTest,LOCDIM)]*dvdX02_1 +
                        DTDv[INDEX3(s,1,q,numTest,LOCDIM)]*dvdX12_1;
                }
                volume[INDEX2(q,e,numQuad)]=(sqrt(D_0)+sqrt(D_1))/2.*QuadWeights[q];
            }
        }
    }
}

} // namespace finley

