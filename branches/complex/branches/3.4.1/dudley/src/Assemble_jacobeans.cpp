
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* Unless the loops in here get complicated again, this file should be compiled with loop unrolling */

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
index_t* element_id[numElements]

output:

double* dTdX[DIM*numTest*NUMSIDES*numQuad*numElements]
double* volume[numQuad*numElements]

*/

#include "ShapeTable.h"

#define SCALING(_nsub_,_dim_) pow(1./(double)(_nsub_),1./(double)(_dim_))

/************************************************************************************/
/*                                                            */
/*  Jacobean 2D with area element                             */
/*                                                            */
void Dudley_Assemble_jacobeans_2D(double *coordinates, dim_t numQuad, dim_t numElements, dim_t numNodes, index_t * nodes,
			   double *dTdX, double *absD, double *quadweight, index_t * element_id)
{
#define DIM 2
#define LOCDIM 2
    register int e, q;
    char error_msg[LenErrorMsg_MAX];
    const dim_t numTest = 3;	/* hoping this is used in constant folding */
    *quadweight = (numQuad == 1) ? 1. / 2 : 1. / 6;	/* numQuad is 1 or 3 */
#pragma omp parallel
    {
	register double dXdv00, dXdv10, dXdv01, dXdv11, dvdX00, dvdX10, dvdX01, dvdX11, D, invD;
#pragma omp for private(e,q, dXdv00,dXdv10,dXdv01,dXdv11,dvdX00,dvdX10,dvdX01,dvdX11, D,invD) schedule(static)
	for (e = 0; e < numElements; e++)
	{
#define COMPDXDV0(P)  coordinates[INDEX2(P,nodes[INDEX2(0,e,numNodes)],DIM)]*(-1)+\
coordinates[INDEX2(P,nodes[INDEX2(1,e,numNodes)],DIM)]*1+\
coordinates[INDEX2(P,nodes[INDEX2(2,e,numNodes)],DIM)]*(0)

#define COMPDXDV1(P)  coordinates[INDEX2(P,nodes[INDEX2(0,e,numNodes)],DIM)]*(-1)+\
coordinates[INDEX2(P,nodes[INDEX2(1,e,numNodes)],DIM)]*(0)+\
coordinates[INDEX2(P,nodes[INDEX2(2,e,numNodes)],DIM)]*(1)

	    dXdv00 = 0;
	    dXdv10 = 0;
	    dXdv01 = 0;
	    dXdv11 = 0;
	    dXdv00 = COMPDXDV0(0);
	    dXdv10 = COMPDXDV0(1);
	    dXdv01 = COMPDXDV1(0);
	    dXdv11 = COMPDXDV1(1);
	    D = dXdv00 * dXdv11 - dXdv01 * dXdv10;
	    absD[e] = ABS(D);
	    if (D == 0.)
	    {
		sprintf(error_msg, "Dudley_Assemble_jacobeans_2D: element %d (id %d) has area zero.", e, element_id[e]);
		Dudley_setError(ZERO_DIVISION_ERROR, error_msg);
	    }
	    else
	    {
		invD = 1. / D;
		dvdX00 = dXdv11 * invD;
		dvdX10 = -dXdv10 * invD;
		dvdX01 = -dXdv01 * invD;
		dvdX11 = dXdv00 * invD;
		if (numQuad == 1)
		{
		    dTdX[INDEX4(0, 0, 0, e, numTest, DIM, numQuad)] = DTDV_2D[0][0] * dvdX00 + DTDV_2D[1][1] * dvdX10;
		    dTdX[INDEX4(1, 0, 0, e, numTest, DIM, numQuad)] = DTDV_2D[0][1] * dvdX00 + DTDV_2D[1][0] * dvdX10;
		    dTdX[INDEX4(2, 0, 0, e, numTest, DIM, numQuad)] = DTDV_2D[2][0] * dvdX00 + DTDV_2D[2][1] * dvdX10;

		    dTdX[INDEX4(0, 1, 0, e, numTest, DIM, numQuad)] = DTDV_2D[0][0] * dvdX01 + DTDV_2D[1][1] * dvdX11;
		    dTdX[INDEX4(1, 1, 0, e, numTest, DIM, numQuad)] = DTDV_2D[0][1] * dvdX01 + DTDV_2D[1][0] * dvdX11;
		    dTdX[INDEX4(2, 1, 0, e, numTest, DIM, numQuad)] = DTDV_2D[2][0] * dvdX01 + DTDV_2D[2][1] * dvdX11;

		}
		else		/* numQuad==3 */
		{
		    for (q = 0; q < numTest; ++q)	/* relying on unroll loops to optimise this */
		    {
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
	}
    }				/* end parallel */
#undef DIM
#undef LOCDIM
#undef DTDXSET
#undef COMPDXDV0
#undef COMPDXDV1
}

/************************************************************************************/
/*                                                            */
/*  Jacobean 1D manifold in 2D and 1D elements                */
/*                                                            */
void Dudley_Assemble_jacobeans_2D_M1D_E1D(double *coordinates, dim_t numQuad,
				   dim_t numElements, dim_t numNodes, index_t * nodes,
				   double *dTdX, double *absD, double *quadweight, index_t * element_id)
{
#define DIM 2
#define LOCDIM 1
    register int e;
    char error_msg[LenErrorMsg_MAX];
    const dim_t numTest = 2;
    *quadweight = (numQuad == 1) ? 1.0 : 0.5;
    /* numQuad is 1 or 2 */
#pragma omp parallel
    {
	register double dXdv00, dXdv10, dvdX00, dvdX01, D, invD;
#pragma omp for private(e,dXdv00,dXdv10,dvdX00,dvdX01,D,invD) schedule(static)
	for (e = 0; e < numElements; e++)
	{
	    dXdv00 = 0;
	    dXdv10 = 0;
	    dXdv00 +=
		coordinates[INDEX2(0, nodes[INDEX2(0, e, numNodes)], DIM)] * (-1.) +
		coordinates[INDEX2(0, nodes[INDEX2(1, e, numNodes)], DIM)];
	    dXdv00 +=
		coordinates[INDEX2(1, nodes[INDEX2(0, e, numNodes)], DIM)] * (-1.) +
		coordinates[INDEX2(1, nodes[INDEX2(1, e, numNodes)], DIM)];
	    D = dXdv00 * dXdv00 + dXdv10 * dXdv10;
	    if (D == 0.)
	    {
		sprintf(error_msg, "Dudley_Assemble_jacobeans_2D_M1D_E1D: element %d (id %d) has length zero.", e,
			element_id[e]);
		Dudley_setError(ZERO_DIVISION_ERROR, error_msg);
	    }
	    else
	    {
		invD = 1. / D;
		dvdX00 = dXdv00 * invD;
		dvdX01 = dXdv10 * invD;
		/* The number of quad points is 1 or 2 */
		dTdX[INDEX4(0, 0, 0, e, numTest, DIM, numQuad)] = -1 * dvdX00;
		dTdX[INDEX4(0, 1, 0, e, numTest, DIM, numQuad)] = -1 * dvdX01;
		dTdX[INDEX4(1, 0, 0, e, numTest, DIM, numQuad)] = -1 * dvdX00;
		dTdX[INDEX4(1, 1, 0, e, numTest, DIM, numQuad)] = -1 * dvdX01;
		absD[e] = sqrt(D);
		if (numQuad == 2)
		{
		    dTdX[INDEX4(0, 0, 1, e, numTest, DIM, numQuad)] = dvdX00;
		    dTdX[INDEX4(0, 1, 1, e, numTest, DIM, numQuad)] = dvdX01;
		    dTdX[INDEX4(1, 0, 1, e, numTest, DIM, numQuad)] = dvdX00;
		    dTdX[INDEX4(1, 1, 1, e, numTest, DIM, numQuad)] = dvdX01;
		}
	    }
	}
    }				/* end parallel */
#undef DIM
#undef LOCDIM
}

/************************************************************************************/
/*                                                            */
/*  Jacobean 3D                                               */
/*                                                            */
void Dudley_Assemble_jacobeans_3D(double *coordinates, dim_t numQuad, dim_t numElements, dim_t numNodes, index_t * nodes,
			   double *dTdX, double *absD, double *quadweight, index_t * element_id)
{
#define DIM 3
#define LOCDIM 3
    int e, q, s;
    char error_msg[LenErrorMsg_MAX];
    /* numQuad is 1 or 4 */
    const dim_t numShape = 4, numTest = 4;
    *quadweight = (numQuad == 1) ? 1. / 6 : 1. / 24;

#pragma omp parallel
    {
	register double dXdv00, dXdv10, dXdv20, dXdv01, dXdv11, dXdv21, dXdv02, dXdv12, dXdv22,
	    dvdX00, dvdX10, dvdX20, dvdX01, dvdX11, dvdX21, dvdX02, dvdX12, dvdX22, D, invD, X0_loc, X1_loc, X2_loc;
#pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,dXdv02,dXdv12,dXdv22,dvdX00,dvdX10,dvdX20,dvdX01,dvdX11,dvdX21,dvdX02,dvdX12,dvdX22,D,invD,X0_loc,X1_loc,X2_loc) schedule(static)
	for (e = 0; e < numElements; e++)
	{
	    dXdv00 = 0;
	    dXdv10 = 0;
	    dXdv20 = 0;
	    dXdv01 = 0;
	    dXdv11 = 0;
	    dXdv21 = 0;
	    dXdv02 = 0;
	    dXdv12 = 0;
	    dXdv22 = 0;
	    for (s = 0; s < numShape; s++)
	    {
		X0_loc = coordinates[INDEX2(0, nodes[INDEX2(s, e, numNodes)], DIM)];
		X1_loc = coordinates[INDEX2(1, nodes[INDEX2(s, e, numNodes)], DIM)];
		X2_loc = coordinates[INDEX2(2, nodes[INDEX2(s, e, numNodes)], DIM)];
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
	    D = dXdv00 * (dXdv11 * dXdv22 - dXdv12 * dXdv21) + dXdv01 * (dXdv20 * dXdv12 - dXdv10 * dXdv22) +
		dXdv02 * (dXdv10 * dXdv21 - dXdv20 * dXdv11);
	    absD[e] = ABS(D);
	    if (D == 0.)
	    {
		sprintf(error_msg, "Dudley_Assemble_jacobeans_3D: element %d (id %d) has volume zero.", e, element_id[e]);
		Dudley_setError(ZERO_DIVISION_ERROR, error_msg);
	    }
	    else
	    {
		invD = 1. / D;
		dvdX00 = (dXdv11 * dXdv22 - dXdv12 * dXdv21) * invD;
		dvdX10 = (dXdv20 * dXdv12 - dXdv10 * dXdv22) * invD;
		dvdX20 = (dXdv10 * dXdv21 - dXdv20 * dXdv11) * invD;
		dvdX01 = (dXdv02 * dXdv21 - dXdv01 * dXdv22) * invD;
		dvdX11 = (dXdv00 * dXdv22 - dXdv20 * dXdv02) * invD;
		dvdX21 = (dXdv01 * dXdv20 - dXdv00 * dXdv21) * invD;
		dvdX02 = (dXdv01 * dXdv12 - dXdv02 * dXdv11) * invD;
		dvdX12 = (dXdv02 * dXdv10 - dXdv00 * dXdv12) * invD;
		dvdX22 = (dXdv00 * dXdv11 - dXdv01 * dXdv10) * invD;
		for (q = 0; q < numQuad; q++)
		{
		    for (s = 0; s < numTest; s++)
		    {
			dTdX[INDEX4(s, 0, q, e, numTest, DIM, numQuad)] =
			    DTDV_3D[s][0] * dvdX00 + DTDV_3D[s][1] * dvdX10 + DTDV_3D[s][2] * dvdX20;
			dTdX[INDEX4(s, 1, q, e, numTest, DIM, numQuad)] =
			    DTDV_3D[s][0] * dvdX01 + DTDV_3D[s][1] * dvdX11 + DTDV_3D[s][2] * dvdX21;
			dTdX[INDEX4(s, 2, q, e, numTest, DIM, numQuad)] =
			    DTDV_3D[s][0] * dvdX02 + DTDV_3D[s][1] * dvdX12 + DTDV_3D[s][2] * dvdX22;
		    }
		}
	    }
	}
    }				/* end parallel */
#undef DIM
#undef LOCDIM
}

/************************************************************************************/
/*                                                            */
/*  Jacobean 2D manifold in 3D with 2D elements               */
/*                                                            */
void Dudley_Assemble_jacobeans_3D_M2D_E2D(double *coordinates, dim_t numQuad, dim_t numElements, dim_t numNodes,
				   index_t * nodes, double *dTdX, double *absD, double *quadweight,
				   index_t * element_id)
{
#define DIM 3
#define LOCDIM 2
    register int e, q, s;
    char error_msg[LenErrorMsg_MAX];
    const double DTDV[3][2] = { {-1., -1.}, {1., 0.}, {0., 1.} };
    const dim_t numShape = 3, numTest = 3;
    /* numQuad is 1 or 3 */
    *quadweight = (numQuad == 1) ? 1. / 2 : 1. / 6;
#pragma omp parallel
    {
	register double dXdv00, dXdv10, dXdv20, dXdv01, dXdv11, dXdv21, m00, m01, m11,
	    dvdX00, dvdX01, dvdX02, dvdX10, dvdX11, dvdX12, D, invD, X0_loc, X1_loc, X2_loc;
#pragma omp for private(e,q,s,dXdv00,dXdv10,dXdv20,dXdv01,dXdv11,dXdv21,m00,m01,m11,dvdX00,dvdX01,dvdX02,dvdX10,dvdX11,dvdX12,D,invD, X0_loc, X1_loc, X2_loc) schedule(static)
	for (e = 0; e < numElements; e++)
	{
	    dXdv00 = 0;
	    dXdv10 = 0;
	    dXdv20 = 0;
	    dXdv01 = 0;
	    dXdv11 = 0;
	    dXdv21 = 0;
	    for (s = 0; s < numShape; s++)
	    {
		X0_loc = coordinates[INDEX2(0, nodes[INDEX2(s, e, numNodes)], DIM)];
		X1_loc = coordinates[INDEX2(1, nodes[INDEX2(s, e, numNodes)], DIM)];
		X2_loc = coordinates[INDEX2(2, nodes[INDEX2(s, e, numNodes)], DIM)];
		dXdv00 += X0_loc * DTDV[s][0];
		dXdv10 += X1_loc * DTDV[s][0];
		dXdv20 += X2_loc * DTDV[s][0];
		dXdv01 += X0_loc * DTDV[s][1];
		dXdv11 += X1_loc * DTDV[s][1];
		dXdv21 += X2_loc * DTDV[s][1];
	    }
	    m00 = dXdv00 * dXdv00 + dXdv10 * dXdv10 + dXdv20 * dXdv20;
	    m01 = dXdv00 * dXdv01 + dXdv10 * dXdv11 + dXdv20 * dXdv21;
	    m11 = dXdv01 * dXdv01 + dXdv11 * dXdv11 + dXdv21 * dXdv21;
	    D = m00 * m11 - m01 * m01;
	    absD[e] = sqrt(D);
	    if (D == 0.)
	    {
		sprintf(error_msg, "Dudley_Assemble_jacobeans_3D_M2D: element %d (id %d) has area zero.", e, element_id[e]);
		Dudley_setError(ZERO_DIVISION_ERROR, error_msg);
	    }
	    else
	    {
		invD = 1. / D;
		dvdX00 = (m00 * dXdv00 - m01 * dXdv01) * invD;
		dvdX01 = (m00 * dXdv10 - m01 * dXdv11) * invD;
		dvdX02 = (m00 * dXdv20 - m01 * dXdv21) * invD;
		dvdX10 = (-m01 * dXdv00 + m11 * dXdv01) * invD;
		dvdX11 = (-m01 * dXdv10 + m11 * dXdv11) * invD;
		dvdX12 = (-m01 * dXdv20 + m11 * dXdv21) * invD;
		for (q = 0; q < numQuad; q++)
		{
		    for (s = 0; s < numTest; s++)
		    {
			dTdX[INDEX4(s, 0, q, e, numTest, DIM, numQuad)] = DTDV[s][0] * dvdX00 + DTDV[s][1] * dvdX10;
			dTdX[INDEX4(s, 1, q, e, numTest, DIM, numQuad)] = DTDV[s][0] * dvdX01 + DTDV[s][1] * dvdX11;
			dTdX[INDEX4(s, 2, q, e, numTest, DIM, numQuad)] = DTDV[s][0] * dvdX02 + DTDV[s][1] * dvdX12;
		    }
		}
	    }
	}
    }				/* end parallel section */
#undef DIM
#undef LOCDIM
}
