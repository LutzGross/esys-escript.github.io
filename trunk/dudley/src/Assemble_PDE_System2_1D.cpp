
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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

/************************************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F  */
/*    the shape functions for test and solution must be identical */

/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  and -(X_{k,i})_i + Y_k */

/*    u has p.numComp components in a 1D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                                                  */

/*    Shape of the coefficients: */

/*      A = p.numEqu x 1 x p.numComp x 1 */
/*      B = 1 x numEqu x p.numComp  */
/*      C = p.numEqu x 1 x p.numComp  */
/*      D = p.numEqu x p.numComp  */
/*      X = p.numEqu x 1  */
/*      Y = p.numEqu   */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Dudley_Assemble_PDE_System2_1D(Dudley_Assemble_Parameters p, Dudley_ElementFile * elements,
				    paso::SystemMatrix_ptr Mat, escriptDataC * F,
				    escriptDataC * A, escriptDataC * B, escriptDataC * C, escriptDataC * D,
				    escriptDataC * X, escriptDataC * Y)
{

#define DIM 1
    index_t color;
    dim_t e;
    __const double *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
    double *EM_S, *EM_F, *DSDX;
    index_t *row_index;
    register dim_t q, s, r, k, m;
    register double rtmp;
    bool add_EM_F, add_EM_S;

    bool extendedA = isExpanded(A);
    bool extendedB = isExpanded(B);
    bool extendedC = isExpanded(C);
    bool extendedD = isExpanded(D);
    bool extendedX = isExpanded(X);
    bool extendedY = isExpanded(Y);
    double *F_p = (requireWrite(F), getSampleDataRW(F, 0));	/* use comma, to get around the mixed code and declarations thing */
    double *S = p.row_jac->BasisFunctions->S;
    dim_t len_EM_S = p.row_numShapes * p.row_numShapes * p.numEqu * p.numComp;
    dim_t len_EM_F = p.row_numShapes * p.numEqu;

#pragma omp parallel private(color, EM_S, EM_F, Vol, DSDX, A_p, B_p, C_p, D_p, X_p, Y_p, A_q, B_q, C_q, D_q, X_q, Y_q,row_index, q, s,r,k,m,rtmp,add_EM_F, add_EM_S)
    {
	EM_S = new  double[len_EM_S];
	EM_F = new  double[len_EM_F];
	row_index = new  index_t[p.row_numShapes];

	if (!Dudley_checkPtr(EM_S) && !Dudley_checkPtr(EM_F) && !Dudley_checkPtr(row_index))
	{

	    for (color = elements->minColor; color <= elements->maxColor; color++)
	    {
		/*  open loop over all elements: */
#pragma omp for private(e) schedule(static)
		for (e = 0; e < elements->numElements; e++)
		{
		    if (elements->Color[e] == color)
		    {

			A_p = getSampleDataRO(A, e);
			B_p = getSampleDataRO(B, e);
			C_p = getSampleDataRO(C, e);
			D_p = getSampleDataRO(D, e);
			X_p = getSampleDataRO(X, e);
			Y_p = getSampleDataRO(Y, e);
			double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
			Vol = &(p.row_jac->volume[INDEX3(0, 0, e, p.numQuadTotal, 1)]);
			DSDX = &(p.row_jac->DSDX[INDEX5(0, 0, 0, 0, e, p.row_numShapes, DIM, p.numQuadTotal, 1)]);
			for (q = 0; q < len_EM_S; ++q)
			    EM_S[q] = 0;
			for (q = 0; q < len_EM_F; ++q)
			    EM_F[q] = 0;
			add_EM_F = FALSE;
			add_EM_S = FALSE;

		      /************************************************************************************/
			/*   process A: */
		      /************************************************************************************/

			if (NULL != A_p)
			{
			    add_EM_S = TRUE;
			    if (extendedA)
			    {
				A_q = &(A_p[INDEX6(0, 0, 0, 0, 0, 0, p.numEqu, DIM, p.numComp, DIM, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0.;
						for (q = 0; q < p.numQuadTotal; q++)
						{
						    rtmp += vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)] *
							A_q[INDEX5(k, 0, m, 0, q, p.numEqu, DIM, p.numComp, DIM)] *
							DSDX[INDEX3(r, 0, q, p.row_numShapes, DIM)];
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)] *
						DSDX[INDEX3(r, 0, q, p.row_numShapes, DIM)];
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] +=
						    rtmp * A_p[INDEX4(k, 0, m, 0, p.numEqu, DIM, p.numComp)];
					    }
					}
				    }
				}
			    }
			}
		      /************************************************************************************/
			/*   process B: */
		      /************************************************************************************/

			if (NULL != B_p)
			{
			    add_EM_S = TRUE;
			    if (extendedB)
			    {
				B_q = &(B_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, DIM, p.numComp, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0.;
						for (q = 0; q < p.numQuadTotal; q++)
						{
						    rtmp += vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)] *
							B_q[INDEX4(k, 0, m, q, p.numEqu, DIM, p.numComp)] *
							S[INDEX2(r, q, p.row_numShapes)];
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)] *
						S[INDEX2(r, q, p.row_numShapes)];
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] +=
						    rtmp * B_p[INDEX3(k, 0, m, p.numEqu, DIM)];
					    }
					}
				    }
				}
			    }
			}
		      /************************************************************************************/
			/*   process C: */
		      /************************************************************************************/

			if (NULL != C_p)
			{
			    add_EM_S = TRUE;
			    if (extendedC)
			    {
				C_q = &(C_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, p.numComp, DIM, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuadTotal; q++)
						{
						    rtmp +=
							vol * S[INDEX2(s, q, p.row_numShapes)] *
							C_q[INDEX4(k, m, 0, q, p.numEqu, p.numComp, DIM)] *
							DSDX[INDEX3(r, 0, q, p.row_numShapes, DIM)];
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * S[INDEX2(s, q, p.row_numShapes)] *
						DSDX[INDEX3(r, 0, q, p.row_numShapes, DIM)];
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] +=
						    rtmp * C_p[INDEX3(k, m, 0, p.numEqu, p.numComp)];
					    }
					}
				    }
				}
			    }
			}
		      /*********************************************************************************** */
			/* process D */
		      /************************************************************************************/

			if (NULL != D_p)
			{
			    add_EM_S = TRUE;
			    if (extendedD)
			    {
				D_q = &(D_p[INDEX4(0, 0, 0, 0, p.numEqu, p.numComp, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuadTotal; q++)
						{
						    rtmp +=
							vol * S[INDEX2(s, q, p.row_numShapes)] *
							D_q[INDEX3(k, m, q, p.numEqu, p.numComp)] *
							S[INDEX2(r, q, p.row_numShapes)];
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] += rtmp;

					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (r = 0; r < p.row_numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * S[INDEX2(s, q, p.row_numShapes)] *
						S[INDEX2(r, q, p.row_numShapes)];
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.row_numShapes)] +=
						    rtmp * D_p[INDEX2(k, m, p.numEqu)];
					    }
					}
				    }
				}
			    }
			}
		      /************************************************************************************/
			/*   process X: */
		      /************************************************************************************/

			if (NULL != X_p)
			{
			    add_EM_F = TRUE;
			    if (extendedX)
			    {
				X_q = &(X_p[INDEX4(0, 0, 0, 0, p.numEqu, DIM, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (k = 0; k < p.numEqu; k++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)] *
						X_q[INDEX3(k, 0, q, p.numEqu, DIM)];
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp;
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    rtmp = 0;
				    for (q = 0; q < p.numQuadTotal; q++)
					rtmp += vol * DSDX[INDEX3(s, 0, q, p.row_numShapes, DIM)];
				    for (k = 0; k < p.numEqu; k++)
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp * X_p[INDEX2(k, 0, p.numEqu)];
				}
			    }
			}
		     /************************************************************************************/
			/*   process Y: */
		     /************************************************************************************/

			if (NULL != Y_p)
			{
			    add_EM_F = TRUE;
			    if (extendedY)
			    {
				Y_q = &(Y_p[INDEX3(0, 0, 0, p.numEqu, p.numQuadTotal)]);
				for (s = 0; s < p.row_numShapes; s++)
				{
				    for (k = 0; k < p.numEqu; k++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuadTotal; q++)
					    rtmp +=
						vol * S[INDEX2(s, q, p.row_numShapes)] * Y_q[INDEX2(k, q, p.numEqu)];
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp;
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.row_numShapes; s++)
				{
				    rtmp = 0;
				    for (q = 0; q < p.numQuadTotal; q++)
					rtmp += vol * S[INDEX2(s, q, p.row_numShapes)];
				    for (k = 0; k < p.numEqu; k++)
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp * Y_p[k];
				}
			    }
			}
		       /*********************************************************************************************************************/
			/* add the element matrices onto the matrix and right hand side                                */
		       /*********************************************************************************************************************/

			for (q = 0; q < p.row_numShapes; q++)
			    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];

			if (add_EM_F)
			    Dudley_Util_AddScatter(p.row_numShapes, row_index, p.numEqu, EM_F, F_p,
						   p.row_DOF_UpperBound);
			if (add_EM_S)
			    Dudley_Assemble_addToSystemMatrix(Mat, p.row_numShapes, row_index, p.numEqu,
							      p.row_numShapes, row_index, p.numComp, EM_S);

		    }		/* end color check */
		}		/* end element loop */
	    }			/* end color loop */

	    delete[] EM_S;
	    delete[] EM_F;
	    delete[] row_index;

	}			/* end of pointer check */
    }				/* end parallel region */
}

/*
 * $Log$
 */
