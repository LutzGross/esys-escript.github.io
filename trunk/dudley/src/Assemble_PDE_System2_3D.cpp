
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

/************************************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F  */
/*    the shape functions for test and solution must be identical */

/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  and -(X_{k,i})_i + Y_k */

/*    u has p.numComp components in a 3D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                                                  */

/*    Shape of the coefficients: */

/*      A = p.numEqu x 3 x p.numComp x 3 */
/*      B = 3 x p.numEqu x p.numComp  */
/*      C = p.numEqu x 3 x p.numComp  */
/*      D = p.numEqu x p.numComp  */
/*      X = p.numEqu x 3  */
/*      Y = p.numEqu   */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Dudley_Assemble_PDE_System2_3D(Dudley_Assemble_Parameters p, Dudley_ElementFile * elements,
				    escript::ASM_ptr mat, escript::Data* F,
				    const escript::Data* A, const escript::Data* B, const escript::Data* C, const escript::Data* D,
				    const escript::Data* X, const escript::Data* Y)
{

#define DIM 3
    index_t color;
    dim_t e;
    __const double *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
    double *EM_S, *EM_F, *DSDX;
    index_t *row_index;
    register dim_t q, s, r, k, m;
    register double rtmp, rtmp0, rtmp1, rtmp2, rtmp00, rtmp01, rtmp02, rtmp10, rtmp11, rtmp12, rtmp20, rtmp21, rtmp22;
    bool add_EM_F, add_EM_S;

    bool extendedA = isExpanded(A);
    bool extendedB = isExpanded(B);
    bool extendedC = isExpanded(C);
    bool extendedD = isExpanded(D);
    bool extendedX = isExpanded(X);
    bool extendedY = isExpanded(Y);
    double *F_p = (requireWrite(F), getSampleDataRW(F, 0));	/* use comma, to get around the mixed code and declarations thing */
    const double *S = p.shapeFns;
    dim_t len_EM_S = p.numShapes * p.numShapes * p.numEqu * p.numComp;
    dim_t len_EM_F = p.numShapes * p.numEqu;

#pragma omp parallel private(color,EM_S, EM_F, DSDX, A_p, B_p, C_p, D_p, X_p, Y_p, A_q, B_q, C_q, D_q, X_q, Y_q,row_index,q, s,r,k,m,rtmp, rtmp0, rtmp1, rtmp2, rtmp00, rtmp01, rtmp02, rtmp10, rtmp11, rtmp12, rtmp20, rtmp21, rtmp22,add_EM_F, add_EM_S)
    {
	EM_S = new  double[len_EM_S];
	EM_F = new  double[len_EM_F];
	row_index = new  index_t[p.numShapes];

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
			double vol = p.row_jac->absD[e] * p.row_jac->quadweight;

			A_p = getSampleDataRO(A, e);
			B_p = getSampleDataRO(B, e);
			C_p = getSampleDataRO(C, e);
			D_p = getSampleDataRO(D, e);
			X_p = getSampleDataRO(X, e);
			Y_p = getSampleDataRO(Y, e);


			DSDX = &(p.row_jac->DSDX[INDEX5(0, 0, 0, 0, e, p.numShapes, DIM, p.numQuad, 1)]);
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
				A_q = &(A_p[INDEX6(0, 0, 0, 0, 0, 0, p.numEqu, DIM, p.numComp, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						{
						    rtmp +=
							vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 0, m, 0, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 0, m, 1, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 0, m, 2, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 1, m, 0, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 1, m, 1, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 1, m, 2, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 2, m, 0, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 2, m, 1, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
							       DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
							       A_q[INDEX5(k, 2, m, 2, q, p.numEqu, DIM, p.numComp, DIM)]
							       * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)]);

						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp00 = 0;
					rtmp01 = 0;
					rtmp02 = 0;
					rtmp10 = 0;
					rtmp11 = 0;
					rtmp12 = 0;
					rtmp20 = 0;
					rtmp21 = 0;
					rtmp22 = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp0 = vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					    rtmp00 += rtmp0 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp01 += rtmp0 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					    rtmp02 += rtmp0 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];

					    rtmp1 = vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
					    rtmp10 += rtmp1 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp11 += rtmp1 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					    rtmp12 += rtmp1 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];

					    rtmp2 = vol * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
					    rtmp20 += rtmp2 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp21 += rtmp2 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					    rtmp22 += rtmp2 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];
					}
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] +=
						    rtmp00 * A_p[INDEX4(k, 0, m, 0, p.numEqu, DIM, p.numComp)] +
						    rtmp01 * A_p[INDEX4(k, 0, m, 1, p.numEqu, DIM, p.numComp)] +
						    rtmp02 * A_p[INDEX4(k, 0, m, 2, p.numEqu, DIM, p.numComp)] +
						    rtmp10 * A_p[INDEX4(k, 1, m, 0, p.numEqu, DIM, p.numComp)] +
						    rtmp11 * A_p[INDEX4(k, 1, m, 1, p.numEqu, DIM, p.numComp)] +
						    rtmp12 * A_p[INDEX4(k, 1, m, 2, p.numEqu, DIM, p.numComp)] +
						    rtmp20 * A_p[INDEX4(k, 2, m, 0, p.numEqu, DIM, p.numComp)] +
						    rtmp21 * A_p[INDEX4(k, 2, m, 1, p.numEqu, DIM, p.numComp)] +
						    rtmp22 * A_p[INDEX4(k, 2, m, 2, p.numEqu, DIM, p.numComp)];
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
				B_q = &(B_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, DIM, p.numComp, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						{
						    rtmp +=
							vol * S[INDEX2(r, q, p.numShapes)] *
							(DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
							 B_q[INDEX4(k, 0, m, q, p.numEqu, DIM, p.numComp)] +
							 DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
							 B_q[INDEX4(k, 1, m, q, p.numEqu, DIM, p.numComp)] +
							 DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
							 B_q[INDEX4(k, 2, m, q, p.numEqu, DIM, p.numComp)]);
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp0 = 0;
					rtmp1 = 0;
					rtmp2 = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp = vol * S[INDEX2(r, q, p.numShapes)];
					    rtmp0 += rtmp * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					    rtmp1 += rtmp * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
					    rtmp2 += rtmp * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
					}
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] +=
						    rtmp0 * B_p[INDEX3(k, 0, m, p.numEqu, DIM)] +
						    rtmp1 * B_p[INDEX3(k, 1, m, p.numEqu, DIM)] +
						    rtmp2 * B_p[INDEX3(k, 2, m, p.numEqu, DIM)];
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
				C_q = &(C_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, p.numComp, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						{
						    rtmp +=
							vol * S[INDEX2(s, q, p.numShapes)] *
							(C_q[INDEX4(k, m, 0, q, p.numEqu, p.numComp, DIM)] *
							 DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
							 C_q[INDEX4(k, m, 1, q, p.numEqu, p.numComp, DIM)] *
							 DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
							 C_q[INDEX4(k, m, 2, q, p.numEqu, p.numComp, DIM)] *
							 DSDX[INDEX3(r, 2, q, p.numShapes, DIM)]);
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp0 = 0;
					rtmp1 = 0;
					rtmp2 = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp = vol * S[INDEX2(s, q, p.numShapes)];
					    rtmp0 += rtmp * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp1 += rtmp * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					    rtmp2 += rtmp * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];
					}
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] +=
						    rtmp0 * C_p[INDEX3(k, m, 0, p.numEqu, p.numComp)] +
						    rtmp1 * C_p[INDEX3(k, m, 1, p.numEqu, p.numComp)] +
						    rtmp2 * C_p[INDEX3(k, m, 2, p.numEqu, p.numComp)];
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
				D_q = &(D_p[INDEX4(0, 0, 0, 0, p.numEqu, p.numComp, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						{
						    rtmp +=
							vol * S[INDEX2(s, q, p.numShapes)] *
							D_q[INDEX3(k, m, q, p.numEqu, p.numComp)] *
							S[INDEX2(r, q, p.numShapes)];
						}
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
					    }
					}
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuad; q++)
					    rtmp += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(r, q, p.numShapes)];
					for (k = 0; k < p.numEqu; k++)
					{
					    for (m = 0; m < p.numComp; m++)
					    {
						EM_S[INDEX4(k, m, s, r, p.numEqu, p.numComp, p.numShapes)] +=
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
				X_q = &(X_p[INDEX4(0, 0, 0, 0, p.numEqu, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (k = 0; k < p.numEqu; k++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp +=
						vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
						       X_q[INDEX3(k, 0, q, p.numEqu, DIM)] +
						       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
						       X_q[INDEX3(k, 1, q, p.numEqu, DIM)] +
						       DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
						       X_q[INDEX3(k, 2, q, p.numEqu, DIM)]);
					}
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp;
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp0 = 0;
				    rtmp1 = 0;
				    rtmp2 = 0;
				    for (q = 0; q < p.numQuad; q++)
				    {
					rtmp0 += vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					rtmp1 += vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
					rtmp2 += vol * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
				    }
				    for (k = 0; k < p.numEqu; k++)
				    {
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp0 * X_p[INDEX2(k, 0, p.numEqu)]
					    + rtmp1 * X_p[INDEX2(k, 1, p.numEqu)] + rtmp2 * X_p[INDEX2(k, 2, p.numEqu)];
				    }
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
				Y_q = &(Y_p[INDEX3(0, 0, 0, p.numEqu, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (k = 0; k < p.numEqu; k++)
				    {
					rtmp = 0.;
					for (q = 0; q < p.numQuad; q++)
					    rtmp += vol * S[INDEX2(s, q, p.numShapes)] * Y_q[INDEX2(k, q, p.numEqu)];
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp;
				    }
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp = 0;
				    for (q = 0; q < p.numQuad; q++)
					rtmp += vol * S[INDEX2(s, q, p.numShapes)];
				    for (k = 0; k < p.numEqu; k++)
					EM_F[INDEX2(k, s, p.numEqu)] += rtmp * Y_p[k];
				}
			    }
			}

		       /*********************************************************************************************************************/
			/* add the element matrices onto the matrix and right hand side                                */
		       /*********************************************************************************************************************/
			for (q = 0; q < p.numShapes; q++)
			    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];

			if (add_EM_F)
			    Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_F, F_p, p.row_DOF_UpperBound);
			if (add_EM_S)
			    Dudley_Assemble_addToSystemMatrix(mat, p.numShapes, row_index, p.numEqu,
							      p.numShapes, row_index, p.numComp, EM_S);
		    }		/* end color check */
		}		/* end element loop */
	    }			/* end color loop */

	    delete[] EM_S;
	    delete[] EM_F;
	    delete[] row_index;

	}			/* end of pointer check */
    }				/* end parallel region */
}

