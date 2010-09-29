
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/**************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F  */
/*    the shape functions for test and solution must be identical */

/*      -(A_{i,j} u_,j)_i-(B_{i} u)_i+C_{j} u_,j-D u_m  and -(X_,i)_i + Y */

/*    in a 2D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                         */

/*    Shape of the coefficients: */

/*      A = 2 x 2 */
/*      B = 2   */
/*      C = 2   */
/*      D = scalar  */
/*      X = 2  */
/*      Y = scalar   */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ShapeTable.h"

/**************************************************************/

void Dudley_Assemble_PDE_Single2_2D(Assemble_Parameters p, Dudley_ElementFile * elements,
				    Paso_SystemMatrix * Mat, escriptDataC * F,
				    escriptDataC * A, escriptDataC * B, escriptDataC * C, escriptDataC * D,
				    escriptDataC * X, escriptDataC * Y)
{

#define DIM 2
    index_t color;
    dim_t e;
    __const double *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
    double *EM_S, *EM_F, *DSDX;
    index_t *row_index;
    register dim_t q, s, r;
    register double rtmp00, rtmp01, rtmp10, rtmp11, rtmp, rtmp0, rtmp1;
    bool_t add_EM_F, add_EM_S;

    bool_t extendedA = isExpanded(A);
    bool_t extendedB = isExpanded(B);
    bool_t extendedC = isExpanded(C);
    bool_t extendedD = isExpanded(D);
    bool_t extendedX = isExpanded(X);
    bool_t extendedY = isExpanded(Y);
    double *F_p = (requireWrite(F), getSampleDataRW(F, 0));	/* use comma, to get around the mixed code and declarations thing */
    const double *S = p.shapeFns;
    dim_t len_EM_S = p.numShapes * p.numShapes;
    dim_t len_EM_F = p.numShapes;

#pragma omp parallel private(color,EM_S, EM_F, Vol, DSDX, A_p, B_p, C_p, D_p, X_p, Y_p, A_q, B_q, C_q, D_q, X_q, Y_q,row_index,q, s,r,rtmp00, rtmp01, rtmp10, rtmp11, rtmp, rtmp0, rtmp1,add_EM_F, add_EM_S)
    {
	EM_S = THREAD_MEMALLOC(len_EM_S, double);
	EM_F = THREAD_MEMALLOC(len_EM_F, double);
	row_index = THREAD_MEMALLOC(p.numShapes, index_t);

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

			DSDX = &(p.row_jac->DSDX[INDEX5(0, 0, 0, 0, e, p.numShapes, DIM, p.numQuad, 1)]);
			for (q = 0; q < len_EM_S; ++q)
			    EM_S[q] = 0;
			for (q = 0; q < len_EM_F; ++q)
			    EM_F[q] = 0;
			add_EM_F = FALSE;
			add_EM_S = FALSE;
		     /**************************************************************/
			/*   process A: */
		     /**************************************************************/
			if (NULL != A_p)
			{
			    add_EM_S = TRUE;
			    if (extendedA)
			    {
				A_q = &(A_p[INDEX4(0, 0, 0, 0, DIM, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp +=
						vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
						       A_q[INDEX3(0, 0, q, DIM, DIM)] *
						       DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
						       DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
						       A_q[INDEX3(0, 1, q, DIM, DIM)] *
						       DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
						       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
						       A_q[INDEX3(1, 0, q, DIM, DIM)] *
						       DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
						       DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
						       A_q[INDEX3(1, 1, q, DIM, DIM)] *
						       DSDX[INDEX3(r, 1, q, p.numShapes, DIM)]);

					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
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
					rtmp10 = 0;
					rtmp11 = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp0 = vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					    rtmp1 = vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
					    rtmp00 += rtmp0 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp01 += rtmp0 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					    rtmp10 += rtmp1 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp11 += rtmp1 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] +=
					    rtmp00 * A_p[INDEX2(0, 0, DIM)] + rtmp01 * A_p[INDEX2(0, 1, DIM)] +
					    rtmp10 * A_p[INDEX2(1, 0, DIM)] + rtmp11 * A_p[INDEX2(1, 1, DIM)];
				    }
				}
			    }
			}
		     /**************************************************************/
			/*   process B: */
		     /**************************************************************/
			if (NULL != B_p)
			{
			    add_EM_S = TRUE;
			    if (extendedB)
			    {
				B_q = &(B_p[INDEX3(0, 0, 0, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp = 0.;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp +=
						vol * S[INDEX2(r, q, p.numShapes)] *
						(DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
						 B_q[INDEX2(0, q, DIM)] +
						 DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] * B_q[INDEX2(1, q, DIM)]);
					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
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
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp = vol * S[INDEX2(r, q, p.numShapes)];
					    rtmp0 += rtmp * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					    rtmp1 += rtmp * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] +=
					    rtmp0 * B_p[0] + rtmp1 * B_p[1];
				    }
				}
			    }
			}
		     /**************************************************************/
			/*   process C: */
		     /**************************************************************/
			if (NULL != C_p)
			{
			    add_EM_S = TRUE;
			    if (extendedC)
			    {
				C_q = &(C_p[INDEX3(0, 0, 0, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp +=
						vol * S[INDEX2(s, q, p.numShapes)] * (C_q[INDEX2(0, q, DIM)] *
										      DSDX[INDEX3
											   (r, 0, q,
											    p.numShapes,
											    DIM)] + C_q[INDEX2(1, q,
													       DIM)]
										      *
										      DSDX[INDEX3
											   (r, 1, q,
											    p.numShapes, DIM)]);
					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
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
					for (q = 0; q < p.numQuad; q++)
					{
					    rtmp = vol * S[INDEX2(s, q, p.numShapes)];
					    rtmp0 += rtmp * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
					    rtmp1 += rtmp * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
					}
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] +=
					    rtmp0 * C_p[0] + rtmp1 * C_p[1];
				    }
				}
			    }
			}
		     /************************************************************* */
			/* process D */
		     /**************************************************************/
			if (NULL != D_p)
			{
			    add_EM_S = TRUE;
			    if (extendedD)
			    {
				D_q = &(D_p[INDEX2(0, 0, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    for (r = 0; r < p.numShapes; r++)
				    {
					rtmp = 0;
					for (q = 0; q < p.numQuad; q++)
					    rtmp +=
						vol * S[INDEX2(s, q, p.numShapes)] * D_q[q] *
						S[INDEX2(r, q, p.numShapes)];
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp;
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
					EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numComp, p.numShapes)] += rtmp * D_p[0];
				    }
				}
			    }
			}
		     /**************************************************************/
			/*   process X: */
		     /**************************************************************/
			if (NULL != X_p)
			{
			    add_EM_F = TRUE;
			    if (extendedX)
			    {
				X_q = &(X_p[INDEX3(0, 0, 0, DIM, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp = 0.;
				    for (q = 0; q < p.numQuad; q++)
				    {
					rtmp +=
					    vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
						   X_q[INDEX2(0, q, DIM)] +
						   DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] * X_q[INDEX2(1, q, DIM)]);
				    }
				    EM_F[INDEX2(0, s, p.numEqu)] += rtmp;
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp0 = 0.;
				    rtmp1 = 0.;
				    for (q = 0; q < p.numQuad; q++)
				    {
					rtmp0 += vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
					rtmp1 += vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
				    }
				    EM_F[INDEX2(0, s, p.numEqu)] += rtmp0 * X_p[0] + rtmp1 * X_p[1];
				}
			    }
			}
		    /**************************************************************/
			/*   process Y: */
		    /**************************************************************/
			if (NULL != Y_p)
			{
			    add_EM_F = TRUE;
			    if (extendedY)
			    {
				Y_q = &(Y_p[INDEX2(0, 0, p.numQuad)]);
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp = 0;
				    for (q = 0; q < p.numQuad; q++)
					rtmp += vol * S[INDEX2(s, q, p.numShapes)] * Y_q[q];
				    EM_F[INDEX2(0, s, p.numEqu)] += rtmp;
				}
			    }
			    else
			    {
				for (s = 0; s < p.numShapes; s++)
				{
				    rtmp = 0;
				    for (q = 0; q < p.numQuad; q++)
					rtmp += vol * S[INDEX2(s, q, p.numShapes)];
				    EM_F[INDEX2(0, s, p.numEqu)] += rtmp * Y_p[0];
				}
			    }
			}
		      /***********************************************************************************************/
			/* add the element matrices onto the matrix and right hand side                                */
		      /***********************************************************************************************/

			for (q = 0; q < p.numShapes; q++)
			    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
			if (add_EM_F)
			    Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_F, F_p, p.row_DOF_UpperBound);
			if (add_EM_S)
			    Dudley_Assemble_addToSystemMatrix(Mat, p.numShapes, row_index, p.numEqu,
							      p.numShapes, row_index, p.numComp, EM_S);
		    }		/* end color check */
		}		/* end element loop */
	    }			/* end color loop */

	    THREAD_MEMFREE(EM_S);
	    THREAD_MEMFREE(EM_F);
	    THREAD_MEMFREE(row_index);

	}			/* end of pointer check */
    }				/* end parallel region */
}
