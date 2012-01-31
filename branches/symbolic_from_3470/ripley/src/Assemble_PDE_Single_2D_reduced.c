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
/*                                                                                                            */
/*    assembles single PDEs into the stiffness matrix S right hand side F                                     */
/*    for coefficients on the reduced integration scheme                                                      */
/*                                                                                                            */
/*      -(A_{i,j} u=,j)_i-(B_{i} u)_i+C_{j} u,j-D u  and -(X_{,i})_i + Y                                      */
/*                                                                                                            */
/*    u has p.numComp components in a 2D domain. The shape functions for test and solution must be identical  */
/*                                                                                                            */
/*    Shape of the coefficients:                                                                              */
/*      A = DIM x  DIM                                                                                        */
/*      B =  DIM                                                                                              */
/*      C =  DIM                                                                                              */
/*      D = scalar                                                                                            */
/*      X = DIM                                                                                               */
/*      Y = scalar                                                                                            */
/*                                                                                                            */
/**************************************************************/
#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/**************************************************************/
void  Finley_Assemble_PDE_Single_2D_reduced(Finley_Assemble_Parameters p,
Finley_ElementFile* elements,
Paso_SystemMatrix* Mat, escriptDataC* F,
escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y) {
   #define DIM 2
   index_t color;
   dim_t e, isub;
   __const double  *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
   double *EM_S, *EM_F, *Vol, *DSDX;
   index_t *row_index;
   register dim_t q, s,r,k,m;
   register double rtmp, rtmp0, rtmp1, rtmp2, rtmp00, rtmp01, rtmp02, rtmp10, rtmp11, rtmp12, rtmp20, rtmp21, rtmp22;
   bool_t add_EM_F, add_EM_S;
   bool_t extendedA=isExpanded(A);
   bool_t extendedB=isExpanded(B);
   bool_t extendedC=isExpanded(C);
   bool_t extendedD=isExpanded(D);
   bool_t extendedX=isExpanded(X);
   bool_t extendedY=isExpanded(Y);
   double *F_p=(requireWrite(F), getSampleDataRW(F,0));	/* use comma, to get around the mixed code and declarations thing */
   dim_t len_EM_S=p.row_numShapesTotal*p.col_numShapesTotal*p.numEqu*p.numComp;
   dim_t len_EM_F=p.row_numShapesTotal*p.numEqu;
   {
      /* GENERATOR SNIP_PRE TOP */
      const double w12 = -0.5*h0;
      const double w2 = -0.25000000000000000000;
      const double w1 = 0.25000000000000000000;
      const double w0 = -0.25*h1/h0;
      const double w15 = 0.25*h0*h1;
      const double w7 = -0.125*h0;
      const double w9 = 0.125*h0;
      const double w14 = 0.5*h0;
      const double w11 = -0.5*h1;
      const double w4 = -0.25*h0/h1;
      const double w10 = 0.0625*h0*h1;
      const double w13 = 0.5*h1;
      const double w8 = 0.125*h1;
      const double w6 = -0.125*h1;
      const double w3 = 0.25*h0/h1;
      const double w5 = 0.25*h1/h0;
      /* GENERATOR SNIP_PRE BOTTOM */
      #pragma omp parallel private(EM_S, EM_F,k2_0, k0, k1, k2, *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p)
      {
         EM_S=THREAD_MEMALLOC(len_EM_S,double);
         EM_F=THREAD_MEMALLOC(len_EM_F,double);
         if (!Finley_checkPtr(EM_S) && !Finley_checkPtr(EM_F) && !Finley_checkPtr(row_index) ) {
            for (k2_0 = 0; k2_0 <2; k2_0++) { /* coloring */
               #pragma omp parallel for private(i2, i1,i0)
               for (k2 = k2_0; k2< N0; k2=k2+2) {
                  for (k1 = 0; k1< N1; ++k1) {
                     for (k0 = 0; k0< N0; ++k0)  {
                        bool_t add_EM_F=FALSE;
                        bool_t add_EM_S=FALSE;
                        index_t e = k0 + M0 * k1 + M0*M1 * k2;
                        A_p=getSampleDataRO(A,e);
                        B_p=getSampleDataRO(B,e);
                        C_p=getSampleDataRO(C,e);
                        D_p=getSampleDataRO(D,e);
                        X_p=getSampleDataRO(X,e);
                        Y_p=getSampleDataRO(Y,e);
                        for (q=0;q<len_EM_S;++q) EM_S[q]=0;
                        for (q=0;q<len_EM_F;++q) EM_F[q]=0;
                        /* GENERATOR SNIP TOP */
                        /**************************************************************/
                        /*   process A: */
                        /**************************************************************/
                        if (NULL!=A_p) {
                           add_EM_S=TRUE;
                           const register double A_00 = A_p[INDEX2(0,0,2)];
                           const register double A_01 = A_p[INDEX2(0,1,2)];
                           const register double A_10 = A_p[INDEX2(1,0,2)];
                           const register double A_11 = A_p[INDEX2(1,1,2)];
                           const register double tmp0_0 = A_01 + A_10;
                           const register double tmp2_1 = A_01*w1;
                           const register double tmp7_1 = tmp0_0*w2;
                           const register double tmp3_1 = A_10*w2;
                           const register double tmp6_1 = A_00*w5;
                           const register double tmp1_1 = A_00*w0;
                           const register double tmp4_1 = tmp0_0*w1;
                           const register double tmp9_1 = A_01*w2;
                           const register double tmp5_1 = A_11*w4;
                           const register double tmp8_1 = A_10*w1;
                           const register double tmp0_1 = A_11*w3;
                           EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                           EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                           EM_S[INDEX2(3,2,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                           EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
                           EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp4_1 + tmp6_1;
                           EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                           EM_S[INDEX2(3,1,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                           EM_S[INDEX2(2,1,4)]+=tmp1_1 + tmp4_1 + tmp5_1;
                           EM_S[INDEX2(0,2,4)]+=tmp5_1 + tmp6_1 + tmp8_1 + tmp9_1;
                           EM_S[INDEX2(2,0,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                           EM_S[INDEX2(1,3,4)]+=tmp2_1 + tmp3_1 + tmp5_1 + tmp6_1;
                           EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                           EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                           EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1 + tmp8_1 + tmp9_1;
                           EM_S[INDEX2(0,3,4)]+=tmp1_1 + tmp5_1 + tmp7_1;
                           EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp6_1 + tmp7_1;
                        }
                        /**************************************************************/
                        /*   process B: */
                        /**************************************************************/
                        if (NULL!=B_p) {
                           add_EM_S=TRUE;
                           const register double B_0 = B_p[0];
                           const register double B_1 = B_p[1];
                           const register double tmp2_1 = B_0*w8;
                           const register double tmp0_1 = B_0*w6;
                           const register double tmp3_1 = B_1*w9;
                           const register double tmp1_1 = B_1*w7;
                           EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp2_1;
                           EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(3,3,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(3,0,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(3,1,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(2,1,4)]+=tmp0_1 + tmp3_1;
                           EM_S[INDEX2(0,2,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(2,0,4)]+=tmp0_1 + tmp3_1;
                           EM_S[INDEX2(1,3,4)]+=tmp1_1 + tmp2_1;
                           EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp3_1;
                           EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp3_1;
                           EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp2_1;
                           EM_S[INDEX2(0,3,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(1,1,4)]+=tmp1_1 + tmp2_1;
                        }
                        /**************************************************************/
                        /*   process C: */
                        /**************************************************************/
                        if (NULL!=C_p) {
                           add_EM_S=TRUE;
                           const register double C_0 = C_p[0];
                           const register double C_1 = C_p[1];
                           const register double tmp1_1 = C_1*w7;
                           const register double tmp0_1 = C_0*w8;
                           const register double tmp3_1 = C_0*w6;
                           const register double tmp2_1 = C_1*w9;
                           EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(1,2,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(0,0,4)]+=tmp1_1 + tmp3_1;
                           EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp2_1;
                           EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp3_1;
                           EM_S[INDEX2(3,1,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(2,1,4)]+=tmp0_1 + tmp1_1;
                           EM_S[INDEX2(0,2,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(2,0,4)]+=tmp1_1 + tmp3_1;
                           EM_S[INDEX2(1,3,4)]+=tmp0_1 + tmp2_1;
                           EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp2_1;
                           EM_S[INDEX2(2,2,4)]+=tmp2_1 + tmp3_1;
                           EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp3_1;
                           EM_S[INDEX2(0,3,4)]+=tmp0_1 + tmp2_1;
                           EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp1_1;
                        }
                        /**************************************************************/
                        /*   process D: */
                        /**************************************************************/
                        if (NULL!=D_p) {
                           add_EM_S=TRUE;
                           const register double D_0 = D_p[0];
                           const register double tmp0_1 = D_0*w10;
                           EM_S[INDEX2(0,1,4)]+=tmp0_1;
                           EM_S[INDEX2(1,2,4)]+=tmp0_1;
                           EM_S[INDEX2(3,2,4)]+=tmp0_1;
                           EM_S[INDEX2(0,0,4)]+=tmp0_1;
                           EM_S[INDEX2(3,3,4)]+=tmp0_1;
                           EM_S[INDEX2(3,0,4)]+=tmp0_1;
                           EM_S[INDEX2(3,1,4)]+=tmp0_1;
                           EM_S[INDEX2(2,1,4)]+=tmp0_1;
                           EM_S[INDEX2(0,2,4)]+=tmp0_1;
                           EM_S[INDEX2(2,0,4)]+=tmp0_1;
                           EM_S[INDEX2(1,3,4)]+=tmp0_1;
                           EM_S[INDEX2(2,3,4)]+=tmp0_1;
                           EM_S[INDEX2(2,2,4)]+=tmp0_1;
                           EM_S[INDEX2(1,0,4)]+=tmp0_1;
                           EM_S[INDEX2(0,3,4)]+=tmp0_1;
                           EM_S[INDEX2(1,1,4)]+=tmp0_1;
                        }
                        /**************************************************************/
                        /*   process X: */
                        /**************************************************************/
                        if (NULL!=X_p) {
                           add_EM_F=TRUE;
                           const register double X_0 = X_p[0];
                           const register double X_1 = X_p[1];
                           const register double tmp0_1 = X_0*w11;
                           const register double tmp2_1 = X_0*w13;
                           const register double tmp1_1 = X_1*w12;
                           const register double tmp3_1 = X_1*w14;
                           EM_F[0]+=tmp0_1 + tmp1_1;
                           EM_F[1]+=tmp1_1 + tmp2_1;
                           EM_F[2]+=tmp0_1 + tmp3_1;
                           EM_F[3]+=tmp2_1 + tmp3_1;
                        }
                        /**************************************************************/
                        /*   process Y: */
                        /**************************************************************/
                        if (NULL!=Y_p) {
                           add_EM_F=TRUE;
                           const register double Y_0 = Y_p[0];
                           const register double tmp0_1 = Y_0*w15;
                           EM_F[0]+=tmp0_1;
                           EM_F[1]+=tmp0_1;
                           EM_F[2]+=tmp0_1;
                           EM_F[3]+=tmp0_1;
                        }
                        /* GENERATOR SNIP BOTTOM */
                        /* add to matrix (if add_EM_S) and RHS (if add_EM_F)*/
                     } /* end k0 */
                  } /* end k1  */
               } /* end k2 */
            } /* end coloring */
         } /* end of pointer check */
         THREAD_MEMFREE(EM_S);
         THREAD_MEMFREE(EM_F);
      } /* end parallel region */
   } /* end of presettings */
}
