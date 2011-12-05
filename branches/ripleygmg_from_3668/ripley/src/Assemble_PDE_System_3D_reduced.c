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
/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F                        */
/*    for coefficients on the reduced integration scheme                                                      */
/*                                                                                                            */
/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  and -(X_{k,i})_i + Y_k          */
/*                                                                                                            */
/*    u has p.numComp components in a 2D domain. The shape functions for test and solution must be identical  */
/*                                                                                                            */
/*      A = p.numEqu x DIM x p.numComp x DIM                                                                  */
/*      B = p.numEqu  x p.numComp x DIM                                                                       */
/*      C = p.numEqu x DIM x p.numComp                                                                        */
/*      D = p.numEqu x p.numComp                                                                              */
/*      X = p.numEqu x DIM                                                                                    */
/*      Y = p.numEqu                                                                                          */
/*                                                                                                            */
/**************************************************************/
#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/**************************************************************/
void  Finley_Assemble_PDE_System_3D_reduced(Finley_Assemble_Parameters p,
Finley_ElementFile* elements,
Paso_SystemMatrix* Mat, escriptDataC* F,
escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y) {
   #define DIM 3
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
      const double w15 = -0.03125*h1*h2;
      const double w17 = -0.03125*h0*h1;
      const double w0 = 0.0625*h1*h2/h0;
      const double w20 = -0.25*h0*h2;
      const double w2 = -0.0625*h1;
      const double w18 = 0.015625*h0*h1*h2;
      const double w3 = 0.0625*h0*h2/h1;
      const double w5 = 0.0625*h1;
      const double w19 = -0.25*h1*h2;
      const double w10 = -0.0625*h0*h2/h1;
      const double w12 = 0.03125*h1*h2;
      const double w21 = -0.25*h0*h1;
      const double w6 = 0.0625*h0;
      const double w14 = 0.03125*h0*h1;
      const double w24 = 0.25*h0*h1;
      const double w4 = -0.0625*h0;
      const double w7 = -0.0625*h0*h1/h2;
      const double w1 = 0.0625*h2;
      const double w13 = 0.03125*h0*h2;
      const double w11 = 0.0625*h0*h1/h2;
      const double w23 = 0.25*h0*h2;
      const double w25 = 0.125*h0*h1*h2;
      const double w9 = -0.0625*h2;
      const double w16 = -0.03125*h0*h2;
      const double w8 = -0.0625*h1*h2/h0;
      const double w22 = 0.25*h1*h2;
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
                           for (k=0;k<p.numEqu;k++) {
                              for (m=0;m<p.numComp;m++) {
                                 const register double A_00 = A_p[INDEX4(k,0,m,0 p.numEqu,3, p.numComp)];
                                 const register double A_01 = A_p[INDEX4(k,0,m,1 p.numEqu,3, p.numComp)];
                                 const register double A_02 = A_p[INDEX4(k,0,m,2 p.numEqu,3, p.numComp)];
                                 const register double A_10 = A_p[INDEX4(k,1,m,0 p.numEqu,3, p.numComp)];
                                 const register double A_11 = A_p[INDEX4(k,1,m,1 p.numEqu,3, p.numComp)];
                                 const register double A_12 = A_p[INDEX4(k,1,m,2 p.numEqu,3, p.numComp)];
                                 const register double A_20 = A_p[INDEX4(k,2,m,0 p.numEqu,3, p.numComp)];
                                 const register double A_21 = A_p[INDEX4(k,2,m,1 p.numEqu,3, p.numComp)];
                                 const register double A_22 = A_p[INDEX4(k,2,m,2 p.numEqu,3, p.numComp)];
                                 const register double tmp0_0 = A_01 + A_10;
                                 const register double tmp1_0 = A_02 + A_20;
                                 const register double tmp2_0 = A_12 + A_21;
                                 const register double tmp3_1 = A_22*w7;
                                 const register double tmp10_1 = A_11*w10;
                                 const register double tmp21_1 = A_02*w5;
                                 const register double tmp2_1 = A_00*w0;
                                 const register double tmp23_1 = tmp2_0*w6;
                                 const register double tmp19_1 = A_20*w2;
                                 const register double tmp4_1 = A_11*w3;
                                 const register double tmp22_1 = tmp1_0*w5;
                                 const register double tmp13_1 = A_21*w4;
                                 const register double tmp5_1 = A_21*w6;
                                 const register double tmp16_1 = A_10*w9;
                                 const register double tmp7_1 = A_20*w5;
                                 const register double tmp18_1 = tmp2_0*w4;
                                 const register double tmp6_1 = A_02*w2;
                                 const register double tmp9_1 = A_22*w11;
                                 const register double tmp15_1 = tmp1_0*w2;
                                 const register double tmp12_1 = A_01*w1;
                                 const register double tmp0_1 = tmp0_0*w1;
                                 const register double tmp20_1 = A_01*w9;
                                 const register double tmp14_1 = A_12*w6;
                                 const register double tmp1_1 = A_12*w4;
                                 const register double tmp8_1 = A_00*w8;
                                 const register double tmp11_1 = tmp0_0*w9;
                                 const register double tmp17_1 = A_10*w1;
                                 EM_S[INDEX4(k,m,7,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,4,7,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,6,4,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp2_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp1_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,5,4,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,0,7,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,5,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,2,6,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,1,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,5,1,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp13_1 + tmp14_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,3,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,2,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,7,2,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,4,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,6,7,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,7,6,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,4,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp15_1 + tmp18_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,6,3,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,1,5,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,3,6,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,7,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,5,7,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp1_1 + tmp22_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,5,3,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,4,1,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp16_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,2,7,p.numEqu,p.numComp,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,6,6,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp15_1 + tmp23_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,5,0,p.numEqu,p.numComp,8)]+=tmp13_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp20_1 + tmp3_1 + tmp4_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,7,1,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,4,5,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp16_1 + tmp18_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,0,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,5,5,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp18_1 + tmp22_1 + tmp2_1 + tmp4_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,1,4,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp1_1 + tmp20_1 + tmp22_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,6,0,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,7,5,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp5_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,4,2,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,8)]+=tmp17_1 + tmp20_1 + tmp23_1 + tmp4_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,6,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp14_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,3,5,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp23_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp16_1 + tmp19_1 + tmp21_1 + tmp23_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,7,0,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp15_1 + tmp18_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,4,6,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,5,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp15_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,6,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp10_1 + tmp18_1 + tmp22_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp20_1 + tmp2_1 + tmp5_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp13_1 + tmp14_1 + tmp17_1 + tmp20_1 + tmp22_1 + tmp2_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,7,4,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp19_1 + tmp21_1 + tmp8_1 + tmp9_1;
                                 EM_S[INDEX4(k,m,0,6,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp17_1 + tmp18_1 + tmp20_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                                 EM_S[INDEX4(k,m,6,2,p.numEqu,p.numComp,8)]+=tmp11_1 + tmp19_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,3,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,1,7,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp18_1 + tmp19_1 + tmp21_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,0,5,p.numEqu,p.numComp,8)]+=tmp12_1 + tmp15_1 + tmp16_1 + tmp1_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,3,4,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp11_1 + tmp22_1 + tmp23_1 + tmp3_1 + tmp8_1;
                                 EM_S[INDEX4(k,m,2,4,p.numEqu,p.numComp,8)]+=tmp10_1 + tmp12_1 + tmp16_1 + tmp23_1 + tmp2_1 + tmp3_1 + tmp6_1 + tmp7_1;
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process B: */
                        /**************************************************************/
                        if (NULL!=B_p) {
                           add_EM_S=TRUE;
                           for (k=0;k<p.numEqu;k++) {
                              for (m=0;m<p.numComp;m++) {
                                 const register double B_0 = B_p[INDEX3(k,0,m, p.numEqu,3)];
                                 const register double B_1 = B_p[INDEX3(k,1,m, p.numEqu,3)];
                                 const register double B_2 = B_p[INDEX3(k,2,m, p.numEqu,3)];
                                 const register double tmp4_1 = B_0*w15;
                                 const register double tmp3_1 = B_1*w16;
                                 const register double tmp2_1 = B_0*w12;
                                 const register double tmp5_1 = B_2*w17;
                                 const register double tmp1_1 = B_2*w14;
                                 const register double tmp0_1 = B_1*w13;
                                 EM_S[INDEX4(k,m,7,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,7,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,4,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,0,7,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,6,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,2,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,1,6,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,3,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,4,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,6,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,3,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,5,7,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,5,3,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,4,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,5,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,7,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,5,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,0,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,5,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,1,4,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,7,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,2,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,5,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,0,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,6,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,5,2,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,6,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,0,6,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,4,3,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,7,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,5,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,3,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,4,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp4_1 + tmp5_1;
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process C: */
                        /**************************************************************/
                        if (NULL!=C_p) {
                           add_EM_S=TRUE;
                           for (k=0;k<p.numEqu;k++) {
                              for (m=0;m<p.numComp;m++) {
                                 const register double C_0 = C_p[INDEX3(k,m,0, p.numEqu,p.numComp)];
                                 const register double C_1 = C_p[INDEX3(k,m,1, p.numEqu,p.numComp)];
                                 const register double C_2 = C_p[INDEX3(k,m,2, p.numEqu,p.numComp)];
                                 const register double tmp5_1 = C_0*w15;
                                 const register double tmp2_1 = C_0*w12;
                                 const register double tmp4_1 = C_1*w16;
                                 const register double tmp1_1 = C_2*w17;
                                 const register double tmp3_1 = C_2*w14;
                                 const register double tmp0_1 = C_1*w13;
                                 EM_S[INDEX4(k,m,7,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,6,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,5,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,1,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,2,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,7,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,1,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,5,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,5,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,4,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,2,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,4,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,0,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,1,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,4,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,7,0,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,5,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,8)]+=tmp1_1 + tmp2_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,7,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,0,6,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp3_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,6,2,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,4,3,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp1_1 + tmp2_1;
                                 EM_S[INDEX4(k,m,1,7,p.numEqu,p.numComp,8)]+=tmp0_1 + tmp2_1 + tmp3_1;
                                 EM_S[INDEX4(k,m,0,5,p.numEqu,p.numComp,8)]+=tmp2_1 + tmp3_1 + tmp4_1;
                                 EM_S[INDEX4(k,m,3,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_S[INDEX4(k,m,2,4,p.numEqu,p.numComp,8)]+=tmp3_1 + tmp4_1 + tmp5_1;
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process D: */
                        /**************************************************************/
                        if (NULL!=D_p) {
                           add_EM_S=TRUE;
                           for (k=0;k<p.numEqu;k++) {
                              for (m=0;m<p.numComp;m++) {
                                 const register double D_0 = D_p[INDEX2(k,m, p.numEqu)];
                                 const register double tmp0_1 = D_0*w18;
                                 EM_S[INDEX4(k,m,7,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,0,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,5,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,7,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,6,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,6,2,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,4,3,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,1,7,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,0,5,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,3,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                                 EM_S[INDEX4(k,m,2,4,p.numEqu,p.numComp,8)]+=tmp0_1;
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process X: */
                        /**************************************************************/
                        if (NULL!=X_p) {
                           add_EM_F=TRUE;
                           for (k=0;k<p.numEqu;k++) {
                              const register double X_0 = X_p[INDEX2(k,0, p.numEqu)];
                              const register double X_1 = X_p[INDEX2(k,1, p.numEqu)];
                              const register double X_2 = X_p[INDEX2(k,2, p.numEqu)];
                              const register double tmp1_1 = X_0*w19;
                              const register double tmp2_1 = X_1*w20;
                              const register double tmp3_1 = X_0*w22;
                              const register double tmp4_1 = X_1*w23;
                              const register double tmp5_1 = X_2*w24;
                              const register double tmp0_1 = X_2*w21;
                              EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1 + tmp1_1 + tmp2_1;
                              EM_F[INDEX2(k,1,p.numEqu)]+=tmp0_1 + tmp2_1 + tmp3_1;
                              EM_F[INDEX2(k,2,p.numEqu)]+=tmp0_1 + tmp1_1 + tmp4_1;
                              EM_F[INDEX2(k,3,p.numEqu)]+=tmp0_1 + tmp3_1 + tmp4_1;
                              EM_F[INDEX2(k,4,p.numEqu)]+=tmp1_1 + tmp2_1 + tmp5_1;
                              EM_F[INDEX2(k,5,p.numEqu)]+=tmp2_1 + tmp3_1 + tmp5_1;
                              EM_F[INDEX2(k,6,p.numEqu)]+=tmp1_1 + tmp4_1 + tmp5_1;
                              EM_F[INDEX2(k,7,p.numEqu)]+=tmp3_1 + tmp4_1 + tmp5_1;
                           }
                        }
                        /**************************************************************/
                        /*   process Y: */
                        /**************************************************************/
                        if (NULL!=Y_p) {
                           add_EM_F=TRUE;
                           for (k=0;k<p.numEqu;k++) {
                              const register double Y_0 = Y_p[k];
                              const register double tmp0_1 = Y_0*w25;
                              EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,1,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,2,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,3,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,4,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,5,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,6,p.numEqu)]+=tmp0_1;
                              EM_F[INDEX2(k,7,p.numEqu)]+=tmp0_1;
                           }
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
