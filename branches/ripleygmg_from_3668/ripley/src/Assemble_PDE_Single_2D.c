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
void  Finley_Assemble_PDE_Single_2D(Finley_Assemble_Parameters p,
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
      const double w42 = 0.0023593469594139828636*h0;
      const double w29 = -0.032861463941450536761*h0;
      const double w44 = -0.16666666666666666667*h1;
      const double w67 = 0.052831216351296779436*h0;
      const double w9 = -0.041666666666666666667*h1/h0;
      const double w0 = -0.1555021169820365539*h1/h0;
      const double w23 = -0.16666666666666666667*h0/h1;
      const double w75 = 0.25*h0*h1;
      const double w71 = 0.5*h0;
      const double w21 = 0.16666666666666666667*h0/h1;
      const double w68 = -0.5*h1;
      const double w65 = 0.052831216351296779436*h1;
      const double w35 = 0.008805202725216129906*h1;
      const double w5 = -0.041666666666666666667;
      const double w55 = 0.09672363354357992482*h0*h1;
      const double w11 = 0.1555021169820365539*h1/h0;
      const double w39 = 0.032861463941450536761*h0;
      const double w66 = 0.19716878364870322056*h0;
      const double w28 = -0.032861463941450536761*h1;
      const double w62 = -0.052831216351296779436*h0;
      const double w69 = -0.5*h0;
      const double w38 = 0.12264065304058601714*h1;
      const double w3 = 0.041666666666666666667*h0/h1;
      const double w48 = 0.083333333333333333333*h0;
      const double w40 = -0.12264065304058601714*h0;
      const double w73 = 0.041666666666666666667*h0*h1;
      const double w70 = 0.5*h1;
      const double w4 = 0.15550211698203655390;
      const double w25 = 0.33333333333333333333*h0/h1;
      const double w47 = 0.16666666666666666667*h1;
      const double w30 = -0.12264065304058601714*h1;
      const double w1 = 0.041666666666666666667;
      const double w41 = -0.0023593469594139828636*h0;
      const double w59 = 0.11111111111111111111*h0*h1;
      const double w63 = -0.052831216351296779436*h1;
      const double w49 = -0.16666666666666666667*h0;
      const double w34 = 0.032861463941450536761*h1;
      const double w20 = -0.25000000000000000000;
      const double w53 = 0.0018607582807716854616*h0*h1;
      const double w72 = 0.1555021169820365539*h0*h1;
      const double w13 = 0.01116454968463011277*h0/h1;
      const double w54 = 0.0069444444444444444444*h0*h1;
      const double w2 = -0.15550211698203655390;
      const double w58 = 0.027777777777777777778*h0*h1;
      const double w22 = -0.16666666666666666667*h1/h0;
      const double w16 = -0.01116454968463011277*h0/h1;
      const double w7 = 0.011164549684630112770;
      const double w36 = 0.008805202725216129906*h0;
      const double w33 = -0.008805202725216129906*h1;
      const double w64 = 0.19716878364870322056*h1;
      const double w26 = 0.16666666666666666667*h1/h0;
      const double w46 = 0.083333333333333333333*h1;
      const double w74 = 0.01116454968463011277*h0*h1;
      const double w43 = 0.12264065304058601714*h0;
      const double w31 = -0.0023593469594139828636*h1;
      const double w8 = -0.011164549684630112770;
      const double w45 = -0.083333333333333333333*h0;
      const double w17 = -0.1555021169820365539*h0/h1;
      const double w50 = 0.16666666666666666667*h0;
      const double w19 = 0.25000000000000000000;
      const double w57 = 0.055555555555555555556*h0*h1;
      const double w24 = 0.33333333333333333333*h1/h0;
      const double w18 = -0.33333333333333333333*h1/h0;
      const double w12 = 0.1555021169820365539*h0/h1;
      const double w61 = -0.19716878364870322056*h0;
      const double w37 = 0.0023593469594139828636*h1;
      const double w14 = 0.01116454968463011277*h1/h0;
      const double w6 = -0.01116454968463011277*h1/h0;
      const double w27 = -0.33333333333333333333*h0/h1;
      const double w32 = -0.008805202725216129906*h0;
      const double w56 = 0.00049858867864229740201*h0*h1;
      const double w52 = 0.025917019497006092316*h0*h1;
      const double w60 = -0.19716878364870322056*h1;
      const double w51 = -0.083333333333333333333*h1;
      const double w15 = 0.041666666666666666667*h1/h0;
      const double w10 = -0.041666666666666666667*h0/h1;
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
                           if (extendedA) {
                              const register double A_00_0 = A_p[INDEX3(0,0,0,2,2)];
                              const register double A_01_0 = A_p[INDEX3(0,1,0,2,2)];
                              const register double A_10_0 = A_p[INDEX3(1,0,0,2,2)];
                              const register double A_11_0 = A_p[INDEX3(1,1,0,2,2)];
                              const register double A_00_1 = A_p[INDEX3(0,0,1,2,2)];
                              const register double A_01_1 = A_p[INDEX3(0,1,1,2,2)];
                              const register double A_10_1 = A_p[INDEX3(1,0,1,2,2)];
                              const register double A_11_1 = A_p[INDEX3(1,1,1,2,2)];
                              const register double A_00_2 = A_p[INDEX3(0,0,2,2,2)];
                              const register double A_01_2 = A_p[INDEX3(0,1,2,2,2)];
                              const register double A_10_2 = A_p[INDEX3(1,0,2,2,2)];
                              const register double A_11_2 = A_p[INDEX3(1,1,2,2,2)];
                              const register double A_00_3 = A_p[INDEX3(0,0,3,2,2)];
                              const register double A_01_3 = A_p[INDEX3(0,1,3,2,2)];
                              const register double A_10_3 = A_p[INDEX3(1,0,3,2,2)];
                              const register double A_11_3 = A_p[INDEX3(1,1,3,2,2)];
                              const register double tmp4_0 = A_10_1 + A_10_2;
                              const register double tmp12_0 = A_11_0 + A_11_2;
                              const register double tmp18_0 = A_01_1 + A_10_1;
                              const register double tmp2_0 = A_11_0 + A_11_1 + A_11_2 + A_11_3;
                              const register double tmp6_0 = A_01_3 + A_10_0;
                              const register double tmp10_0 = A_01_3 + A_10_3;
                              const register double tmp0_0 = A_01_0 + A_01_3;
                              const register double tmp3_0 = A_00_2 + A_00_3;
                              const register double tmp11_0 = A_11_1 + A_11_3;
                              const register double tmp14_0 = A_01_0 + A_01_3 + A_10_0 + A_10_3;
                              const register double tmp1_0 = A_00_0 + A_00_1;
                              const register double tmp8_0 = A_01_1 + A_01_2 + A_10_1 + A_10_2;
                              const register double tmp7_0 = A_01_0 + A_10_3;
                              const register double tmp5_0 = A_00_0 + A_00_1 + A_00_2 + A_00_3;
                              const register double tmp19_0 = A_01_2 + A_10_2;
                              const register double tmp13_0 = A_01_2 + A_10_1;
                              const register double tmp17_0 = A_01_1 + A_01_2;
                              const register double tmp15_0 = A_01_1 + A_10_2;
                              const register double tmp16_0 = A_10_0 + A_10_3;
                              const register double tmp9_0 = A_01_0 + A_10_0;
                              const register double tmp14_1 = A_10_0*w8;
                              const register double tmp23_1 = tmp3_0*w14;
                              const register double tmp35_1 = A_01_0*w8;
                              const register double tmp54_1 = tmp13_0*w8;
                              const register double tmp20_1 = tmp9_0*w4;
                              const register double tmp25_1 = tmp12_0*w12;
                              const register double tmp2_1 = A_01_1*w4;
                              const register double tmp44_1 = tmp7_0*w7;
                              const register double tmp26_1 = tmp10_0*w4;
                              const register double tmp52_1 = tmp18_0*w8;
                              const register double tmp48_1 = A_10_1*w7;
                              const register double tmp46_1 = A_01_3*w8;
                              const register double tmp50_1 = A_01_0*w2;
                              const register double tmp8_1 = tmp4_0*w5;
                              const register double tmp56_1 = tmp19_0*w8;
                              const register double tmp9_1 = tmp2_0*w10;
                              const register double tmp19_1 = A_10_3*w2;
                              const register double tmp47_1 = A_10_2*w4;
                              const register double tmp16_1 = tmp3_0*w0;
                              const register double tmp18_1 = tmp1_0*w6;
                              const register double tmp31_1 = tmp11_0*w12;
                              const register double tmp55_1 = tmp15_0*w2;
                              const register double tmp39_1 = A_10_2*w7;
                              const register double tmp11_1 = tmp6_0*w7;
                              const register double tmp40_1 = tmp11_0*w17;
                              const register double tmp34_1 = tmp15_0*w8;
                              const register double tmp33_1 = tmp14_0*w5;
                              const register double tmp24_1 = tmp11_0*w13;
                              const register double tmp3_1 = tmp1_0*w0;
                              const register double tmp5_1 = tmp2_0*w3;
                              const register double tmp43_1 = tmp17_0*w5;
                              const register double tmp15_1 = A_01_2*w4;
                              const register double tmp53_1 = tmp19_0*w2;
                              const register double tmp27_1 = tmp3_0*w11;
                              const register double tmp32_1 = tmp13_0*w2;
                              const register double tmp10_1 = tmp5_0*w9;
                              const register double tmp37_1 = A_10_1*w4;
                              const register double tmp38_1 = tmp5_0*w15;
                              const register double tmp17_1 = A_01_1*w7;
                              const register double tmp12_1 = tmp7_0*w4;
                              const register double tmp22_1 = tmp10_0*w7;
                              const register double tmp57_1 = tmp18_0*w2;
                              const register double tmp28_1 = tmp9_0*w7;
                              const register double tmp29_1 = tmp1_0*w14;
                              const register double tmp51_1 = tmp11_0*w16;
                              const register double tmp42_1 = tmp12_0*w16;
                              const register double tmp49_1 = tmp12_0*w17;
                              const register double tmp21_1 = tmp1_0*w11;
                              const register double tmp1_1 = tmp0_0*w1;
                              const register double tmp45_1 = tmp6_0*w4;
                              const register double tmp7_1 = A_10_0*w2;
                              const register double tmp6_1 = tmp3_0*w6;
                              const register double tmp13_1 = tmp8_0*w1;
                              const register double tmp36_1 = tmp16_0*w1;
                              const register double tmp41_1 = A_01_3*w2;
                              const register double tmp30_1 = tmp12_0*w13;
                              const register double tmp4_1 = A_01_2*w7;
                              const register double tmp0_1 = A_10_3*w8;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1;
                              EM_S[INDEX2(1,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp9_1;
                              EM_S[INDEX2(3,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp5_1 + tmp8_1;
                              EM_S[INDEX2(0,0,4)]+=tmp13_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1;
                              EM_S[INDEX2(3,3,4)]+=tmp13_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                              EM_S[INDEX2(3,0,4)]+=tmp10_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp9_1;
                              EM_S[INDEX2(3,1,4)]+=tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1;
                              EM_S[INDEX2(2,1,4)]+=tmp10_1 + tmp13_1 + tmp44_1 + tmp45_1 + tmp9_1;
                              EM_S[INDEX2(0,2,4)]+=tmp36_1 + tmp38_1 + tmp43_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                              EM_S[INDEX2(2,0,4)]+=tmp0_1 + tmp15_1 + tmp17_1 + tmp1_1 + tmp38_1 + tmp49_1 + tmp51_1 + tmp7_1 + tmp8_1;
                              EM_S[INDEX2(1,3,4)]+=tmp14_1 + tmp19_1 + tmp1_1 + tmp2_1 + tmp38_1 + tmp40_1 + tmp42_1 + tmp4_1 + tmp8_1;
                              EM_S[INDEX2(2,3,4)]+=tmp16_1 + tmp18_1 + tmp35_1 + tmp36_1 + tmp41_1 + tmp43_1 + tmp47_1 + tmp48_1 + tmp5_1;
                              EM_S[INDEX2(2,2,4)]+=tmp24_1 + tmp25_1 + tmp27_1 + tmp29_1 + tmp33_1 + tmp52_1 + tmp53_1;
                              EM_S[INDEX2(1,0,4)]+=tmp36_1 + tmp37_1 + tmp39_1 + tmp3_1 + tmp43_1 + tmp46_1 + tmp50_1 + tmp5_1 + tmp6_1;
                              EM_S[INDEX2(0,3,4)]+=tmp10_1 + tmp33_1 + tmp54_1 + tmp55_1 + tmp9_1;
                              EM_S[INDEX2(1,1,4)]+=tmp21_1 + tmp23_1 + tmp30_1 + tmp31_1 + tmp33_1 + tmp56_1 + tmp57_1;
                           } else { /* constant data */
                              const register double A_00 = A_p[INDEX2(0,0,2)];
                              const register double A_01 = A_p[INDEX2(0,1,2)];
                              const register double A_10 = A_p[INDEX2(1,0,2)];
                              const register double A_11 = A_p[INDEX2(1,1,2)];
                              const register double tmp0_0 = A_01 + A_10;
                              const register double tmp0_1 = A_00*w18;
                              const register double tmp10_1 = A_01*w20;
                              const register double tmp12_1 = A_00*w26;
                              const register double tmp4_1 = A_00*w22;
                              const register double tmp8_1 = A_00*w24;
                              const register double tmp13_1 = A_10*w19;
                              const register double tmp9_1 = tmp0_0*w20;
                              const register double tmp3_1 = A_11*w21;
                              const register double tmp11_1 = A_11*w27;
                              const register double tmp1_1 = A_01*w19;
                              const register double tmp6_1 = A_11*w23;
                              const register double tmp7_1 = A_11*w25;
                              const register double tmp2_1 = A_10*w20;
                              const register double tmp5_1 = tmp0_0*w19;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                              EM_S[INDEX2(1,2,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                              EM_S[INDEX2(3,2,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                              EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                              EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp7_1 + tmp8_1;
                              EM_S[INDEX2(3,0,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                              EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                              EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                              EM_S[INDEX2(0,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                              EM_S[INDEX2(2,0,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                              EM_S[INDEX2(1,3,4)]+=tmp11_1 + tmp12_1 + tmp1_1 + tmp2_1;
                              EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                              EM_S[INDEX2(2,2,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                              EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp3_1;
                              EM_S[INDEX2(0,3,4)]+=tmp4_1 + tmp6_1 + tmp9_1;
                              EM_S[INDEX2(1,1,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                           }
                        }
                        /**************************************************************/
                        /*   process B: */
                        /**************************************************************/
                        if (NULL!=B_p) {
                           add_EM_S=TRUE;
                           if (extendedB) {
                              const register double B_0_0 = B_p[INDEX2(0,0,2)];
                              const register double B_1_0 = B_p[INDEX2(1,0,2)];
                              const register double B_0_1 = B_p[INDEX2(0,1,2)];
                              const register double B_1_1 = B_p[INDEX2(1,1,2)];
                              const register double B_0_2 = B_p[INDEX2(0,2,2)];
                              const register double B_1_2 = B_p[INDEX2(1,2,2)];
                              const register double B_0_3 = B_p[INDEX2(0,3,2)];
                              const register double B_1_3 = B_p[INDEX2(1,3,2)];
                              const register double tmp3_0 = B_0_0 + B_0_2;
                              const register double tmp1_0 = B_1_2 + B_1_3;
                              const register double tmp2_0 = B_0_1 + B_0_3;
                              const register double tmp0_0 = B_1_0 + B_1_1;
                              const register double tmp63_1 = B_1_1*w42;
                              const register double tmp79_1 = B_1_1*w40;
                              const register double tmp37_1 = tmp3_0*w35;
                              const register double tmp8_1 = tmp0_0*w32;
                              const register double tmp71_1 = B_0_1*w34;
                              const register double tmp19_1 = B_0_3*w31;
                              const register double tmp15_1 = B_0_3*w34;
                              const register double tmp9_1 = tmp3_0*w34;
                              const register double tmp35_1 = B_1_0*w36;
                              const register double tmp66_1 = B_0_3*w28;
                              const register double tmp28_1 = B_1_0*w42;
                              const register double tmp22_1 = B_1_0*w40;
                              const register double tmp16_1 = B_1_2*w29;
                              const register double tmp6_1 = tmp2_0*w35;
                              const register double tmp55_1 = B_1_3*w40;
                              const register double tmp50_1 = B_1_3*w42;
                              const register double tmp7_1 = tmp1_0*w29;
                              const register double tmp1_1 = tmp1_0*w32;
                              const register double tmp57_1 = B_0_3*w30;
                              const register double tmp18_1 = B_1_1*w32;
                              const register double tmp53_1 = B_1_0*w41;
                              const register double tmp61_1 = B_1_3*w36;
                              const register double tmp27_1 = B_0_3*w38;
                              const register double tmp64_1 = B_0_2*w30;
                              const register double tmp76_1 = B_0_1*w38;
                              const register double tmp39_1 = tmp2_0*w34;
                              const register double tmp62_1 = B_0_1*w31;
                              const register double tmp56_1 = B_0_0*w31;
                              const register double tmp49_1 = B_1_1*w36;
                              const register double tmp2_1 = B_0_2*w31;
                              const register double tmp23_1 = B_0_2*w33;
                              const register double tmp38_1 = B_1_1*w43;
                              const register double tmp74_1 = B_1_2*w41;
                              const register double tmp43_1 = B_1_1*w41;
                              const register double tmp58_1 = B_0_2*w28;
                              const register double tmp67_1 = B_0_0*w33;
                              const register double tmp33_1 = tmp0_0*w39;
                              const register double tmp4_1 = B_0_0*w28;
                              const register double tmp20_1 = B_0_0*w30;
                              const register double tmp13_1 = B_0_2*w38;
                              const register double tmp65_1 = B_1_2*w43;
                              const register double tmp0_1 = tmp0_0*w29;
                              const register double tmp41_1 = tmp3_0*w33;
                              const register double tmp73_1 = B_0_2*w37;
                              const register double tmp69_1 = B_0_0*w38;
                              const register double tmp48_1 = B_1_2*w39;
                              const register double tmp59_1 = B_0_1*w33;
                              const register double tmp17_1 = B_1_3*w41;
                              const register double tmp5_1 = B_0_3*w33;
                              const register double tmp3_1 = B_0_1*w30;
                              const register double tmp21_1 = B_0_1*w28;
                              const register double tmp42_1 = B_1_0*w29;
                              const register double tmp54_1 = B_1_2*w32;
                              const register double tmp60_1 = B_1_0*w39;
                              const register double tmp32_1 = tmp1_0*w36;
                              const register double tmp10_1 = B_0_1*w37;
                              const register double tmp14_1 = B_0_0*w35;
                              const register double tmp29_1 = B_0_1*w35;
                              const register double tmp26_1 = B_1_2*w36;
                              const register double tmp30_1 = B_1_3*w43;
                              const register double tmp70_1 = B_0_2*w35;
                              const register double tmp34_1 = B_1_3*w39;
                              const register double tmp51_1 = B_1_0*w43;
                              const register double tmp31_1 = B_0_2*w34;
                              const register double tmp45_1 = tmp3_0*w28;
                              const register double tmp11_1 = tmp1_0*w39;
                              const register double tmp52_1 = B_1_1*w29;
                              const register double tmp44_1 = B_1_3*w32;
                              const register double tmp25_1 = B_1_1*w39;
                              const register double tmp47_1 = tmp2_0*w33;
                              const register double tmp72_1 = B_1_3*w29;
                              const register double tmp40_1 = tmp2_0*w28;
                              const register double tmp46_1 = B_1_2*w40;
                              const register double tmp36_1 = B_1_2*w42;
                              const register double tmp24_1 = B_0_0*w37;
                              const register double tmp77_1 = B_0_3*w35;
                              const register double tmp68_1 = B_0_3*w37;
                              const register double tmp78_1 = B_0_0*w34;
                              const register double tmp12_1 = tmp0_0*w36;
                              const register double tmp75_1 = B_1_0*w32;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                              EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                              EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                              EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                              EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                              EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp6_1 + tmp9_1;
                              EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                              EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                              EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                              EM_S[INDEX2(2,0,4)]+=tmp45_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                              EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp39_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                              EM_S[INDEX2(2,3,4)]+=tmp11_1 + tmp12_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                              EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                              EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                              EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp7_1 + tmp8_1;
                              EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                           } else { /* constant data */
                              const register double B_0 = B_p[0];
                              const register double B_1 = B_p[1];
                              const register double tmp6_1 = B_1*w50;
                              const register double tmp1_1 = B_1*w45;
                              const register double tmp5_1 = B_1*w49;
                              const register double tmp4_1 = B_1*w48;
                              const register double tmp0_1 = B_0*w44;
                              const register double tmp2_1 = B_0*w46;
                              const register double tmp7_1 = B_0*w51;
                              const register double tmp3_1 = B_0*w47;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                              EM_S[INDEX2(1,2,4)]+=tmp1_1 + tmp2_1;
                              EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                              EM_S[INDEX2(0,0,4)]+=tmp0_1 + tmp5_1;
                              EM_S[INDEX2(3,3,4)]+=tmp3_1 + tmp6_1;
                              EM_S[INDEX2(3,0,4)]+=tmp2_1 + tmp4_1;
                              EM_S[INDEX2(3,1,4)]+=tmp2_1 + tmp6_1;
                              EM_S[INDEX2(2,1,4)]+=tmp4_1 + tmp7_1;
                              EM_S[INDEX2(0,2,4)]+=tmp5_1 + tmp7_1;
                              EM_S[INDEX2(2,0,4)]+=tmp6_1 + tmp7_1;
                              EM_S[INDEX2(1,3,4)]+=tmp2_1 + tmp5_1;
                              EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp4_1;
                              EM_S[INDEX2(2,2,4)]+=tmp0_1 + tmp6_1;
                              EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp3_1;
                              EM_S[INDEX2(0,3,4)]+=tmp1_1 + tmp7_1;
                              EM_S[INDEX2(1,1,4)]+=tmp3_1 + tmp5_1;
                           }
                        }
                        /**************************************************************/
                        /*   process C: */
                        /**************************************************************/
                        if (NULL!=C_p) {
                           add_EM_S=TRUE;
                           if (extendedC) {
                              const register double C_0_0 = C_p[INDEX2(0,0,2)];
                              const register double C_1_0 = C_p[INDEX2(1,0,2)];
                              const register double C_0_1 = C_p[INDEX2(0,1,2)];
                              const register double C_1_1 = C_p[INDEX2(1,1,2)];
                              const register double C_0_2 = C_p[INDEX2(0,2,2)];
                              const register double C_1_2 = C_p[INDEX2(1,2,2)];
                              const register double C_0_3 = C_p[INDEX2(0,3,2)];
                              const register double C_1_3 = C_p[INDEX2(1,3,2)];
                              const register double tmp2_0 = C_0_1 + C_0_3;
                              const register double tmp1_0 = C_1_2 + C_1_3;
                              const register double tmp3_0 = C_0_0 + C_0_2;
                              const register double tmp0_0 = C_1_0 + C_1_1;
                              const register double tmp64_1 = C_0_2*w30;
                              const register double tmp14_1 = C_0_2*w28;
                              const register double tmp19_1 = C_0_3*w31;
                              const register double tmp22_1 = C_1_0*w40;
                              const register double tmp37_1 = tmp3_0*w35;
                              const register double tmp29_1 = C_0_1*w35;
                              const register double tmp73_1 = C_0_2*w37;
                              const register double tmp74_1 = C_1_2*w41;
                              const register double tmp52_1 = C_1_3*w39;
                              const register double tmp25_1 = C_1_1*w39;
                              const register double tmp62_1 = C_0_1*w31;
                              const register double tmp79_1 = C_1_1*w40;
                              const register double tmp43_1 = C_1_1*w36;
                              const register double tmp27_1 = C_0_3*w38;
                              const register double tmp28_1 = C_1_0*w42;
                              const register double tmp63_1 = C_1_1*w42;
                              const register double tmp59_1 = C_0_3*w34;
                              const register double tmp72_1 = C_1_3*w29;
                              const register double tmp40_1 = tmp2_0*w35;
                              const register double tmp13_1 = C_0_3*w30;
                              const register double tmp51_1 = C_1_2*w40;
                              const register double tmp54_1 = C_1_2*w42;
                              const register double tmp12_1 = C_0_0*w31;
                              const register double tmp2_1 = tmp1_0*w32;
                              const register double tmp68_1 = C_0_2*w31;
                              const register double tmp75_1 = C_1_0*w32;
                              const register double tmp49_1 = C_1_1*w41;
                              const register double tmp4_1 = C_0_2*w35;
                              const register double tmp66_1 = C_0_3*w28;
                              const register double tmp56_1 = C_0_1*w37;
                              const register double tmp5_1 = C_0_1*w34;
                              const register double tmp38_1 = tmp2_0*w34;
                              const register double tmp21_1 = C_0_1*w28;
                              const register double tmp69_1 = C_0_1*w30;
                              const register double tmp53_1 = C_1_0*w36;
                              const register double tmp42_1 = C_1_2*w39;
                              const register double tmp76_1 = C_0_1*w38;
                              const register double tmp32_1 = tmp1_0*w29;
                              const register double tmp45_1 = C_1_0*w43;
                              const register double tmp33_1 = tmp0_0*w32;
                              const register double tmp35_1 = C_1_0*w41;
                              const register double tmp26_1 = C_1_2*w36;
                              const register double tmp67_1 = C_0_0*w33;
                              const register double tmp31_1 = C_0_2*w34;
                              const register double tmp20_1 = C_0_0*w30;
                              const register double tmp70_1 = C_0_0*w28;
                              const register double tmp8_1 = tmp0_0*w39;
                              const register double tmp30_1 = C_1_3*w43;
                              const register double tmp0_1 = tmp0_0*w29;
                              const register double tmp17_1 = C_1_3*w41;
                              const register double tmp58_1 = C_0_0*w35;
                              const register double tmp9_1 = tmp3_0*w33;
                              const register double tmp61_1 = C_1_3*w36;
                              const register double tmp41_1 = tmp3_0*w34;
                              const register double tmp50_1 = C_1_3*w32;
                              const register double tmp18_1 = C_1_1*w32;
                              const register double tmp6_1 = tmp1_0*w36;
                              const register double tmp3_1 = C_0_0*w38;
                              const register double tmp34_1 = C_1_1*w29;
                              const register double tmp77_1 = C_0_3*w35;
                              const register double tmp65_1 = C_1_2*w43;
                              const register double tmp71_1 = C_0_3*w33;
                              const register double tmp55_1 = C_1_1*w43;
                              const register double tmp46_1 = tmp3_0*w28;
                              const register double tmp24_1 = C_0_0*w37;
                              const register double tmp10_1 = tmp1_0*w39;
                              const register double tmp48_1 = C_1_0*w29;
                              const register double tmp15_1 = C_0_1*w33;
                              const register double tmp36_1 = C_1_2*w32;
                              const register double tmp60_1 = C_1_0*w39;
                              const register double tmp47_1 = tmp2_0*w33;
                              const register double tmp16_1 = C_1_2*w29;
                              const register double tmp1_1 = C_0_3*w37;
                              const register double tmp7_1 = tmp2_0*w28;
                              const register double tmp39_1 = C_1_3*w40;
                              const register double tmp44_1 = C_1_3*w42;
                              const register double tmp57_1 = C_0_2*w38;
                              const register double tmp78_1 = C_0_0*w34;
                              const register double tmp11_1 = tmp0_0*w36;
                              const register double tmp23_1 = C_0_2*w33;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                              EM_S[INDEX2(1,2,4)]+=tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                              EM_S[INDEX2(3,2,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                              EM_S[INDEX2(0,0,4)]+=tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1;
                              EM_S[INDEX2(3,3,4)]+=tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1;
                              EM_S[INDEX2(3,0,4)]+=tmp32_1 + tmp33_1 + tmp7_1 + tmp9_1;
                              EM_S[INDEX2(3,1,4)]+=tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                              EM_S[INDEX2(2,1,4)]+=tmp32_1 + tmp33_1 + tmp40_1 + tmp41_1;
                              EM_S[INDEX2(0,2,4)]+=tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1;
                              EM_S[INDEX2(2,0,4)]+=tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1;
                              EM_S[INDEX2(1,3,4)]+=tmp37_1 + tmp38_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1;
                              EM_S[INDEX2(2,3,4)]+=tmp10_1 + tmp11_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1;
                              EM_S[INDEX2(2,2,4)]+=tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1;
                              EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp2_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1;
                              EM_S[INDEX2(0,3,4)]+=tmp40_1 + tmp41_1 + tmp6_1 + tmp8_1;
                              EM_S[INDEX2(1,1,4)]+=tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1;
                           } else { /* constant data */
                              const register double C_0 = C_p[0];
                              const register double C_1 = C_p[1];
                              const register double tmp1_1 = C_1*w45;
                              const register double tmp3_1 = C_0*w51;
                              const register double tmp4_1 = C_0*w44;
                              const register double tmp7_1 = C_0*w46;
                              const register double tmp5_1 = C_1*w49;
                              const register double tmp2_1 = C_1*w48;
                              const register double tmp0_1 = C_0*w47;
                              const register double tmp6_1 = C_1*w50;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                              EM_S[INDEX2(1,2,4)]+=tmp2_1 + tmp3_1;
                              EM_S[INDEX2(3,2,4)]+=tmp2_1 + tmp4_1;
                              EM_S[INDEX2(0,0,4)]+=tmp4_1 + tmp5_1;
                              EM_S[INDEX2(3,3,4)]+=tmp0_1 + tmp6_1;
                              EM_S[INDEX2(3,0,4)]+=tmp1_1 + tmp3_1;
                              EM_S[INDEX2(3,1,4)]+=tmp5_1 + tmp7_1;
                              EM_S[INDEX2(2,1,4)]+=tmp1_1 + tmp7_1;
                              EM_S[INDEX2(0,2,4)]+=tmp3_1 + tmp6_1;
                              EM_S[INDEX2(2,0,4)]+=tmp3_1 + tmp5_1;
                              EM_S[INDEX2(1,3,4)]+=tmp6_1 + tmp7_1;
                              EM_S[INDEX2(2,3,4)]+=tmp0_1 + tmp2_1;
                              EM_S[INDEX2(2,2,4)]+=tmp4_1 + tmp6_1;
                              EM_S[INDEX2(1,0,4)]+=tmp1_1 + tmp4_1;
                              EM_S[INDEX2(0,3,4)]+=tmp2_1 + tmp7_1;
                              EM_S[INDEX2(1,1,4)]+=tmp0_1 + tmp5_1;
                           }
                        }
                        /**************************************************************/
                        /*   process D: */
                        /**************************************************************/
                        if (NULL!=D_p) {
                           add_EM_S=TRUE;
                           if (extendedD) {
                              const register double D_0 = D_p[0];
                              const register double D_1 = D_p[1];
                              const register double D_2 = D_p[2];
                              const register double D_3 = D_p[3];
                              const register double tmp4_0 = D_1 + D_3;
                              const register double tmp2_0 = D_0 + D_1 + D_2 + D_3;
                              const register double tmp5_0 = D_0 + D_2;
                              const register double tmp0_0 = D_0 + D_1;
                              const register double tmp6_0 = D_0 + D_3;
                              const register double tmp1_0 = D_2 + D_3;
                              const register double tmp3_0 = D_1 + D_2;
                              const register double tmp16_1 = D_1*w56;
                              const register double tmp14_1 = tmp6_0*w54;
                              const register double tmp8_1 = D_3*w55;
                              const register double tmp2_1 = tmp2_0*w54;
                              const register double tmp12_1 = tmp5_0*w52;
                              const register double tmp4_1 = tmp0_0*w53;
                              const register double tmp3_1 = tmp1_0*w52;
                              const register double tmp13_1 = tmp4_0*w53;
                              const register double tmp10_1 = tmp4_0*w52;
                              const register double tmp1_1 = tmp1_0*w53;
                              const register double tmp7_1 = D_3*w56;
                              const register double tmp0_1 = tmp0_0*w52;
                              const register double tmp11_1 = tmp5_0*w53;
                              const register double tmp9_1 = D_0*w56;
                              const register double tmp5_1 = tmp3_0*w54;
                              const register double tmp18_1 = D_2*w56;
                              const register double tmp17_1 = D_1*w55;
                              const register double tmp6_1 = D_0*w55;
                              const register double tmp15_1 = D_2*w55;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1 + tmp1_1;
                              EM_S[INDEX2(1,2,4)]+=tmp2_1;
                              EM_S[INDEX2(3,2,4)]+=tmp3_1 + tmp4_1;
                              EM_S[INDEX2(0,0,4)]+=tmp5_1 + tmp6_1 + tmp7_1;
                              EM_S[INDEX2(3,3,4)]+=tmp5_1 + tmp8_1 + tmp9_1;
                              EM_S[INDEX2(3,0,4)]+=tmp2_1;
                              EM_S[INDEX2(3,1,4)]+=tmp10_1 + tmp11_1;
                              EM_S[INDEX2(2,1,4)]+=tmp2_1;
                              EM_S[INDEX2(0,2,4)]+=tmp12_1 + tmp13_1;
                              EM_S[INDEX2(2,0,4)]+=tmp12_1 + tmp13_1;
                              EM_S[INDEX2(1,3,4)]+=tmp10_1 + tmp11_1;
                              EM_S[INDEX2(2,3,4)]+=tmp3_1 + tmp4_1;
                              EM_S[INDEX2(2,2,4)]+=tmp14_1 + tmp15_1 + tmp16_1;
                              EM_S[INDEX2(1,0,4)]+=tmp0_1 + tmp1_1;
                              EM_S[INDEX2(0,3,4)]+=tmp2_1;
                              EM_S[INDEX2(1,1,4)]+=tmp14_1 + tmp17_1 + tmp18_1;
                           } else { /* constant data */
                              const register double D_0 = D_p[0];
                              const register double tmp2_1 = D_0*w59;
                              const register double tmp1_1 = D_0*w58;
                              const register double tmp0_1 = D_0*w57;
                              EM_S[INDEX2(0,1,4)]+=tmp0_1;
                              EM_S[INDEX2(1,2,4)]+=tmp1_1;
                              EM_S[INDEX2(3,2,4)]+=tmp0_1;
                              EM_S[INDEX2(0,0,4)]+=tmp2_1;
                              EM_S[INDEX2(3,3,4)]+=tmp2_1;
                              EM_S[INDEX2(3,0,4)]+=tmp1_1;
                              EM_S[INDEX2(3,1,4)]+=tmp0_1;
                              EM_S[INDEX2(2,1,4)]+=tmp1_1;
                              EM_S[INDEX2(0,2,4)]+=tmp0_1;
                              EM_S[INDEX2(2,0,4)]+=tmp0_1;
                              EM_S[INDEX2(1,3,4)]+=tmp0_1;
                              EM_S[INDEX2(2,3,4)]+=tmp0_1;
                              EM_S[INDEX2(2,2,4)]+=tmp2_1;
                              EM_S[INDEX2(1,0,4)]+=tmp0_1;
                              EM_S[INDEX2(0,3,4)]+=tmp1_1;
                              EM_S[INDEX2(1,1,4)]+=tmp2_1;
                           }
                        }
                        /**************************************************************/
                        /*   process X: */
                        /**************************************************************/
                        if (NULL!=X_p) {
                           add_EM_F=TRUE;
                           if (extendedX) {
                              const register double X_0_0 = X_p[INDEX2(0,0,2)];
                              const register double X_1_0 = X_p[INDEX2(1,0,2)];
                              const register double X_0_1 = X_p[INDEX2(0,1,2)];
                              const register double X_1_1 = X_p[INDEX2(1,1,2)];
                              const register double X_0_2 = X_p[INDEX2(0,2,2)];
                              const register double X_1_2 = X_p[INDEX2(1,2,2)];
                              const register double X_0_3 = X_p[INDEX2(0,3,2)];
                              const register double X_1_3 = X_p[INDEX2(1,3,2)];
                              const register double tmp1_0 = X_1_1 + X_1_3;
                              const register double tmp3_0 = X_0_0 + X_0_1;
                              const register double tmp2_0 = X_1_0 + X_1_2;
                              const register double tmp0_0 = X_0_2 + X_0_3;
                              const register double tmp8_1 = tmp2_0*w66;
                              const register double tmp5_1 = tmp3_0*w64;
                              const register double tmp14_1 = tmp0_0*w64;
                              const register double tmp3_1 = tmp3_0*w60;
                              const register double tmp9_1 = tmp3_0*w63;
                              const register double tmp13_1 = tmp3_0*w65;
                              const register double tmp12_1 = tmp1_0*w66;
                              const register double tmp10_1 = tmp0_0*w60;
                              const register double tmp2_1 = tmp2_0*w61;
                              const register double tmp6_1 = tmp2_0*w62;
                              const register double tmp4_1 = tmp0_0*w65;
                              const register double tmp11_1 = tmp1_0*w67;
                              const register double tmp1_1 = tmp1_0*w62;
                              const register double tmp7_1 = tmp1_0*w61;
                              const register double tmp0_1 = tmp0_0*w63;
                              const register double tmp15_1 = tmp2_0*w67;
                              EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                              EM_F[1]+=tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1;
                              EM_F[2]+=tmp10_1 + tmp11_1 + tmp8_1 + tmp9_1;
                              EM_F[3]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1;
                           } else { /* constant data */
                              const register double X_0 = X_p[0];
                              const register double X_1 = X_p[1];
                              const register double tmp3_1 = X_1*w71;
                              const register double tmp2_1 = X_0*w70;
                              const register double tmp1_1 = X_0*w68;
                              const register double tmp0_1 = X_1*w69;
                              EM_F[0]+=tmp0_1 + tmp1_1;
                              EM_F[1]+=tmp0_1 + tmp2_1;
                              EM_F[2]+=tmp1_1 + tmp3_1;
                              EM_F[3]+=tmp2_1 + tmp3_1;
                           }
                        }
                        /**************************************************************/
                        /*   process Y: */
                        /**************************************************************/
                        if (NULL!=Y_p) {
                           add_EM_F=TRUE;
                           if (extendedY) {
                              const register double Y_0 = Y_p[0];
                              const register double Y_1 = Y_p[1];
                              const register double Y_2 = Y_p[2];
                              const register double Y_3 = Y_p[3];
                              const register double tmp0_0 = Y_1 + Y_2;
                              const register double tmp1_0 = Y_0 + Y_3;
                              const register double tmp9_1 = Y_0*w74;
                              const register double tmp4_1 = tmp1_0*w73;
                              const register double tmp5_1 = Y_2*w74;
                              const register double tmp7_1 = Y_1*w74;
                              const register double tmp6_1 = Y_2*w72;
                              const register double tmp2_1 = Y_3*w74;
                              const register double tmp8_1 = Y_3*w72;
                              const register double tmp3_1 = Y_1*w72;
                              const register double tmp0_1 = Y_0*w72;
                              const register double tmp1_1 = tmp0_0*w73;
                              EM_F[0]+=tmp0_1 + tmp1_1 + tmp2_1;
                              EM_F[1]+=tmp3_1 + tmp4_1 + tmp5_1;
                              EM_F[2]+=tmp4_1 + tmp6_1 + tmp7_1;
                              EM_F[3]+=tmp1_1 + tmp8_1 + tmp9_1;
                           } else { /* constant data */
                              const register double Y_0 = Y_p[0];
                              const register double tmp0_1 = Y_0*w75;
                              EM_F[0]+=tmp0_1;
                              EM_F[1]+=tmp0_1;
                              EM_F[2]+=tmp0_1;
                              EM_F[3]+=tmp0_1;
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
