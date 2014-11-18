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
/*                                                                                                            */
/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  and -(X_{k,i})_i + Y_k          */
/*                                                                                                            */
/*    u has p.numComp components in a 3D domain. The shape functions for test and solution must be identical  */
/*                                                                                                            */
/*    Shape of the coefficients:                                                                              */
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
void  Finley_Assemble_PDE_System_3D(Finley_Assemble_Parameters p,
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
      const double w56 = -0.048598662630909980773*h1;
      const double w13 = -0.049382716049382716049*h1/h0;
      const double w67 = -0.0034784464623227873914*h0;
      const double w23 = -0.0077160493827160493827*h0/h1;
      const double w60 = -0.0061728395061728395062*h0;
      const double w10 = -0.030864197530864197531*h1/h0;
      const double w74 = 0.0034784464623227873914*h1;
      const double w95 = -0.083333333333333333333*h0;
      const double w127 = 0.10954300427416564056*h1;
      const double w65 = -0.00086961161558069684786*h0;
      const double w98 = 0.083333333333333333333*h0;
      const double w68 = -0.00086961161558069684786*h1;
      const double w128 = 0.061728395061728395062*h1;
      const double w130 = 0.0086961161558069684786*h1;
      const double w75 = 0.00086961161558069684786*h0;
      const double w91 = 0.048598662630909980773*h0;
      const double w121 = -0.098765432098765432099*h1;
      const double w102 = 0.0060748328288637475966*h0*h1;
      const double w27 = -0.049382716049382716049*h0/h1;
      const double w29 = 0.060748328288637475966*h0/h1;
      const double w28 = 0.060748328288637475966*h1/h0;
      const double w148 = 0.25*h0*h1;
      const double w107 = 0.00039202670923636763834*h0*h1;
      const double w86 = -0.048598662630909980773*h0;
      const double w2 = -0.060748328288637475966;
      const double w45 = 0.25000000000000000000;
      const double w55 = -0.0068464377671353525349*h0;
      const double w12 = 0.012345679012345679012*h0/h1;
      const double w25 = -0.030864197530864197531*h0/h1;
      const double w8 = 0.060748328288637475966;
      const double w62 = -0.024691358024691358025*h0;
      const double w125 = -0.013913785849291149566*h1;
      const double w104 = 0.0030864197530864197531*h0*h1;
      const double w109 = 0.047827057692638375834*h0*h1;
      const double w48 = -0.16666666666666666667*h1/h0;
      const double w40 = 0.012345679012345679012*h1/h0;
      const double w119 = -0.061728395061728395062*h1;
      const double w32 = 0.030864197530864197531*h1/h0;
      const double w81 = 0.053901890521502123431*h1;
      const double w19 = 0.00098006677309091909584;
      const double w39 = -0.00098006677309091909584*h0/h1;
      const double w77 = 0.0034784464623227873914*h0;
      const double w14 = 0.049382716049382716049;
      const double w138 = -0.5*h1;
      const double w110 = 0.0000124484802011303384*h0*h1;
      const double w1 = 0.0077160493827160493827;
      const double w42 = -0.0015681068369454705533*h0/h1;
      const double w99 = -0.16666666666666666667*h0;
      const double w115 = -0.068464377671353525349*h0;
      const double w58 = -0.053901890521502123431*h1;
      const double w118 = -0.0086961161558069684786*h0;
      const double w103 = 0.024299331315454990386*h0*h1;
      const double w140 = 0.5*h1;
      const double w22 = -0.0077160493827160493827*h1/h0;
      const double w97 = 0.16666666666666666667*h1;
      const double w146 = 0.0069568929246455747829*h0*h1;
      const double w87 = -0.00078405341847273527667*h0;
      const double w51 = 0.33333333333333333333*h0/h1;
      const double w61 = -0.024691358024691358025*h1;
      const double w83 = 0.048598662630909980773*h1;
      const double w124 = -0.0086961161558069684786*h1;
      const double w117 = -0.061728395061728395062*h0;
      const double w20 = -0.0015681068369454705533*h1/h0;
      const double w84 = 0.027385751068541410139*h0;
      const double w96 = 0.083333333333333333333*h1;
      const double w137 = 0.013913785849291149566*h0;
      const double w11 = 0.0069568929246455747828;
      const double w135 = 0.10954300427416564056*h0;
      const double w26 = -0.012345679012345679012*h0/h1;
      const double w116 = -0.10954300427416564056*h1;
      const double w92 = 0.053901890521502123431*h0;
      const double w38 = 0.0077160493827160493827*h1/h0;
      const double w17 = -0.0069568929246455747828;
      const double w69 = 0.0068464377671353525349*h1;
      const double w113 = 0.11111111111111111111*h0*h1;
      const double w141 = 0.5*h0;
      const double w36 = 0.00098006677309091909584*h1/h0;
      const double w131 = 0.013913785849291149566*h1;
      const double w50 = 0.33333333333333333333*h1/h0;
      const double w16 = 0.049382716049382716049*h0/h1;
      const double w143 = 0.054771502137082820279*h0*h1;
      const double w0 = -0.060748328288637475966*h1/h0;
      const double w144 = 0.0077160493827160493827*h0*h1;
      const double w85 = -0.053901890521502123431*h0;
      const double w79 = 0.0061728395061728395062*h0;
      const double w3 = 0.0077160493827160493827*h0/h1;
      const double w35 = 0.0015681068369454705533*h0/h1;
      const double w53 = -0.33333333333333333333*h0/h1;
      const double w114 = -0.068464377671353525349*h1;
      const double w49 = -0.16666666666666666667*h0/h1;
      const double w126 = 0.068464377671353525349*h1;
      const double w44 = -0.33333333333333333333*h1/h0;
      const double w4 = -0.097197325261819961546*h1/h0;
      const double w82 = 0.0068464377671353525349*h0;
      const double w46 = -0.25000000000000000000;
      const double w139 = -0.5*h0;
      const double w94 = -0.16666666666666666667*h1;
      const double w37 = 0.0015681068369454705533*h1/h0;
      const double w64 = -0.00011045515751022224798*h1;
      const double w78 = 0.00011045515751022224798*h1;
      const double w89 = 0.00011045515751022224798*h0;
      const double w111 = 0.055555555555555555556*h0*h1;
      const double w101 = -0.083333333333333333333*h1;
      const double w136 = 0.098765432098765432099*h0;
      const double w9 = -0.0077160493827160493827;
      const double w31 = 0.00098006677309091909584*h0/h1;
      const double w18 = -0.00098006677309091909584*h1/h0;
      const double w73 = 0.024691358024691358025*h1;
      const double w132 = 0.068464377671353525349*h0;
      const double w122 = -0.098765432098765432099*h0;
      const double w34 = 0.049382716049382716049*h1/h0;
      const double w76 = 0.00078405341847273527667*h1;
      const double w108 = 0.00077160493827160493827*h0*h1;
      const double w54 = -0.0068464377671353525349*h1;
      const double w33 = 0.097197325261819961546*h0/h1;
      const double w47 = 0.16666666666666666667*h0/h1;
      const double w80 = 0.024691358024691358025*h0;
      const double w21 = -0.00098006677309091909584;
      const double w88 = -0.00011045515751022224798*h0;
      const double w129 = 0.098765432098765432099*h1;
      const double w63 = -0.027385751068541410139*h1;
      const double w123 = -0.013913785849291149566*h0;
      const double w30 = 0.097197325261819961546*h1/h0;
      const double w90 = 0.00078405341847273527667*h0;
      const double w105 = 0.012345679012345679012*h0*h1;
      const double w7 = 0.030864197530864197531*h0/h1;
      const double w147 = 0.00098006677309091909584*h0*h1;
      const double w142 = 0.060748328288637475966*h0*h1;
      const double w120 = -0.10954300427416564056*h0;
      const double w66 = -0.00078405341847273527667*h1;
      const double w145 = 0.049382716049382716049*h0*h1;
      const double w71 = 0.00086961161558069684786*h1;
      const double w24 = -0.012345679012345679012*h1/h0;
      const double w43 = -0.097197325261819961546*h0/h1;
      const double w72 = 0.027385751068541410139*h1;
      const double w106 = 0.000098006677309091909584*h0*h1;
      const double w6 = -0.054771502137082820279;
      const double w41 = -0.060748328288637475966*h0/h1;
      const double w59 = -0.0034784464623227873914*h1;
      const double w134 = 0.0086961161558069684786*h0;
      const double w112 = 0.027777777777777777778*h0*h1;
      const double w57 = -0.027385751068541410139*h0;
      const double w52 = 0.16666666666666666667*h1/h0;
      const double w5 = 0.054771502137082820279;
      const double w70 = 0.0061728395061728395062*h1;
      const double w93 = -0.0061728395061728395062*h1;
      const double w100 = 0.16666666666666666667*h0;
      const double w133 = 0.061728395061728395062*h0;
      const double w15 = -0.049382716049382716049;
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
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double A_00_0 = A_p[INDEX5(k,0,m,0,0, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_0 = A_p[INDEX5(k,0,m,1,0, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_0 = A_p[INDEX5(k,1,m,0,0, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_0 = A_p[INDEX5(k,1,m,1,0, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_1 = A_p[INDEX5(k,0,m,0,1, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_1 = A_p[INDEX5(k,0,m,1,1, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_1 = A_p[INDEX5(k,1,m,0,1, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_1 = A_p[INDEX5(k,1,m,1,1, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_2 = A_p[INDEX5(k,0,m,0,2, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_2 = A_p[INDEX5(k,0,m,1,2, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_2 = A_p[INDEX5(k,1,m,0,2, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_2 = A_p[INDEX5(k,1,m,1,2, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_3 = A_p[INDEX5(k,0,m,0,3, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_3 = A_p[INDEX5(k,0,m,1,3, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_3 = A_p[INDEX5(k,1,m,0,3, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_3 = A_p[INDEX5(k,1,m,1,3, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_4 = A_p[INDEX5(k,0,m,0,4, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_4 = A_p[INDEX5(k,0,m,1,4, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_4 = A_p[INDEX5(k,1,m,0,4, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_4 = A_p[INDEX5(k,1,m,1,4, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_5 = A_p[INDEX5(k,0,m,0,5, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_5 = A_p[INDEX5(k,0,m,1,5, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_5 = A_p[INDEX5(k,1,m,0,5, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_5 = A_p[INDEX5(k,1,m,1,5, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_6 = A_p[INDEX5(k,0,m,0,6, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_6 = A_p[INDEX5(k,0,m,1,6, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_6 = A_p[INDEX5(k,1,m,0,6, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_6 = A_p[INDEX5(k,1,m,1,6, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_7 = A_p[INDEX5(k,0,m,0,7, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_7 = A_p[INDEX5(k,0,m,1,7, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_7 = A_p[INDEX5(k,1,m,0,7, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_7 = A_p[INDEX5(k,1,m,1,7, p.numEqu,2,p.numComp,2)];
                                    const register double A_00_8 = A_p[INDEX5(k,0,m,0,8, p.numEqu,2,p.numComp,2)];
                                    const register double A_01_8 = A_p[INDEX5(k,0,m,1,8, p.numEqu,2,p.numComp,2)];
                                    const register double A_10_8 = A_p[INDEX5(k,1,m,0,8, p.numEqu,2,p.numComp,2)];
                                    const register double A_11_8 = A_p[INDEX5(k,1,m,1,8, p.numEqu,2,p.numComp,2)];
                                    const register double tmp34_0 = A_01_5 + A_01_7;
                                    const register double tmp14_0 = A_01_5 + A_01_7 + A_10_1 + A_10_3;
                                    const register double tmp38_0 = A_01_3 + A_01_7 + A_10_3 + A_10_7;
                                    const register double tmp1_0 = A_00_0 + A_00_2;
                                    const register double tmp29_0 = A_01_1 + A_01_5 + A_10_3 + A_10_7;
                                    const register double tmp8_0 = A_10_2 + A_10_6;
                                    const register double tmp27_0 = A_01_6 + A_10_2;
                                    const register double tmp19_0 = A_01_1 + A_01_3 + A_10_5 + A_10_7;
                                    const register double tmp26_0 = A_01_2 + A_10_6;
                                    const register double tmp16_0 = A_00_0 + A_00_2 + A_00_6 + A_00_8;
                                    const register double tmp12_0 = A_01_0 + A_10_8;
                                    const register double tmp0_0 = A_01_0 + A_01_8;
                                    const register double tmp24_0 = A_01_8 + A_10_8;
                                    const register double tmp23_0 = A_01_5 + A_01_7 + A_10_5 + A_10_7;
                                    const register double tmp21_0 = A_11_0 + A_11_6;
                                    const register double tmp32_0 = A_10_1 + A_10_5;
                                    const register double tmp17_0 = A_01_8 + A_10_0;
                                    const register double tmp35_0 = A_01_2 + A_01_6;
                                    const register double tmp20_0 = A_01_0 + A_10_0;
                                    const register double tmp11_0 = A_00_6 + A_00_8;
                                    const register double tmp7_0 = A_10_1 + A_10_3;
                                    const register double tmp37_0 = A_01_6 + A_10_6;
                                    const register double tmp18_0 = A_00_1 + A_00_7;
                                    const register double tmp10_0 = A_10_5 + A_10_7;
                                    const register double tmp9_0 = A_11_3 + A_11_5;
                                    const register double tmp5_0 = A_11_1 + A_11_7;
                                    const register double tmp30_0 = A_01_0 + A_01_8 + A_10_0 + A_10_8;
                                    const register double tmp40_0 = A_01_2 + A_10_2;
                                    const register double tmp22_0 = A_11_2 + A_11_8;
                                    const register double tmp4_0 = A_01_3 + A_01_7;
                                    const register double tmp28_0 = A_01_3 + A_01_7 + A_10_1 + A_10_5;
                                    const register double tmp13_0 = A_01_2 + A_01_6 + A_10_2 + A_10_6;
                                    const register double tmp33_0 = A_10_3 + A_10_7;
                                    const register double tmp6_0 = A_00_3 + A_00_5;
                                    const register double tmp15_0 = A_01_4 + A_10_4;
                                    const register double tmp3_0 = A_01_1 + A_01_5;
                                    const register double tmp25_0 = A_01_1 + A_01_3 + A_10_1 + A_10_3;
                                    const register double tmp31_0 = A_10_0 + A_10_8;
                                    const register double tmp2_0 = A_11_0 + A_11_2 + A_11_6 + A_11_8;
                                    const register double tmp36_0 = A_01_1 + A_01_3;
                                    const register double tmp39_0 = A_01_1 + A_01_5 + A_10_1 + A_10_5;
                                    const register double tmp29_1 = tmp17_0*w19;
                                    const register double tmp100_1 = tmp22_0*w39;
                                    const register double tmp115_1 = tmp40_0*w21;
                                    const register double tmp89_1 = A_11_3*w42;
                                    const register double tmp83_1 = tmp33_0*w11;
                                    const register double tmp75_1 = tmp28_0*w6;
                                    const register double tmp30_1 = A_11_4*w27;
                                    const register double tmp4_1 = tmp3_0*w5;
                                    const register double tmp87_1 = A_01_4*w15;
                                    const register double tmp93_1 = tmp22_0*w41;
                                    const register double tmp81_1 = A_01_8*w2;
                                    const register double tmp99_1 = tmp14_0*w5;
                                    const register double tmp1_1 = tmp1_0*w0;
                                    const register double tmp70_1 = tmp21_0*w31;
                                    const register double tmp95_1 = A_10_6*w19;
                                    const register double tmp38_1 = tmp3_0*w11;
                                    const register double tmp79_1 = tmp21_0*w39;
                                    const register double tmp77_1 = tmp30_0*w9;
                                    const register double tmp25_1 = tmp14_0*w11;
                                    const register double tmp117_1 = tmp26_0*w2;
                                    const register double tmp2_1 = tmp2_0*w3;
                                    const register double tmp72_1 = tmp15_0*w15;
                                    const register double tmp18_1 = A_00_7*w20;
                                    const register double tmp44_1 = A_01_2*w19;
                                    const register double tmp52_1 = tmp24_0*w19;
                                    const register double tmp15_1 = A_01_4*w14;
                                    const register double tmp50_1 = tmp23_0*w11;
                                    const register double tmp102_1 = tmp33_0*w5;
                                    const register double tmp21_1 = tmp11_0*w18;
                                    const register double tmp64_1 = tmp25_0*w11;
                                    const register double tmp109_1 = tmp21_0*w41;
                                    const register double tmp106_1 = A_11_5*w42;
                                    const register double tmp35_1 = A_10_8*w2;
                                    const register double tmp16_1 = tmp10_0*w17;
                                    const register double tmp49_1 = tmp22_0*w31;
                                    const register double tmp118_1 = tmp29_0*w6;
                                    const register double tmp114_1 = tmp39_0*w17;
                                    const register double tmp19_1 = A_10_8*w21;
                                    const register double tmp68_1 = A_11_5*w33;
                                    const register double tmp73_1 = tmp26_0*w21;
                                    const register double tmp46_1 = tmp20_0*w8;
                                    const register double tmp36_1 = tmp4_0*w5;
                                    const register double tmp12_1 = tmp9_0*w12;
                                    const register double tmp111_1 = A_10_2*w19;
                                    const register double tmp31_1 = tmp18_0*w24;
                                    const register double tmp105_1 = A_10_6*w8;
                                    const register double tmp10_1 = A_00_4*w13;
                                    const register double tmp82_1 = tmp32_0*w5;
                                    const register double tmp51_1 = A_11_5*w35;
                                    const register double tmp57_1 = tmp1_0*w28;
                                    const register double tmp80_1 = tmp16_0*w38;
                                    const register double tmp48_1 = tmp21_0*w29;
                                    const register double tmp47_1 = tmp6_0*w32;
                                    const register double tmp110_1 = A_01_8*w21;
                                    const register double tmp85_1 = tmp35_0*w9;
                                    const register double tmp28_1 = tmp16_0*w22;
                                    const register double tmp34_1 = tmp11_0*w0;
                                    const register double tmp55_1 = A_00_4*w34;
                                    const register double tmp123_1 = tmp40_0*w2;
                                    const register double tmp91_1 = A_11_5*w43;
                                    const register double tmp7_1 = tmp5_0*w7;
                                    const register double tmp41_1 = tmp7_0*w17;
                                    const register double tmp0_1 = tmp0_0*w1;
                                    const register double tmp32_1 = tmp2_0*w23;
                                    const register double tmp98_1 = tmp12_0*w19;
                                    const register double tmp107_1 = tmp34_0*w17;
                                    const register double tmp76_1 = tmp29_0*w17;
                                    const register double tmp120_1 = tmp37_0*w21;
                                    const register double tmp40_1 = A_01_6*w8;
                                    const register double tmp53_1 = A_00_1*w30;
                                    const register double tmp116_1 = tmp27_0*w21;
                                    const register double tmp101_1 = A_01_0*w2;
                                    const register double tmp45_1 = tmp1_0*w18;
                                    const register double tmp67_1 = A_00_7*w30;
                                    const register double tmp61_1 = A_00_1*w37;
                                    const register double tmp17_1 = A_11_4*w16;
                                    const register double tmp66_1 = tmp20_0*w19;
                                    const register double tmp84_1 = tmp34_0*w6;
                                    const register double tmp92_1 = tmp18_0*w40;
                                    const register double tmp23_1 = tmp13_0*w1;
                                    const register double tmp97_1 = tmp19_0*w11;
                                    const register double tmp5_1 = A_00_1*w4;
                                    const register double tmp62_1 = tmp24_0*w8;
                                    const register double tmp39_1 = tmp10_0*w6;
                                    const register double tmp9_1 = tmp7_0*w6;
                                    const register double tmp63_1 = tmp1_0*w36;
                                    const register double tmp65_1 = A_11_3*w35;
                                    const register double tmp69_1 = tmp11_0*w28;
                                    const register double tmp96_1 = tmp17_0*w8;
                                    const register double tmp108_1 = A_11_3*w43;
                                    const register double tmp14_1 = A_10_4*w15;
                                    const register double tmp6_1 = tmp4_0*w11;
                                    const register double tmp3_1 = A_10_0*w2;
                                    const register double tmp22_1 = tmp12_0*w8;
                                    const register double tmp26_1 = tmp9_0*w26;
                                    const register double tmp8_1 = tmp6_0*w10;
                                    const register double tmp113_1 = tmp38_0*w6;
                                    const register double tmp33_1 = tmp19_0*w5;
                                    const register double tmp60_1 = tmp22_0*w29;
                                    const register double tmp27_1 = tmp15_0*w14;
                                    const register double tmp54_1 = A_11_3*w33;
                                    const register double tmp71_1 = tmp23_0*w5;
                                    const register double tmp37_1 = A_00_7*w4;
                                    const register double tmp43_1 = A_10_0*w21;
                                    const register double tmp59_1 = A_00_7*w37;
                                    const register double tmp58_1 = tmp25_0*w5;
                                    const register double tmp78_1 = tmp31_0*w1;
                                    const register double tmp20_1 = A_01_6*w19;
                                    const register double tmp119_1 = tmp28_0*w17;
                                    const register double tmp121_1 = tmp39_0*w6;
                                    const register double tmp122_1 = tmp38_0*w17;
                                    const register double tmp88_1 = A_10_4*w14;
                                    const register double tmp56_1 = tmp11_0*w36;
                                    const register double tmp103_1 = tmp32_0*w11;
                                    const register double tmp104_1 = tmp36_0*w6;
                                    const register double tmp94_1 = A_01_0*w21;
                                    const register double tmp11_1 = tmp8_0*w9;
                                    const register double tmp86_1 = A_10_2*w8;
                                    const register double tmp42_1 = A_00_1*w20;
                                    const register double tmp74_1 = tmp27_0*w2;
                                    const register double tmp24_1 = tmp5_0*w25;
                                    const register double tmp112_1 = tmp37_0*w2;
                                    const register double tmp90_1 = tmp36_0*w17;
                                    const register double tmp13_1 = A_01_2*w8;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp21_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp14_1 + tmp15_1 + tmp17_1 + tmp2_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp17_1 + tmp23_1 + tmp27_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp17_1 + tmp23_1 + tmp27_1 + tmp47_1 + tmp55_1 + tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp24_1 + tmp26_1 + tmp28_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp24_1 + tmp30_1 + tmp47_1 + tmp55_1 + tmp78_1 + tmp79_1 + tmp80_1 + tmp81_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp23_1 + tmp24_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp8_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp100_1 + tmp101_1 + tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp24_1 + tmp30_1 + tmp47_1 + tmp55_1 + tmp78_1 + tmp80_1 + tmp85_1 + tmp87_1 + tmp88_1 + tmp92_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp100_1 + tmp106_1 + tmp108_1 + tmp109_1 + tmp11_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp19_1 + tmp24_1 + tmp30_1 + tmp36_1 + tmp38_1 + tmp3_1 + tmp40_1 + tmp44_1 + tmp47_1 + tmp55_1 + tmp80_1 + tmp92_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp11_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp20_1 + tmp24_1 + tmp30_1 + tmp35_1 + tmp39_1 + tmp41_1 + tmp43_1 + tmp47_1 + tmp4_1 + tmp55_1 + tmp6_1 + tmp79_1 + tmp80_1 + tmp89_1 + tmp91_1 + tmp92_1 + tmp93_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp102_1 + tmp103_1 + tmp105_1 + tmp10_1 + tmp111_1 + tmp12_1 + tmp17_1 + tmp2_1 + tmp34_1 + tmp37_1 + tmp42_1 + tmp45_1 + tmp78_1 + tmp7_1 + tmp81_1 + tmp84_1 + tmp85_1 + tmp87_1 + tmp88_1 + tmp8_1 + tmp90_1 + tmp94_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp17_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp51_1 + tmp54_1 + tmp55_1 + tmp61_1 + tmp63_1 + tmp67_1 + tmp69_1 + tmp72_1 + tmp77_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp101_1 + tmp104_1 + tmp107_1 + tmp10_1 + tmp110_1 + tmp12_1 + tmp17_1 + tmp18_1 + tmp1_1 + tmp21_1 + tmp2_1 + tmp5_1 + tmp78_1 + tmp7_1 + tmp82_1 + tmp83_1 + tmp85_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp8_1 + tmp95_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp24_1 + tmp26_1 + tmp28_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp72_1 + tmp77_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp120_1 + tmp121_1 + tmp122_1 + tmp123_1 + tmp17_1 + tmp47_1 + tmp53_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp59_1 + tmp60_1 + tmp65_1 + tmp68_1 + tmp70_1 + tmp72_1 + tmp77_1 + tmp7_1;
                                 }
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double A_00 = A_p[INDEX4(k,0,m,0 p.numEqu,2, p.numComp)];
                                    const register double A_01 = A_p[INDEX4(k,0,m,1 p.numEqu,2, p.numComp)];
                                    const register double A_10 = A_p[INDEX4(k,1,m,0 p.numEqu,2, p.numComp)];
                                    const register double A_11 = A_p[INDEX4(k,1,m,1 p.numEqu,2, p.numComp)];
                                    const register double tmp0_0 = A_01 + A_10;
                                    const register double tmp5_1 = A_00*w48;
                                    const register double tmp3_1 = A_01*w45;
                                    const register double tmp4_1 = A_11*w49;
                                    const register double tmp11_1 = A_00*w52;
                                    const register double tmp7_1 = A_00*w50;
                                    const register double tmp2_1 = A_10*w46;
                                    const register double tmp0_1 = A_00*w44;
                                    const register double tmp13_1 = A_10*w45;
                                    const register double tmp10_1 = A_01*w46;
                                    const register double tmp6_1 = tmp0_0*w45;
                                    const register double tmp8_1 = A_11*w51;
                                    const register double tmp12_1 = A_11*w53;
                                    const register double tmp9_1 = tmp0_0*w46;
                                    const register double tmp1_1 = A_11*w47;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp6_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp6_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp11_1 + tmp12_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp11_1 + tmp12_1 + tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp13_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp7_1 + tmp8_1 + tmp9_1;
                                 }
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process B: */
                        /**************************************************************/
                        if (NULL!=B_p) {
                           add_EM_S=TRUE;
                           if (extendedB) {
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double B_0_0 = B_p[INDEX4(k,0,m,0, p.numEqu,2,p.numComp)];
                                    const register double B_1_0 = B_p[INDEX4(k,1,m,0, p.numEqu,2,p.numComp)];
                                    const register double B_0_1 = B_p[INDEX4(k,0,m,1, p.numEqu,2,p.numComp)];
                                    const register double B_1_1 = B_p[INDEX4(k,1,m,1, p.numEqu,2,p.numComp)];
                                    const register double B_0_2 = B_p[INDEX4(k,0,m,2, p.numEqu,2,p.numComp)];
                                    const register double B_1_2 = B_p[INDEX4(k,1,m,2, p.numEqu,2,p.numComp)];
                                    const register double B_0_3 = B_p[INDEX4(k,0,m,3, p.numEqu,2,p.numComp)];
                                    const register double B_1_3 = B_p[INDEX4(k,1,m,3, p.numEqu,2,p.numComp)];
                                    const register double B_0_4 = B_p[INDEX4(k,0,m,4, p.numEqu,2,p.numComp)];
                                    const register double B_1_4 = B_p[INDEX4(k,1,m,4, p.numEqu,2,p.numComp)];
                                    const register double B_0_5 = B_p[INDEX4(k,0,m,5, p.numEqu,2,p.numComp)];
                                    const register double B_1_5 = B_p[INDEX4(k,1,m,5, p.numEqu,2,p.numComp)];
                                    const register double B_0_6 = B_p[INDEX4(k,0,m,6, p.numEqu,2,p.numComp)];
                                    const register double B_1_6 = B_p[INDEX4(k,1,m,6, p.numEqu,2,p.numComp)];
                                    const register double B_0_7 = B_p[INDEX4(k,0,m,7, p.numEqu,2,p.numComp)];
                                    const register double B_1_7 = B_p[INDEX4(k,1,m,7, p.numEqu,2,p.numComp)];
                                    const register double B_0_8 = B_p[INDEX4(k,0,m,8, p.numEqu,2,p.numComp)];
                                    const register double B_1_8 = B_p[INDEX4(k,1,m,8, p.numEqu,2,p.numComp)];
                                    const register double tmp1_0 = B_1_6 + B_1_8;
                                    const register double tmp2_0 = B_1_0 + B_1_2;
                                    const register double tmp3_0 = B_0_2 + B_0_8;
                                    const register double tmp5_0 = B_0_0 + B_0_6;
                                    const register double tmp4_0 = B_0_1 + B_0_7;
                                    const register double tmp0_0 = B_1_3 + B_1_5;
                                    const register double tmp37_1 = B_1_0*w85;
                                    const register double tmp89_1 = B_1_8*w85;
                                    const register double tmp97_1 = B_0_6*w54;
                                    const register double tmp83_1 = B_1_5*w90;
                                    const register double tmp71_1 = tmp4_0*w93;
                                    const register double tmp3_1 = B_0_1*w56;
                                    const register double tmp4_1 = B_0_6*w64;
                                    const register double tmp95_1 = B_0_8*w58;
                                    const register double tmp38_1 = B_0_5*w59;
                                    const register double tmp17_1 = B_0_5*w74;
                                    const register double tmp64_1 = B_1_1*w84;
                                    const register double tmp119_1 = B_0_6*w78;
                                    const register double tmp102_1 = B_0_0*w68;
                                    const register double tmp24_1 = B_1_1*w67;
                                    const register double tmp74_1 = B_1_6*w85;
                                    const register double tmp36_1 = B_0_2*w78;
                                    const register double tmp77_1 = B_1_0*w55;
                                    const register double tmp46_1 = B_1_8*w88;
                                    const register double tmp2_1 = B_0_3*w59;
                                    const register double tmp57_1 = B_1_3*w90;
                                    const register double tmp52_1 = B_0_8*w81;
                                    const register double tmp111_1 = B_0_7*w76;
                                    const register double tmp92_1 = B_0_7*w56;
                                    const register double tmp9_1 = tmp2_0*w55;
                                    const register double tmp99_1 = B_1_8*w75;
                                    const register double tmp0_1 = B_1_1*w57;
                                    const register double tmp19_1 = tmp1_0*w55;
                                    const register double tmp26_1 = B_1_4*w80;
                                    const register double tmp33_1 = B_0_7*w83;
                                    const register double tmp45_1 = B_0_2*w54;
                                    const register double tmp44_1 = B_1_5*w87;
                                    const register double tmp85_1 = B_1_8*w89;
                                    const register double tmp79_1 = tmp5_0*w54;
                                    const register double tmp15_1 = tmp3_0*w71;
                                    const register double tmp112_1 = B_0_8*w71;
                                    const register double tmp23_1 = B_0_3*w72;
                                    const register double tmp109_1 = B_0_1*w83;
                                    const register double tmp21_1 = B_0_4*w73;
                                    const register double tmp54_1 = B_1_2*w82;
                                    const register double tmp29_1 = tmp1_0*w82;
                                    const register double tmp43_1 = B_1_6*w55;
                                    const register double tmp47_1 = B_0_3*w63;
                                    const register double tmp94_1 = B_0_2*w68;
                                    const register double tmp116_1 = B_1_0*w65;
                                    const register double tmp7_1 = B_0_4*w61;
                                    const register double tmp56_1 = B_0_5*w72;
                                    const register double tmp12_1 = B_0_0*w54;
                                    const register double tmp30_1 = B_1_1*w77;
                                    const register double tmp10_1 = B_1_4*w62;
                                    const register double tmp13_1 = B_0_5*w63;
                                    const register double tmp107_1 = B_0_0*w81;
                                    const register double tmp66_1 = B_1_2*w92;
                                    const register double tmp65_1 = tmp5_0*w71;
                                    const register double tmp14_1 = B_1_7*w67;
                                    const register double tmp106_1 = B_0_6*w71;
                                    const register double tmp113_1 = B_1_2*w85;
                                    const register double tmp63_1 = tmp2_0*w82;
                                    const register double tmp93_1 = B_0_0*w64;
                                    const register double tmp28_1 = B_1_7*w84;
                                    const register double tmp75_1 = tmp3_0*w68;
                                    const register double tmp80_1 = B_1_0*w92;
                                    const register double tmp53_1 = B_1_8*w92;
                                    const register double tmp49_1 = B_0_2*w71;
                                    const register double tmp76_1 = B_1_8*w65;
                                    const register double tmp40_1 = B_0_0*w58;
                                    const register double tmp31_1 = B_0_1*w76;
                                    const register double tmp86_1 = B_1_5*w86;
                                    const register double tmp50_1 = B_0_6*w69;
                                    const register double tmp22_1 = tmp5_0*w69;
                                    const register double tmp90_1 = B_1_3*w87;
                                    const register double tmp87_1 = B_1_6*w65;
                                    const register double tmp88_1 = B_1_2*w55;
                                    const register double tmp98_1 = B_1_6*w92;
                                    const register double tmp6_1 = B_0_8*w68;
                                    const register double tmp11_1 = B_0_7*w66;
                                    const register double tmp58_1 = B_0_0*w78;
                                    const register double tmp69_1 = tmp3_0*w69;
                                    const register double tmp48_1 = B_0_6*w68;
                                    const register double tmp68_1 = B_1_8*w82;
                                    const register double tmp108_1 = B_0_2*w69;
                                    const register double tmp96_1 = B_0_1*w66;
                                    const register double tmp101_1 = B_1_0*w82;
                                    const register double tmp34_1 = tmp0_0*w79;
                                    const register double tmp70_1 = B_1_6*w89;
                                    const register double tmp18_1 = tmp4_0*w70;
                                    const register double tmp59_1 = B_1_6*w75;
                                    const register double tmp72_1 = tmp5_0*w68;
                                    const register double tmp51_1 = B_0_3*w74;
                                    const register double tmp25_1 = B_0_0*w71;
                                    const register double tmp118_1 = B_1_6*w88;
                                    const register double tmp61_1 = tmp1_0*w75;
                                    const register double tmp78_1 = B_1_2*w88;
                                    const register double tmp16_1 = B_1_7*w57;
                                    const register double tmp35_1 = tmp2_0*w75;
                                    const register double tmp5_1 = B_0_2*w58;
                                    const register double tmp105_1 = B_1_2*w89;
                                    const register double tmp100_1 = B_0_2*w64;
                                    const register double tmp84_1 = B_1_2*w75;
                                    const register double tmp114_1 = B_0_2*w81;
                                    const register double tmp62_1 = B_1_7*w77;
                                    const register double tmp91_1 = B_1_0*w88;
                                    const register double tmp110_1 = B_0_8*w78;
                                    const register double tmp27_1 = B_0_6*w81;
                                    const register double tmp32_1 = B_0_8*w69;
                                    const register double tmp42_1 = B_1_2*w65;
                                    const register double tmp8_1 = tmp1_0*w65;
                                    const register double tmp1_1 = tmp0_0*w60;
                                    const register double tmp104_1 = B_0_8*w54;
                                    const register double tmp117_1 = B_0_0*w69;
                                    const register double tmp67_1 = B_1_0*w75;
                                    const register double tmp20_1 = tmp2_0*w65;
                                    const register double tmp103_1 = B_0_6*w58;
                                    const register double tmp115_1 = B_1_8*w55;
                                    const register double tmp60_1 = B_1_0*w89;
                                    const register double tmp41_1 = B_1_3*w86;
                                    const register double tmp55_1 = B_1_5*w91;
                                    const register double tmp82_1 = B_1_6*w82;
                                    const register double tmp39_1 = B_0_8*w64;
                                    const register double tmp73_1 = tmp3_0*w54;
                                    const register double tmp81_1 = B_1_3*w91;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp17_1 + tmp21_1 + tmp23_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp14_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp3_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp21_1 + tmp26_1 + tmp28_1 + tmp30_1 + tmp31_1 + tmp33_1 + tmp49_1 + tmp50_1 + tmp51_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp60_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp15_1 + tmp17_1 + tmp18_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp26_1 + tmp34_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp18_1 + tmp21_1 + tmp26_1 + tmp51_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp62_1 + tmp64_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp69_1 + tmp70_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp13_1 + tmp26_1 + tmp2_1 + tmp34_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp16_1 + tmp24_1 + tmp38_1 + tmp41_1 + tmp44_1 + tmp47_1 + tmp71_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp78_1 + tmp79_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp26_1 + tmp38_1 + tmp47_1 + tmp62_1 + tmp64_1 + tmp71_1 + tmp75_1 + tmp79_1 + tmp7_1 + tmp80_1 + tmp81_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp16_1 + tmp18_1 + tmp21_1 + tmp24_1 + tmp51_1 + tmp56_1 + tmp65_1 + tmp69_1 + tmp86_1 + tmp87_1 + tmp88_1 + tmp89_1 + tmp90_1 + tmp91_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp13_1 + tmp26_1 + tmp28_1 + tmp29_1 + tmp2_1 + tmp30_1 + tmp34_1 + tmp35_1 + tmp7_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1 + tmp96_1 + tmp97_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp100_1 + tmp101_1 + tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp26_1 + tmp28_1 + tmp30_1 + tmp38_1 + tmp47_1 + tmp7_1 + tmp81_1 + tmp83_1 + tmp92_1 + tmp96_1 + tmp98_1 + tmp99_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp106_1 + tmp107_1 + tmp108_1 + tmp109_1 + tmp10_1 + tmp110_1 + tmp111_1 + tmp14_1 + tmp17_1 + tmp1_1 + tmp21_1 + tmp23_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp13_1 + tmp16_1 + tmp19_1 + tmp1_1 + tmp20_1 + tmp24_1 + tmp2_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp109_1 + tmp10_1 + tmp111_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp14_1 + tmp21_1 + tmp51_1 + tmp56_1 + tmp86_1 + tmp90_1;
                                 }
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double B_0 = B_p[INDEX3(k,0,m, p.numEqu,2)];
                                    const register double B_1 = B_p[INDEX3(k,1,m, p.numEqu,2)];
                                    const register double tmp1_1 = B_1*w95;
                                    const register double tmp2_1 = B_0*w96;
                                    const register double tmp0_1 = B_0*w94;
                                    const register double tmp6_1 = B_1*w100;
                                    const register double tmp7_1 = B_0*w101;
                                    const register double tmp4_1 = B_0*w97;
                                    const register double tmp3_1 = B_1*w98;
                                    const register double tmp5_1 = B_1*w99;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp3_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp3_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp5_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1;
                                 }
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process C: */
                        /**************************************************************/
                        if (NULL!=C_p) {
                           add_EM_S=TRUE;
                           if (extendedC) {
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double C_0_0 = C_p[INDEX4(k,m,0, 0, p.numEqu,p.numComp,2)];
                                    const register double C_1_0 = C_p[INDEX4(k,m,1, 0, p.numEqu,p.numComp,2)];
                                    const register double C_0_1 = C_p[INDEX4(k,m,0, 1, p.numEqu,p.numComp,2)];
                                    const register double C_1_1 = C_p[INDEX4(k,m,1, 1, p.numEqu,p.numComp,2)];
                                    const register double C_0_2 = C_p[INDEX4(k,m,0, 2, p.numEqu,p.numComp,2)];
                                    const register double C_1_2 = C_p[INDEX4(k,m,1, 2, p.numEqu,p.numComp,2)];
                                    const register double C_0_3 = C_p[INDEX4(k,m,0, 3, p.numEqu,p.numComp,2)];
                                    const register double C_1_3 = C_p[INDEX4(k,m,1, 3, p.numEqu,p.numComp,2)];
                                    const register double C_0_4 = C_p[INDEX4(k,m,0, 4, p.numEqu,p.numComp,2)];
                                    const register double C_1_4 = C_p[INDEX4(k,m,1, 4, p.numEqu,p.numComp,2)];
                                    const register double C_0_5 = C_p[INDEX4(k,m,0, 5, p.numEqu,p.numComp,2)];
                                    const register double C_1_5 = C_p[INDEX4(k,m,1, 5, p.numEqu,p.numComp,2)];
                                    const register double C_0_6 = C_p[INDEX4(k,m,0, 6, p.numEqu,p.numComp,2)];
                                    const register double C_1_6 = C_p[INDEX4(k,m,1, 6, p.numEqu,p.numComp,2)];
                                    const register double C_0_7 = C_p[INDEX4(k,m,0, 7, p.numEqu,p.numComp,2)];
                                    const register double C_1_7 = C_p[INDEX4(k,m,1, 7, p.numEqu,p.numComp,2)];
                                    const register double C_0_8 = C_p[INDEX4(k,m,0, 8, p.numEqu,p.numComp,2)];
                                    const register double C_1_8 = C_p[INDEX4(k,m,1, 8, p.numEqu,p.numComp,2)];
                                    const register double tmp5_0 = C_0_2 + C_0_8;
                                    const register double tmp2_0 = C_1_0 + C_1_2;
                                    const register double tmp0_0 = C_1_3 + C_1_5;
                                    const register double tmp3_0 = C_0_1 + C_0_7;
                                    const register double tmp4_0 = C_0_0 + C_0_6;
                                    const register double tmp1_0 = C_1_6 + C_1_8;
                                    const register double tmp68_1 = C_1_1*w67;
                                    const register double tmp110_1 = C_0_8*w68;
                                    const register double tmp85_1 = tmp4_0*w54;
                                    const register double tmp25_1 = C_0_5*w63;
                                    const register double tmp43_1 = C_1_2*w65;
                                    const register double tmp19_1 = C_0_4*w61;
                                    const register double tmp41_1 = C_0_0*w58;
                                    const register double tmp9_1 = C_0_2*w69;
                                    const register double tmp103_1 = C_1_0*w82;
                                    const register double tmp47_1 = C_0_2*w54;
                                    const register double tmp51_1 = C_0_2*w71;
                                    const register double tmp95_1 = C_1_6*w89;
                                    const register double tmp72_1 = C_1_6*w65;
                                    const register double tmp7_1 = tmp2_0*w55;
                                    const register double tmp4_1 = C_0_0*w81;
                                    const register double tmp30_1 = C_0_2*w68;
                                    const register double tmp79_1 = tmp4_0*w69;
                                    const register double tmp82_1 = C_1_3*w91;
                                    const register double tmp63_1 = C_1_6*w75;
                                    const register double tmp71_1 = C_1_5*w86;
                                    const register double tmp67_1 = tmp1_0*w55;
                                    const register double tmp20_1 = tmp4_0*w68;
                                    const register double tmp24_1 = tmp5_0*w54;
                                    const register double tmp118_1 = C_1_6*w88;
                                    const register double tmp97_1 = C_0_6*w81;
                                    const register double tmp96_1 = C_0_0*w71;
                                    const register double tmp70_1 = tmp3_0*w70;
                                    const register double tmp34_1 = C_0_1*w66;
                                    const register double tmp113_1 = C_1_2*w85;
                                    const register double tmp99_1 = C_0_2*w78;
                                    const register double tmp54_1 = C_0_8*w81;
                                    const register double tmp42_1 = C_1_3*w86;
                                    const register double tmp87_1 = C_1_8*w89;
                                    const register double tmp61_1 = C_1_3*w90;
                                    const register double tmp48_1 = C_1_8*w88;
                                    const register double tmp55_1 = C_1_8*w92;
                                    const register double tmp0_1 = C_1_1*w57;
                                    const register double tmp65_1 = C_1_7*w57;
                                    const register double tmp94_1 = C_1_8*w82;
                                    const register double tmp116_1 = C_1_0*w65;
                                    const register double tmp1_1 = C_0_6*w71;
                                    const register double tmp101_1 = C_1_8*w75;
                                    const register double tmp32_1 = tmp1_0*w82;
                                    const register double tmp2_1 = C_0_5*w74;
                                    const register double tmp108_1 = C_0_6*w64;
                                    const register double tmp8_1 = C_1_4*w62;
                                    const register double tmp90_1 = C_1_0*w55;
                                    const register double tmp14_1 = C_0_7*w76;
                                    const register double tmp22_1 = C_1_7*w77;
                                    const register double tmp40_1 = C_0_8*w64;
                                    const register double tmp5_1 = C_0_4*w73;
                                    const register double tmp74_1 = C_1_8*w85;
                                    const register double tmp78_1 = tmp5_0*w71;
                                    const register double tmp76_1 = C_1_3*w87;
                                    const register double tmp57_1 = C_0_7*w83;
                                    const register double tmp21_1 = tmp2_0*w82;
                                    const register double tmp52_1 = C_0_6*w69;
                                    const register double tmp39_1 = C_0_1*w56;
                                    const register double tmp88_1 = C_1_6*w85;
                                    const register double tmp31_1 = C_0_8*w58;
                                    const register double tmp38_1 = C_0_5*w59;
                                    const register double tmp10_1 = C_0_3*w72;
                                    const register double tmp102_1 = C_0_2*w64;
                                    const register double tmp80_1 = C_1_0*w92;
                                    const register double tmp49_1 = C_0_3*w63;
                                    const register double tmp75_1 = tmp5_0*w69;
                                    const register double tmp77_1 = C_1_0*w88;
                                    const register double tmp104_1 = C_0_0*w68;
                                    const register double tmp60_1 = C_0_5*w72;
                                    const register double tmp56_1 = C_1_2*w82;
                                    const register double tmp69_1 = tmp4_0*w71;
                                    const register double tmp91_1 = C_1_2*w88;
                                    const register double tmp93_1 = C_1_0*w75;
                                    const register double tmp107_1 = C_1_2*w89;
                                    const register double tmp106_1 = C_0_8*w54;
                                    const register double tmp13_1 = C_0_8*w78;
                                    const register double tmp64_1 = C_1_0*w89;
                                    const register double tmp33_1 = C_1_1*w77;
                                    const register double tmp115_1 = C_1_8*w55;
                                    const register double tmp11_1 = C_0_1*w83;
                                    const register double tmp86_1 = C_1_2*w75;
                                    const register double tmp83_1 = C_1_6*w82;
                                    const register double tmp109_1 = C_0_2*w58;
                                    const register double tmp98_1 = C_0_8*w69;
                                    const register double tmp23_1 = tmp0_0*w79;
                                    const register double tmp114_1 = C_0_2*w81;
                                    const register double tmp81_1 = tmp5_0*w68;
                                    const register double tmp53_1 = C_0_3*w74;
                                    const register double tmp92_1 = C_1_2*w92;
                                    const register double tmp119_1 = C_0_6*w78;
                                    const register double tmp50_1 = C_0_6*w68;
                                    const register double tmp18_1 = tmp1_0*w75;
                                    const register double tmp15_1 = tmp3_0*w93;
                                    const register double tmp73_1 = C_1_2*w55;
                                    const register double tmp36_1 = tmp2_0*w75;
                                    const register double tmp17_1 = C_0_3*w59;
                                    const register double tmp45_1 = C_1_5*w87;
                                    const register double tmp27_1 = C_0_7*w56;
                                    const register double tmp84_1 = C_1_5*w90;
                                    const register double tmp59_1 = C_1_5*w91;
                                    const register double tmp105_1 = C_0_6*w58;
                                    const register double tmp44_1 = C_1_6*w55;
                                    const register double tmp89_1 = C_1_8*w65;
                                    const register double tmp62_1 = C_0_0*w78;
                                    const register double tmp12_1 = C_1_7*w67;
                                    const register double tmp100_1 = C_1_6*w92;
                                    const register double tmp46_1 = C_0_7*w66;
                                    const register double tmp117_1 = C_0_0*w69;
                                    const register double tmp6_1 = tmp1_0*w65;
                                    const register double tmp3_1 = tmp0_0*w60;
                                    const register double tmp37_1 = C_1_0*w85;
                                    const register double tmp111_1 = C_0_0*w54;
                                    const register double tmp66_1 = tmp2_0*w65;
                                    const register double tmp35_1 = C_0_6*w54;
                                    const register double tmp112_1 = C_0_8*w71;
                                    const register double tmp29_1 = C_0_0*w64;
                                    const register double tmp28_1 = C_1_7*w84;
                                    const register double tmp16_1 = C_1_4*w80;
                                    const register double tmp26_1 = C_1_1*w84;
                                    const register double tmp58_1 = C_0_1*w76;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp16_1 + tmp17_1 + tmp19_1 + tmp23_1 + tmp25_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp12_1 + tmp19_1 + tmp37_1 + tmp38_1 + tmp39_1 + tmp40_1 + tmp41_1 + tmp42_1 + tmp43_1 + tmp44_1 + tmp45_1 + tmp46_1 + tmp47_1 + tmp48_1 + tmp49_1 + tmp50_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp16_1 + tmp28_1 + tmp33_1 + tmp51_1 + tmp52_1 + tmp53_1 + tmp54_1 + tmp55_1 + tmp56_1 + tmp57_1 + tmp58_1 + tmp59_1 + tmp5_1 + tmp60_1 + tmp61_1 + tmp62_1 + tmp63_1 + tmp64_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp15_1 + tmp17_1 + tmp19_1 + tmp20_1 + tmp24_1 + tmp25_1 + tmp3_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp53_1 + tmp5_1 + tmp60_1 + tmp65_1 + tmp68_1 + tmp69_1 + tmp70_1 + tmp71_1 + tmp72_1 + tmp73_1 + tmp74_1 + tmp75_1 + tmp76_1 + tmp77_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp2_1 + tmp3_1 + tmp5_1 + tmp65_1 + tmp66_1 + tmp67_1 + tmp68_1 + tmp70_1 + tmp78_1 + tmp79_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp15_1 + tmp16_1 + tmp19_1 + tmp22_1 + tmp26_1 + tmp38_1 + tmp49_1 + tmp80_1 + tmp81_1 + tmp82_1 + tmp83_1 + tmp84_1 + tmp85_1 + tmp86_1 + tmp87_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp15_1 + tmp19_1 + tmp38_1 + tmp42_1 + tmp45_1 + tmp49_1 + tmp65_1 + tmp68_1 + tmp81_1 + tmp85_1 + tmp88_1 + tmp89_1 + tmp8_1 + tmp90_1 + tmp91_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp16_1 + tmp22_1 + tmp26_1 + tmp53_1 + tmp59_1 + tmp5_1 + tmp60_1 + tmp61_1 + tmp69_1 + tmp70_1 + tmp75_1 + tmp92_1 + tmp93_1 + tmp94_1 + tmp95_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp16_1 + tmp23_1 + tmp28_1 + tmp2_1 + tmp32_1 + tmp33_1 + tmp36_1 + tmp57_1 + tmp58_1 + tmp5_1 + tmp96_1 + tmp97_1 + tmp98_1 + tmp99_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp100_1 + tmp101_1 + tmp102_1 + tmp103_1 + tmp104_1 + tmp105_1 + tmp106_1 + tmp107_1 + tmp16_1 + tmp19_1 + tmp27_1 + tmp28_1 + tmp33_1 + tmp34_1 + tmp38_1 + tmp49_1 + tmp82_1 + tmp84_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp108_1 + tmp109_1 + tmp110_1 + tmp111_1 + tmp12_1 + tmp17_1 + tmp19_1 + tmp25_1 + tmp39_1 + tmp3_1 + tmp46_1 + tmp6_1 + tmp7_1 + tmp8_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp16_1 + tmp18_1 + tmp21_1 + tmp22_1 + tmp23_1 + tmp26_1 + tmp2_1 + tmp5_1 + tmp70_1 + tmp78_1 + tmp79_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp112_1 + tmp113_1 + tmp114_1 + tmp115_1 + tmp116_1 + tmp117_1 + tmp118_1 + tmp119_1 + tmp11_1 + tmp12_1 + tmp14_1 + tmp53_1 + tmp5_1 + tmp60_1 + tmp71_1 + tmp76_1 + tmp8_1;
                                 }
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double C_0 = C_p[INDEX3(k,m,0, p.numEqu,p.numComp)];
                                    const register double C_1 = C_p[INDEX3(k,m,1, p.numEqu,p.numComp)];
                                    const register double tmp1_1 = C_1*w95;
                                    const register double tmp3_1 = C_0*w101;
                                    const register double tmp0_1 = C_0*w97;
                                    const register double tmp6_1 = C_1*w100;
                                    const register double tmp7_1 = C_0*w96;
                                    const register double tmp4_1 = C_0*w94;
                                    const register double tmp2_1 = C_1*w98;
                                    const register double tmp5_1 = C_1*w99;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp3_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp5_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp3_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp3_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp2_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp1_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp2_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp5_1;
                                 }
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process D: */
                        /**************************************************************/
                        if (NULL!=D_p) {
                           add_EM_S=TRUE;
                           if (extendedD) {
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double D_0 = D_p[INDEX3(k,m,0, p.numEqu,p.numComp)];
                                    const register double D_1 = D_p[INDEX3(k,m,1, p.numEqu,p.numComp)];
                                    const register double D_2 = D_p[INDEX3(k,m,2, p.numEqu,p.numComp)];
                                    const register double D_3 = D_p[INDEX3(k,m,3, p.numEqu,p.numComp)];
                                    const register double D_4 = D_p[INDEX3(k,m,4, p.numEqu,p.numComp)];
                                    const register double D_5 = D_p[INDEX3(k,m,5, p.numEqu,p.numComp)];
                                    const register double D_6 = D_p[INDEX3(k,m,6, p.numEqu,p.numComp)];
                                    const register double D_7 = D_p[INDEX3(k,m,7, p.numEqu,p.numComp)];
                                    const register double D_8 = D_p[INDEX3(k,m,8, p.numEqu,p.numComp)];
                                    const register double tmp12_0 = D_1 + D_5;
                                    const register double tmp7_0 = D_1 + D_3;
                                    const register double tmp3_0 = D_1 + D_3 + D_5 + D_7;
                                    const register double tmp11_0 = D_0 + D_8;
                                    const register double tmp2_0 = D_0 + D_2;
                                    const register double tmp9_0 = D_1 + D_7;
                                    const register double tmp0_0 = D_6 + D_8;
                                    const register double tmp8_0 = D_0 + D_6;
                                    const register double tmp10_0 = D_2 + D_8;
                                    const register double tmp1_0 = D_3 + D_5;
                                    const register double tmp6_0 = D_5 + D_7;
                                    const register double tmp5_0 = D_2 + D_6;
                                    const register double tmp4_0 = D_0 + D_2 + D_6 + D_8;
                                    const register double tmp13_0 = D_3 + D_7;
                                    const register double tmp2_1 = D_7*w107;
                                    const register double tmp18_1 = tmp7_0*w107;
                                    const register double tmp5_1 = tmp2_0*w102;
                                    const register double tmp25_1 = tmp10_0*w102;
                                    const register double tmp32_1 = tmp12_0*w107;
                                    const register double tmp30_1 = tmp11_0*w108;
                                    const register double tmp35_1 = D_6*w110;
                                    const register double tmp31_1 = D_2*w110;
                                    const register double tmp22_1 = D_5*w103;
                                    const register double tmp24_1 = tmp9_0*w104;
                                    const register double tmp29_1 = tmp8_0*w102;
                                    const register double tmp38_1 = tmp12_0*w103;
                                    const register double tmp9_1 = D_7*w103;
                                    const register double tmp23_1 = D_3*w107;
                                    const register double tmp1_1 = D_1*w103;
                                    const register double tmp27_1 = D_3*w103;
                                    const register double tmp34_1 = tmp13_0*w103;
                                    const register double tmp15_1 = D_0*w109;
                                    const register double tmp36_1 = tmp13_0*w107;
                                    const register double tmp16_1 = tmp7_0*w103;
                                    const register double tmp19_1 = D_8*w109;
                                    const register double tmp26_1 = tmp10_0*w106;
                                    const register double tmp33_1 = D_6*w109;
                                    const register double tmp12_1 = tmp5_0*w108;
                                    const register double tmp8_1 = tmp2_0*w106;
                                    const register double tmp13_1 = D_8*w110;
                                    const register double tmp6_1 = tmp3_0*w104;
                                    const register double tmp20_1 = D_0*w110;
                                    const register double tmp28_1 = D_5*w107;
                                    const register double tmp0_1 = tmp0_0*w106;
                                    const register double tmp37_1 = D_2*w109;
                                    const register double tmp17_1 = tmp6_0*w103;
                                    const register double tmp14_1 = tmp6_0*w107;
                                    const register double tmp21_1 = tmp8_0*w106;
                                    const register double tmp3_1 = tmp1_0*w104;
                                    const register double tmp7_1 = tmp4_0*w108;
                                    const register double tmp11_1 = tmp0_0*w102;
                                    const register double tmp4_1 = D_4*w105;
                                    const register double tmp10_1 = D_1*w107;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp11_1 + tmp3_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp12_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp24_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp24_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp21_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp10_1 + tmp11_1 + tmp3_1 + tmp4_1 + tmp8_1 + tmp9_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp30_1 + tmp31_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp4_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp4_1 + tmp6_1 + tmp7_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp30_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp4_1;
                                 }
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 for (m=0;m<p.numComp;m++) {
                                    const register double D_0 = D_p[INDEX2(k,m, p.numEqu)];
                                    const register double tmp0_1 = D_0*w111;
                                    const register double tmp2_1 = D_0*w113;
                                    const register double tmp1_1 = D_0*w112;
                                    EM_S[INDEX4(k,m,0,1,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,2,p.numEqu,p.numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,2,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,0,p.numEqu,p.numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,3,p.numEqu,p.numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,3,0,p.numEqu,p.numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,3,1,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,1,p.numEqu,p.numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,0,2,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,0,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,1,3,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,3,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,2,2,p.numEqu,p.numComp,4)]+=tmp2_1;
                                    EM_S[INDEX4(k,m,1,0,p.numEqu,p.numComp,4)]+=tmp0_1;
                                    EM_S[INDEX4(k,m,0,3,p.numEqu,p.numComp,4)]+=tmp1_1;
                                    EM_S[INDEX4(k,m,1,1,p.numEqu,p.numComp,4)]+=tmp2_1;
                                 }
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process X: */
                        /**************************************************************/
                        if (NULL!=X_p) {
                           add_EM_F=TRUE;
                           if (extendedX) {
                              for (k=0;k<p.numEqu;k++)   {
                                 const register double X_0_0 = X_p[INDEX3(k,0, 0, p.numEqu,2)];
                                 const register double X_1_0 = X_p[INDEX3(k,1, 0, p.numEqu,2)];
                                 const register double X_0_1 = X_p[INDEX3(k,0, 1, p.numEqu,2)];
                                 const register double X_1_1 = X_p[INDEX3(k,1, 1, p.numEqu,2)];
                                 const register double X_0_2 = X_p[INDEX3(k,0, 2, p.numEqu,2)];
                                 const register double X_1_2 = X_p[INDEX3(k,1, 2, p.numEqu,2)];
                                 const register double X_0_3 = X_p[INDEX3(k,0, 3, p.numEqu,2)];
                                 const register double X_1_3 = X_p[INDEX3(k,1, 3, p.numEqu,2)];
                                 const register double X_0_4 = X_p[INDEX3(k,0, 4, p.numEqu,2)];
                                 const register double X_1_4 = X_p[INDEX3(k,1, 4, p.numEqu,2)];
                                 const register double X_0_5 = X_p[INDEX3(k,0, 5, p.numEqu,2)];
                                 const register double X_1_5 = X_p[INDEX3(k,1, 5, p.numEqu,2)];
                                 const register double X_0_6 = X_p[INDEX3(k,0, 6, p.numEqu,2)];
                                 const register double X_1_6 = X_p[INDEX3(k,1, 6, p.numEqu,2)];
                                 const register double X_0_7 = X_p[INDEX3(k,0, 7, p.numEqu,2)];
                                 const register double X_1_7 = X_p[INDEX3(k,1, 7, p.numEqu,2)];
                                 const register double X_0_8 = X_p[INDEX3(k,0, 8, p.numEqu,2)];
                                 const register double X_1_8 = X_p[INDEX3(k,1, 8, p.numEqu,2)];
                                 const register double tmp0_0 = X_1_0 + X_1_6;
                                 const register double tmp4_0 = X_0_3 + X_0_5;
                                 const register double tmp2_0 = X_0_0 + X_0_2;
                                 const register double tmp5_0 = X_1_2 + X_1_8;
                                 const register double tmp1_0 = X_0_6 + X_0_8;
                                 const register double tmp3_0 = X_1_1 + X_1_7;
                                 const register double tmp24_1 = X_1_3*w135;
                                 const register double tmp14_1 = tmp2_0*w126;
                                 const register double tmp8_1 = X_0_1*w116;
                                 const register double tmp22_1 = tmp2_0*w124;
                                 const register double tmp34_1 = X_1_5*w135;
                                 const register double tmp15_1 = X_0_4*w129;
                                 const register double tmp20_1 = X_1_5*w120;
                                 const register double tmp36_1 = tmp0_0*w134;
                                 const register double tmp39_1 = tmp5_0*w132;
                                 const register double tmp0_1 = X_0_4*w121;
                                 const register double tmp5_1 = X_1_5*w123;
                                 const register double tmp38_1 = X_1_3*w137;
                                 const register double tmp3_1 = X_1_4*w122;
                                 const register double tmp2_1 = tmp1_0*w124;
                                 const register double tmp11_1 = X_0_7*w125;
                                 const register double tmp19_1 = X_0_7*w131;
                                 const register double tmp1_1 = tmp0_0*w115;
                                 const register double tmp10_1 = tmp5_0*w118;
                                 const register double tmp12_1 = tmp5_0*w115;
                                 const register double tmp16_1 = X_0_1*w127;
                                 const register double tmp23_1 = X_1_4*w136;
                                 const register double tmp31_1 = X_0_1*w125;
                                 const register double tmp17_1 = tmp4_0*w128;
                                 const register double tmp37_1 = tmp1_0*w126;
                                 const register double tmp30_1 = tmp0_0*w132;
                                 const register double tmp26_1 = tmp5_0*w134;
                                 const register double tmp29_1 = X_1_5*w137;
                                 const register double tmp7_1 = tmp4_0*w119;
                                 const register double tmp4_1 = tmp2_0*w114;
                                 const register double tmp18_1 = X_1_3*w123;
                                 const register double tmp33_1 = X_0_1*w131;
                                 const register double tmp9_1 = X_1_3*w120;
                                 const register double tmp32_1 = X_0_7*w127;
                                 const register double tmp13_1 = tmp1_0*w130;
                                 const register double tmp35_1 = tmp2_0*w130;
                                 const register double tmp6_1 = tmp3_0*w117;
                                 const register double tmp28_1 = X_0_7*w116;
                                 const register double tmp21_1 = tmp0_0*w118;
                                 const register double tmp27_1 = tmp3_0*w133;
                                 const register double tmp25_1 = tmp1_0*w114;
                                 EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1 + tmp10_1 + tmp11_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_F[INDEX2(k,1,p.numEqu)]+=tmp12_1 + tmp13_1 + tmp14_1 + tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp19_1 + tmp20_1 + tmp21_1 + tmp3_1 + tmp6_1;
                                 EM_F[INDEX2(k,2,p.numEqu)]+=tmp0_1 + tmp22_1 + tmp23_1 + tmp24_1 + tmp25_1 + tmp26_1 + tmp27_1 + tmp28_1 + tmp29_1 + tmp30_1 + tmp31_1 + tmp7_1;
                                 EM_F[INDEX2(k,3,p.numEqu)]+=tmp15_1 + tmp17_1 + tmp23_1 + tmp27_1 + tmp32_1 + tmp33_1 + tmp34_1 + tmp35_1 + tmp36_1 + tmp37_1 + tmp38_1 + tmp39_1;
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 const register double X_0 = X_p[INDEX2(k,0, p.numEqu)];
                                 const register double X_1 = X_p[INDEX2(k,1, p.numEqu)];
                                 const register double tmp2_1 = X_0*w140;
                                 const register double tmp3_1 = X_1*w141;
                                 const register double tmp1_1 = X_1*w139;
                                 const register double tmp0_1 = X_0*w138;
                                 EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1 + tmp1_1;
                                 EM_F[INDEX2(k,1,p.numEqu)]+=tmp1_1 + tmp2_1;
                                 EM_F[INDEX2(k,2,p.numEqu)]+=tmp0_1 + tmp3_1;
                                 EM_F[INDEX2(k,3,p.numEqu)]+=tmp2_1 + tmp3_1;
                              }
                           }
                        }
                        /**************************************************************/
                        /*   process Y: */
                        /**************************************************************/
                        if (NULL!=Y_p) {
                           add_EM_F=TRUE;
                           if (extendedY) {
                              for (k=0;k<p.numEqu;k++) {
                                 const register double Y_0 = Y_p[INDEX3(k,0, p.numEqu)];
                                 const register double Y_1 = Y_p[INDEX3(k,1, p.numEqu)];
                                 const register double Y_2 = Y_p[INDEX3(k,2, p.numEqu)];
                                 const register double Y_3 = Y_p[INDEX3(k,3, p.numEqu)];
                                 const register double Y_4 = Y_p[INDEX3(k,4, p.numEqu)];
                                 const register double Y_5 = Y_p[INDEX3(k,5, p.numEqu)];
                                 const register double Y_6 = Y_p[INDEX3(k,6, p.numEqu)];
                                 const register double Y_7 = Y_p[INDEX3(k,7, p.numEqu)];
                                 const register double Y_8 = Y_p[INDEX3(k,8, p.numEqu)];
                                 const register double tmp5_0 = Y_0 + Y_8;
                                 const register double tmp0_0 = Y_5 + Y_7;
                                 const register double tmp1_0 = Y_1 + Y_3;
                                 const register double tmp4_0 = Y_1 + Y_5;
                                 const register double tmp3_0 = Y_3 + Y_7;
                                 const register double tmp2_0 = Y_2 + Y_6;
                                 const register double tmp14_1 = tmp3_0*w143;
                                 const register double tmp11_1 = tmp4_0*w146;
                                 const register double tmp6_1 = tmp3_0*w146;
                                 const register double tmp17_1 = Y_0*w147;
                                 const register double tmp16_1 = Y_8*w142;
                                 const register double tmp13_1 = Y_2*w147;
                                 const register double tmp18_1 = tmp0_0*w143;
                                 const register double tmp10_1 = tmp5_0*w144;
                                 const register double tmp3_1 = tmp1_0*w143;
                                 const register double tmp12_1 = Y_6*w142;
                                 const register double tmp0_1 = tmp0_0*w146;
                                 const register double tmp9_1 = tmp4_0*w143;
                                 const register double tmp4_1 = tmp2_0*w144;
                                 const register double tmp2_1 = Y_8*w147;
                                 const register double tmp15_1 = tmp1_0*w146;
                                 const register double tmp8_1 = Y_6*w147;
                                 const register double tmp7_1 = Y_2*w142;
                                 const register double tmp5_1 = Y_4*w145;
                                 const register double tmp1_1 = Y_0*w142;
                                 EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1 + tmp1_1 + tmp2_1 + tmp3_1 + tmp4_1 + tmp5_1;
                                 EM_F[INDEX2(k,1,p.numEqu)]+=tmp10_1 + tmp5_1 + tmp6_1 + tmp7_1 + tmp8_1 + tmp9_1;
                                 EM_F[INDEX2(k,2,p.numEqu)]+=tmp10_1 + tmp11_1 + tmp12_1 + tmp13_1 + tmp14_1 + tmp5_1;
                                 EM_F[INDEX2(k,3,p.numEqu)]+=tmp15_1 + tmp16_1 + tmp17_1 + tmp18_1 + tmp4_1 + tmp5_1;
                              }
                           } else { /* constant data */
                              for (k=0;k<p.numEqu;k++) {
                                 const register double Y_0 = Y_p[k];
                                 const register double tmp0_1 = Y_0*w148;
                                 EM_F[INDEX2(k,0,p.numEqu)]+=tmp0_1;
                                 EM_F[INDEX2(k,1,p.numEqu)]+=tmp0_1;
                                 EM_F[INDEX2(k,2,p.numEqu)]+=tmp0_1;
                                 EM_F[INDEX2(k,3,p.numEqu)]+=tmp0_1;
                              }
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