
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


#ifndef INC_PASO_BLOCKOPS
#define INC_PASO_BLOCKOPS

#include "Common.h"
#include "Paso.h"

void Paso_BlockOps_allMV(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x);


/*  do we need this ? */
#define Paso_BlockOps_Cpy_1(R, V) \
{\
   *(R+0)=*(V+0);\
}\

#define Paso_BlockOps_Cpy_2(R, V) \
{\
   *(R+0)=*(V+0);\
   *(R+1)=*(V+1);\
}\

#define Paso_BlockOps_Cpy_3(R, V) \
{\
   *(R+0)=*(V+0);\
   *(R+1)=*(V+1);\
   *(R+2)=*(V+2);\
}\

#define Paso_BlockOps_Cpy_N(N,R,V) \
{\
   memcpy((void*)R, (void*)V, ( (size_t) N ) * sizeof(double) );\
}\

#define Paso_BlockOps_MV_1( R, MAT, V)   \
{ \
   register double S1=*V;\
   register double A11=*MAT;\
   *R=A11 * S1;\
}\


#define Paso_BlockOps_MV_2(R, MAT, V) \
{ \
   register double S1=*(V+0);\
   register double S2=*(V+1);\
   register double A11=*(MAT+0);\
   register double A12=*(MAT+2);\
   register double A21=*(MAT+1);\
   register double A22=*(MAT+3);\
   *(R+0)  = A11 * S1 + A12 * S2;\
   *(R+1)  = A21 * S1 + A22 * S2;\
}\


#define Paso_BlockOps_MV_3(R, MAT, V) \
{ \
   register double S1=*(V+0);\
   register double S2=*(V+1);\
   register double S3=*(V+2);\
   register double A11=*(MAT+0);\
   register double A21=*(MAT+1);\
   register double A31=*(MAT+2);\
   register double A12=*(MAT+3);\
   register double A22=*(MAT+4);\
   register double A32=*(MAT+5);\
   register double A13=*(MAT+6);\
   register double A23=*(MAT+7);\
   register double A33=*(MAT+8);\
   *(R+0)=A11 * S1 + A12 * S2 + A13 * S3; \
   *(R+1)=A21 * S1 + A22 * S2 + A23 * S3;\
   *(R+2)=A31 * S1 + A32 * S2 + A33 * S3;\
}\

/* this does not work in place */
#define Paso_BlockOps_MV_N(N, R, MAT, V) \
{ \
   dim_t p,q; \
   for (p=0; p<N; p++ ) {    \
      register double S= *(MAT+p) * (*(V+0)); \
      for ( q=1; q<N; q++ ) S+=*(MAT+p+N*q) * (*(V+q));  \
      *(R+p)=S;   \
   }  \
} \

#define Paso_BlockOps_Solve_N(N, R, MAT, pivot, V) \
{ \
   *(pivot+0)=0;\
   Esys_setError(TYPE_ERROR, "Paso_BlockOps_Solve_N: Right now there is support block size less than 4 only"); \
} \

#define Paso_BlockOps_SMV_1(R, MAT, V)   \
{\
   register double S1=*(V+0); \
   register double A11=*(MAT+0);\
   *(R+0)-=A11 * S1;\
}\

#define Paso_BlockOps_SMV_2(R,MAT, V) \
{\
   register double S1=*(V+0); \
   register double S2=*(V+1);\
   register double A11=*(MAT+0);\
   register double A12=*(MAT+2);\
   register double A21=*(MAT+1);\
   register double A22=*(MAT+3);\
   *(R+0)-=A11 * S1 + A12 * S2;\
   *(R+1)-=A21 * S1 + A22 * S2;\
}\

#define Paso_BlockOps_SMV_3(R, MAT, V)  \
{\
   register double S1=*(V+0);\
   register double S2=*(V+1);\
   register double S3=*(V+2);\
   register double A11=*(MAT+0);\
   register double A21=*(MAT+1);\
   register double A31=*(MAT+2);\
   register double A12=*(MAT+3);\
   register double A22=*(MAT+4);\
   register double A32=*(MAT+5);\
   register double A13=*(MAT+6);\
   register double A23=*(MAT+7);\
   register double A33=*(MAT+8);\
   *(R+0)-=A11 * S1 + A12 * S2 + A13 * S3; \
   *(R+1)-=A21 * S1 + A22 * S2 + A23 * S3;\
   *(R+2)-=A31 * S1 + A32 * S2 + A33 * S3;\
}\

/* this does not work in place */
#define Paso_BlockOps_SMV_N(N, R, MAT, V) \
{\
   dim_t p,q; \
   for (p=0; p<N; p++ ) { \
      register double S= *(MAT+p) * *(V+0); \
      for ( q=1; q<N; q++ ) S+=*(MAT+p+N*q) * (*(V+q));  \
      *(R+p)-=S;  \
   }\
} \
	
#endif /* #INC_Paso_BlockOpsOPS */
