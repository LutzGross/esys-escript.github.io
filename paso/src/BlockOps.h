
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

INLINE void Paso_BlockOps_Cpy_1(double* R, const double* V) 
{
   R[0]=V[0];
}

INLINE void Paso_BlockOps_Cpy_2(double* R, const double* V) 
{
   R[0]=V[0];
   R[1]=V[1];
}

INLINE void Paso_BlockOps_Cpy_3(double* R, const double* V) 
{
   R[0]=V[0];
   R[1]=V[1];
   R[2]=V[2];
}

INLINE void Paso_BlockOps_Cpy_N(const dim_t N,double* R,const double* V) 
{
   memcpy((void*)R, (void*)V, ( (size_t) N ) * sizeof(double) );
}

INLINE void Paso_BlockOps_MV_1(double* R, const double* MAT, const double* V)  
{ 
   register double S1=V[0]; 
   register double A11=MAT[0];
   R[0]=A11 * S1;
}


INLINE void Paso_BlockOps_MV_2(double* R, const double* MAT, const double* V)  
{ 
   register double S1=V[0];
   register double S2=V[1];
   register double A11=MAT[0];
   register double A12=MAT[2];
   register double A21=MAT[1];
   register double A22=MAT[3];
   R[0]  =A11 * S1 + A12 * S2;
   R[1]=A21 * S1 + A22 * S2;
}


INLINE void Paso_BlockOps_MV_3(double* R, const double* MAT, const double* V)  
{ 
   register double S1=V[0];
   register double S2=V[1];
   register double S3=V[2];
   register double A11=MAT[0];
   register double A21=MAT[1];
   register double A31=MAT[2];
   register double A12=MAT[3];
   register double A22=MAT[4];
   register double A32=MAT[5];
   register double A13=MAT[6];
   register double A23=MAT[7];
   register double A33=MAT[8];
   R[0]=A11 * S1 + A12 * S2 + A13 * S3; 
   R[1]=A21 * S1 + A22 * S2 + A23 * S3;
   R[2]=A31 * S1 + A32 * S2 + A33 * S3;
}

/* this does not work in place */
INLINE void Paso_BlockOps_MV_N(const dim_t N, double* R, const double* MAT, const double* V)  
{ 
   dim_t p,q; 
   for (p=0; p<N; p++ ) {                               
      register double S= MAT[p] * V[0];  
      for ( q=1; q<N; q++ ) S+=MAT[p+N*q] * V[q];           
      R[p]=S;                                           
   }                                                
} 

INLINE void Paso_BlockOps_Solve_N(const dim_t N, double* R, const double* MAT, const index_t* pivot, const double* V) 
{ 
   Paso_setError(TYPE_ERROR, "Paso_BlockOps_Solve_N: Right now there is support block size less than 4 only");  
} 

INLINE void Paso_BlockOps_SMV_1(double* R, const double* MAT, const double* V)  
{
   register double S1=V[0]; 
   register double A11=MAT[0];
   R[0]-=A11 * S1;
}

INLINE void Paso_BlockOps_SMV_2(double* R, const double* MAT, const double* V)  
{
   register double S1=V[0]; 
   register double S2=V[1];
   register double A11=MAT[0];
   register double A12=MAT[2];
   register double A21=MAT[1];
   register double A22=MAT[3];
   R[0]-=A11 * S1 + A12 * S2;
   R[1]-=A21 * S1 + A22 * S2;
}

INLINE void Paso_BlockOps_SMV_3(double* R, const double* MAT, const double* V)  
{
   register double S1=V[0];
   register double S2=V[1];
   register double S3=V[2];
   register double A11=MAT[0];
   register double A21=MAT[1];
   register double A31=MAT[2];
   register double A12=MAT[3];
   register double A22=MAT[4];
   register double A32=MAT[5];
   register double A13=MAT[6];
   register double A23=MAT[7];
   register double A33=MAT[8];
   R[0]-=A11 * S1 + A12 * S2 + A13 * S3; 
   R[1]-=A21 * S1 + A22 * S2 + A23 * S3;
   R[2]-=A31 * S1 + A32 * S2 + A33 * S3;
}

/* this does not work in place */
INLINE void Paso_BlockOps_SMV_N(const dim_t N, double* R, const double* MAT, const double* V)  
{ 
   dim_t p,q; 
   for (p=0; p<N; p++ ) {                               
      register double S= MAT[p] * V[0];  
      for ( q=1; q<N; q++ ) S+=MAT[p+N*q] * V[q];           
      R[p]-=S;                                           
   }                                                
} 
	 
#endif /* #INC_Paso_BlockOpsOPS */