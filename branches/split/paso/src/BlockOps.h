
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#ifndef INC_PASO_BLOCKOPS
#define INC_PASO_BLOCKOPS

#include "Common.h"
#include "Paso.h"

#ifdef USE_LAPACK
   #ifdef MKL_LAPACK
      #include <mkl_lapack.h>
      #include <mkl_cblas.h>
   #else
      extern "C"
      {
      #include <clapack.h>
      #include <cblas.h>
      }
   #endif
#endif

void Paso_BlockOps_solveAll(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x);

#define PASO_MISSING_CLAPACK Esys_setError(TYPE_ERROR, "You need to install a CLAPACK version to run a block size > 3.")  


/* copy functions: */
#define Paso_BlockOps_Cpy_2(R, V) \
{\
   *(R+0)=*(V+0);\
   *(R+1)=*(V+1);\
}

#define Paso_BlockOps_Cpy_3(R, V) \
{\
   *(R+0)=*(V+0);\
   *(R+1)=*(V+1);\
   *(R+2)=*(V+2);\
}

#define Paso_BlockOps_Cpy_N(N,R,V)  memcpy((void*)R, (void*)V, ( (size_t) N ) * sizeof(double) )


/* operations R=R-MAT*V (V and R are not overlapping) */

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
}

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
}


/* Paso_BlockOps_SMV_N */

#ifdef USE_LAPACK
   #define Paso_BlockOps_SMV_N(_N_, _R_, _MAT_, _V_) cblas_dgemv(CblasColMajor,CblasNoTrans,_N_,_N_, -1., _MAT_, _N_, _V_, 1, (1.), _R_, 1)  
#else
   #define Paso_BlockOps_SMV_N(N, R, MAT, V)   PASO_MISSING_CLAPACK 
#endif

#ifdef USE_LAPACK
     #define Paso_BlockOps_MV_N(_N_, _R_, _MAT_, _V_) cblas_dgemv(CblasColMajor,CblasNoTrans,_N_,_N_, 1., _MAT_, _N_, _V_, 1, (0.), _R_, 1)  
   #else
      #define Paso_BlockOps_MV_N(N, R, MAT, V)   PASO_MISSING_CLAPACK 
#endif

#define Paso_BlockOps_invM_2(invA, A, failed) \
{\
   register double A11=*(A+0);\
   register double A12=*(A+2);\
   register double A21=*(A+1);\
   register double A22=*(A+3);\
   register double D = A11*A22-A12*A21; \
   if (ABS(D) > 0 ){\
      D=1./D;\
      *(invA)= A22*D;\
      *(invA+1)=-A21*D;\
      *(invA+2)=-A12*D;\
      *(invA+3)= A11*D;\
   } else {\
      *failed=1;\
   }\
}

#define Paso_BlockOps_invM_3(invA, A, failed) \
{\
   register double A11=*(A+0);\
   register double A21=*(A+1);\
   register double A31=*(A+2);\
   register double A12=*(A+3);\
   register double A22=*(A+4);\
   register double A32=*(A+5);\
   register double A13=*(A+6);\
   register double A23=*(A+7);\
   register double A33=*(A+8);\
   register double D =  A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);\
   if (ABS(D) > 0 ){\
      D=1./D;\
      *(invA  )=(A22*A33-A23*A32)*D;\
      *(invA+1)=(A31*A23-A21*A33)*D;\
      *(invA+2)=(A21*A32-A31*A22)*D;\
      *(invA+3)=(A13*A32-A12*A33)*D;\
      *(invA+4)=(A11*A33-A31*A13)*D;\
      *(invA+5)=(A12*A31-A11*A32)*D;\
      *(invA+6)=(A12*A23-A13*A22)*D;\
      *(invA+7)=(A13*A21-A11*A23)*D;\
      *(invA+8)=(A11*A22-A12*A21)*D;\
   } else {\
      *failed=1;\
   }\
}

#ifdef USE_LAPACK
   #ifdef MKL_LAPACK
	 #define Paso_BlockOps_invM_N(N, MAT, PIVOT,failed) \
	 {\
	    int NN=N; \
	    int res =0; \
	    dgetrf(&NN,&NN,MAT,&NN,PIVOT,&res); \
	    if (res!=0) *failed=1;\
	 }
   #else
	 #define Paso_BlockOps_invM_N(N, MAT, PIVOT,failed) \
	 {\
	     int res=clapack_dgetrf(CblasColMajor, N,N,MAT,N, PIVOT); \
	     if (res!=0) *failed=1;\
	 }
   #endif
#else
   #define Paso_BlockOps_invM_N(N, MAT, PIVOT,failed) PASO_MISSING_CLAPACK
#endif

#ifdef USE_LAPACK
   #ifdef MKL_LAPACK
      #define Paso_BlockOps_solve_N(N, X, MAT, PIVOT,failed) \
      {\
	 int NN=N; \
	 int res =0; \
	 int ONE=1; \
	 dgetrs("N", &NN, &ONE, MAT, &NN, PIVOT, X, &NN, &res); \
	 if (res!=0) *failed=1;\
      }
   #else
      #define Paso_BlockOps_solve_N(N, X, MAT, PIVOT,failed) \
      {\
      int res=clapack_dgetrs(CblasColMajor, CblasNoTrans, N,1,MAT,N, PIVOT, X, N); \
	 if (res!=0) *failed=1;\
      }
   #endif
#else
   #define Paso_BlockOps_solve_N(N, X, MAT, PIVOT,failed) PASO_MISSING_CLAPACK
#endif
   
 
/* inplace matrix vector product */

#define Paso_BlockOps_MViP_2(MAT, V) \
{ \
   register double S1=*(V+0);\
   register double S2=*(V+1);\
   register double A11=*(MAT+0);\
   register double A12=*(MAT+2);\
   register double A21=*(MAT+1);\
   register double A22=*(MAT+3);\
   *(V+0)  = A11 * S1 + A12 * S2;\
   *(V+1)  = A21 * S1 + A22 * S2;\
}


#define Paso_BlockOps_MViP_3(MAT, V) \
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
   *(V+0)=A11 * S1 + A12 * S2 + A13 * S3; \
   *(V+1)=A21 * S1 + A22 * S2 + A23 * S3;\
   *(V+2)=A31 * S1 + A32 * S2 + A33 * S3;\
}

	
#endif /* #INC_Paso_BlockOpsOPS */
