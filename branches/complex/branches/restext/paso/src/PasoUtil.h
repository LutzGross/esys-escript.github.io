
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_PASO_UTIL
#define INC_PASO_UTIL

/**************************************************************/

/*   Some utility routines: */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
/*   author: l.gross@uq.edu.au */

/**************************************************************/

#include "Common.h"
#include "Paso_MPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

index_t Paso_Util_cumsum(dim_t,index_t*);
bool_t Paso_Util_isAny(dim_t N,index_t* array,index_t value);
void Paso_zeroes(const dim_t n, double* x);
void Paso_Update(const dim_t n, const double a, double* x, const double b, const double* y);
void Paso_LinearCombination(const dim_t n, double*z, const double a,const double* x, const double b, const double* y);
double Paso_InnerProduct(const dim_t n,const double* x, const double* y, Paso_MPIInfo* mpiinfo);
double Paso_l2(const dim_t n, const double* x, Paso_MPIInfo* mpiinfo);
void ApplyGivensRotations(const dim_t n,double* v,const double* c,const double* s);
void Paso_Copy(const dim_t n, double* out, const double* in);
bool_t Paso_fileExists( const char* filename );
double Paso_lsup(const dim_t n, const double* x, Paso_MPIInfo* mpiinfo);

#define Paso_copyShortDouble(n, source, target)  memcpy(target,source,sizeof(double)*(size_t)n)

#endif /* #ifndef INC_PASO_UTIL */
