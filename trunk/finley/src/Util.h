
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Some utility routines

*****************************************************************************/

#ifndef INC_FINLEY_UTIL
#define INC_FINLEY_UTIL

#include "Finley.h"
#include "esysUtils/Esys_MPI.h"

#include <vector>

void Finley_Util_Gather_double(dim_t len,index_t* index,dim_t numData,double* in,double * out);
void Finley_Util_Gather_int(dim_t len,index_t* index,dim_t numData,index_t* in,index_t * out);
void Finley_Util_AddScatter(const dim_t len, const index_t* index, const dim_t numData, const double* in,double * out, const index_t upperBound);
void Finley_Util_SmallMatMult(dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_SmallMatSetMult(dim_t len,dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_SmallMatSetMult1(dim_t len,dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_InvertSmallMat(dim_t len,dim_t dim,double* A,double *invA, double* det);
void Finley_Util_DetOfSmallMat(dim_t len,dim_t dim,double* A,double* det);
void Finley_NormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* Normal);
void Finley_LengthOfNormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* length);
void Finley_Util_InvertMap(dim_t, index_t*,dim_t, index_t*);
index_t Finley_Util_getMaxInt(int dim, int N, const int* values);
index_t Finley_Util_getMinInt(int dim, int N, const int* values);
std::pair<int,int> Finley_Util_getMinMaxInt(int dim, int N, const int* values);
index_t Finley_Util_getFlaggedMaxInt(dim_t dim,dim_t N,index_t* values,index_t ignore);
index_t Finley_Util_getFlaggedMinInt(dim_t dim,dim_t N,index_t* values,index_t ignore);
dim_t Finley_Util_packMask(dim_t N,bool_t* mask,index_t* index);
bool_t Finley_Util_isAny(dim_t N,index_t* array,index_t value);
index_t Finley_Util_cumsum(dim_t,index_t*);
bool_t Finley_Util_anyNonZeroDouble(dim_t N,double* values);
void Finley_Util_setValuesInUse(const int *values, const int numValues,
                                std::vector<int>& valuesInUse,
                                Esys_MPIInfo* mpiinfo);

#ifdef ESYS_MPI
void Finley_printDoubleArray( FILE *fid, dim_t n, double *array, char *name  );
void Finley_printIntArray( FILE *fid, dim_t n, int *array, char *name  );
void Finley_printMaskArray( FILE *fid, dim_t n, int *array, char *name  );
#endif




/* Finley_Util_orderValueAndIndex is used to sort items by a value */
/* index points to the location of the original item array. */
/* it can be used to reorder the array */
struct Finley_Util_ValueAndIndex {
   index_t index;
   index_t value;
};
typedef struct Finley_Util_ValueAndIndex Finley_Util_ValueAndIndex;

void Finley_Util_sortValueAndIndex(dim_t n,Finley_Util_ValueAndIndex* array);
int Finley_Util_ValueAndIndex_compar(const void *, const void *);

#endif /* #ifndef INC_FINLEY_UTIL */

