
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

/************************************************************************************/

/*   Some utility routines: */

/************************************************************************************/

#ifndef INC_DUDLEY_UTIL
#define INC_DUDLEY_UTIL

#include "Dudley.h"

/************************************************************************************/

void Dudley_Util_Gather_double(dim_t len, index_t * index, dim_t numData, double *in, double *out);
void Dudley_Util_Gather_int(dim_t len, index_t * index, dim_t numData, index_t * in, index_t * out);
void Dudley_Util_AddScatter(const dim_t len, const index_t * index, const dim_t numData, const double *in, double *out, const index_t upperBound);
void Dudley_Util_SmallMatMult(dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C);
void Dudley_Util_SmallMatSetMult(dim_t len, dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C);
void Dudley_Util_SmallMatSetMult1(dim_t len, dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C);
void Dudley_Util_InvertSmallMat(dim_t len, dim_t dim, double *A, double *invA, double *det);
void Dudley_Util_DetOfSmallMat(dim_t len, dim_t dim, double *A, double *det);
void Dudley_NormalVector(dim_t len, dim_t dim, dim_t dim1, double *A, double *Normal);
void Dudley_LengthOfNormalVector(dim_t len, dim_t dim, dim_t dim1, double *A, double *length);
void Dudley_Util_InvertMap(dim_t, index_t *, dim_t, index_t *);
index_t Dudley_Util_getMaxInt(dim_t dim, dim_t N, index_t * values);
index_t Dudley_Util_getMinInt(dim_t dim, dim_t N, index_t * values);
index_t Dudley_Util_getFlaggedMaxInt(dim_t dim, dim_t N, index_t * values, index_t ignore);
index_t Dudley_Util_getFlaggedMinInt(dim_t dim, dim_t N, index_t * values, index_t ignore);
dim_t Dudley_Util_packMask(dim_t N, bool_t * mask, index_t * index);
bool_t Dudley_Util_isAny(dim_t N, index_t * array, index_t value);
index_t Dudley_Util_cumsum(dim_t, index_t *);
bool_t Dudley_Util_anyNonZeroDouble(dim_t N, double *values);
void Dudley_Util_setValuesInUse(const index_t * values, const dim_t numValues, dim_t * numValuesInUse,
				index_t ** valuesInUse, Esys_MPIInfo * mpiinfo);

#ifdef ESYS_MPI
void Dudley_printDoubleArray(FILE * fid, dim_t n, double *array, char *name);
void Dudley_printIntArray(FILE * fid, dim_t n, int *array, char *name);
void Dudley_printMaskArray(FILE * fid, dim_t n, int *array, char *name);
#endif

/* Dudley_Util_orderValueAndIndex is used to sort items by a value */
/* index points to the location of the original item array. */
/* it can be used to reorder the array */
struct Dudley_Util_ValueAndIndex {
    index_t index;
    index_t value;
};
typedef struct Dudley_Util_ValueAndIndex Dudley_Util_ValueAndIndex;

void Dudley_Util_sortValueAndIndex(dim_t n, Dudley_Util_ValueAndIndex * array);
int Dudley_Util_ValueAndIndex_compar(const void *, const void *);

#endif				/* #ifndef INC_UTIL_UTIL */
