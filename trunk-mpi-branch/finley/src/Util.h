/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************

*/

/**************************************************************/

/*   Some utility routines: */

/**************************************************************/

/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_UTIL
#define INC_FINLEY_UTIL

#include "Finley.h"

/**************************************************************/

void Finley_Util_Gather_double(dim_t len,index_t* index,dim_t numData,double* in,double * out);
void Finley_Util_Gather_int(dim_t len,index_t* index,dim_t numData,index_t* in,index_t * out);
void Finley_Util_AddScatter(dim_t len,index_t* index,dim_t numData,double* in,double * out, index_t bound);
void Finley_Util_SmallMatMult(dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_SmallMatSetMult(dim_t len,dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_InvertSmallMat(dim_t len,dim_t dim,double* A,double *invA, double* det);
void Finley_Util_DetOfSmallMat(dim_t len,dim_t dim,double* A,double* det);
void Finley_NormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* Normal);
void Finley_LengthOfNormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* length);
void Finley_Util_InvertMap(dim_t, index_t*,dim_t, index_t*);
index_t Finley_Util_getMaxInt(dim_t dim,dim_t N,index_t* values);
index_t Finley_Util_getMinInt(dim_t dim,dim_t N,index_t* values);
index_t Finley_Util_getFlaggedMaxInt(dim_t dim,dim_t N,index_t* values,index_t ignore);
index_t Finley_Util_getFlaggedMinInt(dim_t dim,dim_t N,index_t* values,index_t ignore);
dim_t Finley_Util_packMask(dim_t N,bool_t* mask,index_t* index);
bool_t Finley_Util_isAny(dim_t N,index_t* array,index_t value);
index_t Finley_Util_cumsum(dim_t,index_t*);
bool_t Finley_Util_anyNonZeroDouble(dim_t N,double* values);

#ifdef PASO_MPI
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

