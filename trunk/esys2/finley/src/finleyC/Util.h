/* $Id$ */

#ifndef INC_FINLEY_UTIL
#define INC_FINLEY_UTIL

/**************************************************************/

/*   Some utility routines: */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"

/**************************************************************/

void Finley_Util_Gather_double(int len,maybelong* index,int numData,double* in,double * out);
void Finley_Util_Gather_int(int len,maybelong* index,int numData,maybelong* in,maybelong * out);
void Finley_Util_AddScatter(int len,maybelong* index,int numData,double* in,double * out);
void Finley_Util_SmallMatMult(int A1,int A2, double* A, int B2, double*B, double* C);
void Finley_Util_SmallMatSetMult(int len,int A1,int A2, double* A, int B2, double*B, double* C);
void Finley_Util_InvertSmallMat(int len,int dim,double* A,double *invA, double* det);
void Finley_Util_DetOfSmallMat(int len,int dim,double* A,double* det);
void Finley_NormalVector(int len, int dim, int dim1, double* A,double* Normal);
void Finley_LengthOfNormalVector(int len, int dim, int dim1, double* A,double* length);
void Finley_Util_InvertMap(int, maybelong*,int, maybelong*);
maybelong Finley_Util_getMaxInt(int dim,int N,maybelong* values);
maybelong Finley_Util_getMinInt(int dim,int N,maybelong* values);
maybelong Finley_Util_packMask(maybelong N,maybelong* mask,maybelong* index);
int Finley_Util_isAny(maybelong N,maybelong* array,maybelong value);
void Finley_copyDouble(int n,double* source,double* target);


/* Finley_Util_orderValueAndIndex is used to sort items by a value */
/* index points to the location of the original item array. */
/* it can be used to reorder the array */
struct Finley_Util_ValueAndIndex {
   maybelong index;
   maybelong value;
};
typedef struct Finley_Util_ValueAndIndex Finley_Util_ValueAndIndex;

void Finley_Util_sortValueAndIndex(int n,Finley_Util_ValueAndIndex* array);

int Finley_Util_ValueAndIndex_compar(const void *, const void *);

#endif /* #ifndef INC_FINLEY_UTIL */

/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:47  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.3  2004/08/26 12:03:52  gross
 * Some other bug in Finley_Assemble_gradient fixed.
 *
 * Revision 1.2  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
