/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
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
#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************/

void Finley_Util_Gather_double(dim_t len,index_t* index,dim_t numData,double* in,double * out);
void Finley_Util_Gather_int(dim_t len,index_t* index,dim_t numData,index_t* in,index_t * out);
void Finley_Util_AddScatter(dim_t len,index_t* index,dim_t numData,double* in,double * out);
void Finley_Util_SmallMatMult(dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_SmallMatSetMult(dim_t len,dim_t A1,dim_t A2, double* A, dim_t B2, double*B, double* C);
void Finley_Util_InvertSmallMat(dim_t len,dim_t dim,double* A,double *invA, double* det);
void Finley_Util_DetOfSmallMat(dim_t len,dim_t dim,double* A,double* det);
void Finley_NormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* Normal);
void Finley_LengthOfNormalVector(dim_t len, dim_t dim, dim_t dim1, double* A,double* length);
void Finley_Util_InvertMap(dim_t, index_t*,dim_t, index_t*);
index_t Finley_Util_getMaxInt(dim_t dim,dim_t N,index_t* values);
index_t Finley_Util_getMinInt(dim_t dim,dim_t N,index_t* values);
dim_t Finley_Util_packMask(dim_t N,bool_t* mask,index_t* index);
bool_t Finley_Util_isAny(dim_t N,index_t* array,index_t value);
void Finley_copyDouble(dim_t n,double* source,double* target);
index_t Finley_Util_cumsum(dim_t,index_t*);
bool_t Finley_Util_anyNonZeroDouble(dim_t N,double* values);



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
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* #ifndef INC_FINLEY_UTIL */

/*
 * Revision 1.8  2005/08/12 01:45:43  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.7.2.2  2005/09/07 06:26:22  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.7.2.1  2005/08/04 22:41:11  gross
 * some extra routines for finley that might speed-up RHS assembling in some cases (not actived right now)
 *
 * Revision 1.7  2005/07/08 04:07:59  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.4  2005/06/29 02:34:57  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.3  2005/03/02 23:35:06  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 * Revision 1.1.1.1.2.2  2005/02/18 02:27:31  gross
 * two function that will be used for a reimplementation of the ILU preconditioner
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:19  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
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
