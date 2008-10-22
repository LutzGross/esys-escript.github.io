
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**********************************************************************/

/* Paso: Pattern: Paso_Pattern_coupling 

   searches for a maximal independent set MIS in the matrix pattern 
   vertices in the maximal independent set are marked in mis_marker
   nodes to be considered are marked by -1 on the input in mis_marker

*/
/**********************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: artak@uq.edu.au                                */

/**************************************************************/

#include "mpi_C.h"
#include "Paso.h"
#include "PasoUtil.h"
#include "Pattern.h"
#include "Solver.h"


/***************************************************************/
 
#define IS_AVAILABLE -1
#define IS_IN_MIS_NOW -2
#define IS_IN_MIS -3
#define IS_CONNECTED_TO_MIS -4



void Paso_Pattern_coup(Paso_SparseMatrix* A, index_t* mis_marker) {

  index_t index_offset=(A->pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  dim_t i,j;
  double threshold=0.05;
  index_t iptr,*index,*where_p,diagptr;
  bool_t flag;
  dim_t n=A->pattern->numOutput;
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_mis: symmetric matrix pattern is not supported yet");
    return;
  }
   
     /* is there any vertex available ?*/
     while (Paso_Util_isAny(n,mis_marker,IS_AVAILABLE)) {

           #pragma omp parallel for private(i,iptr,flag) schedule(static) 
           for (i=0;i<n;++i) {
              if (mis_marker[i]==IS_AVAILABLE) {
                 flag=IS_IN_MIS;
                 diagptr=A->pattern->ptr[i];
                 index=&(A->pattern->index[diagptr]);
                 where_p=(index_t*)bsearch(&i,
                                        index,
                                        A->pattern->ptr[i + 1]-A->pattern->ptr[i],
                                        sizeof(index_t),
                                        Paso_comparIndex);
                if (where_p==NULL) {
                    Paso_setError(VALUE_ERROR, "Paso_Solver_getAMG: main diagonal element missing.");
                } else {
                    diagptr+=(index_t)(where_p-index);
                }
                 for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
                     j=A->pattern->index[iptr]-index_offset;
                     if (j!=i && A->val[iptr]>=threshold*A->val[diagptr]) {
                        flag=IS_AVAILABLE;
                        break;
                     }
                 }
                 mis_marker[i]=flag;
                }
            }
           
              #pragma omp parallel for private(i,iptr) schedule(static)
              for (i=0;i<n;i++) {
               if (mis_marker[i]==IS_AVAILABLE) {
                 diagptr=A->pattern->ptr[i];
                 index=&(A->pattern->index[diagptr]);
                 where_p=(index_t*)bsearch(&i,
                                        index,
                                        A->pattern->ptr[i + 1]-A->pattern->ptr[i],
                                        sizeof(index_t),
                                        Paso_comparIndex);
                if (where_p==NULL) {
                    Paso_setError(VALUE_ERROR, "Paso_Solver_getAMG: main diagonal element missing.");
                } else {
                    diagptr+=(index_t)(where_p-index);
                }
                 for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
                     j=A->pattern->index[iptr]-index_offset;
                     if (j!=i && mis_marker[j]==IS_IN_MIS && (A->val[iptr]/A->val[diagptr])>=-threshold){
                         mis_marker[i]=IS_IN_MIS;
                     }
                     else {
                         mis_marker[i]=IS_CONNECTED_TO_MIS;
                     }
                 }
               }
              } 
        }
     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_MIS);
}

/*
 *
 * Return a strength of connection mask using the classical 
 * strength of connection measure by Ruge and Stuben.
 *
 * Specifically, an off-diagonal entry A[i.j] is a strong 
 * connection if:
 *  
 *      -A[i,j] >= theta * max( -A[i,k] )   where k != i
 * 
 * Otherwise, the connection is weak.
 *  
 */  
void Paso_Pattern_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta) 
{
  index_t index_offset=(A->pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  dim_t i,j;
  index_t iptr;
  double threshold,min_offdiagonal;
  bool_t flag;
  dim_t n=A->pattern->numOutput;
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_RS: symmetric matrix pattern is not supported yet");
    return;
  }

/* is there any vertex available ?*/
     if (Paso_Util_isAny(n,mis_marker,IS_AVAILABLE)) {

     #pragma omp parallel for private(i,iptr,min_offdiagonal,threshold) schedule(static) 
      for (i=0;i<n;++i) {
        min_offdiagonal = A->val[A->pattern->ptr[i]-index_offset];
        for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
            if(A->pattern->index[iptr] != i){
                min_offdiagonal = MIN(min_offdiagonal,A->val[iptr-index_offset]);
            }
        }

        threshold = theta*min_offdiagonal;
        mis_marker[i]=IS_CONNECTED_TO_MIS;
        #pragma omp parallel for private(iptr) schedule(static) 
        for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
            if(-1.*(A->val[iptr-index_offset]) <= threshold){
                    mis_marker[i]=IS_IN_MIS;
            }
        }
      }
    }
     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_MIS);
}
#undef IS_AVAILABLE 
#undef IS_IN_MIS_NOW 
#undef IS_IN_MIS 
#undef IS_CONNECTED_TO_MIS 
