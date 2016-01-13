
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

void Paso_Pattern_cop(Paso_SparseMatrix* A, index_t* mis_marker) {

  index_t index_offset=(A->pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  dim_t i;
  double threshold=0.05;
  index_t naib,iptr;
  bool_t flag;
  dim_t n=A->pattern->numOutput;
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_mis: symmetric matrix pattern is not supported yet");
    return;
  }
   
     /* is there any vertex available ?*/
     while (Paso_Util_isAny(n,mis_marker,IS_AVAILABLE)) {

           #pragma omp parallel for private(naib,i,iptr,flag) schedule(static) 
           for (i=0;i<n;++i) {
              if (mis_marker[i]==IS_AVAILABLE) {
                 flag=IS_IN_MIS_NOW;
                 for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
                     naib=A->pattern->index[iptr]-index_offset;
                     if (naib!=i && A->val[naib]<threshold*A->val[i]) {
                        flag=IS_AVAILABLE;
                        break;
                     }
                 }
                 mis_marker[i]=flag;
              }
           }

           #pragma omp parallel for private(naib,i,iptr) schedule(static)
           for (i=0;i<n;i++) {
              if (mis_marker[i]==IS_AVAILABLE) {
                 for (iptr=A->pattern->ptr[i]-index_offset;iptr<A->pattern->ptr[i+1]-index_offset; ++iptr) {
                     naib=A->pattern->index[iptr]-index_offset;
                     if (naib!=i && mis_marker[i]==IS_IN_MIS_NOW && A->val[naib]/A->val[i]>=-threshold)
                         mis_marker[naib]=IS_IN_MIS;
                 }
                 mis_marker[i]=IS_CONNECTED_TO_MIS;
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
