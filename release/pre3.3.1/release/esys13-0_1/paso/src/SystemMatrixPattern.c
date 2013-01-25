/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* allocates a SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type, int n_ptr, index_t* ptr,index_t* index) {
  Paso_SystemMatrixPattern*out;
  index_t index_offset=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t loc_min_index,loc_max_index,min_index=index_offset,max_index=index_offset-1;
  dim_t i;
  Paso_resetError();

  if (type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"symmetric matrix pattern is not supported yet");
    return NULL;
  }
  #pragma omp parallel private(loc_min_index,loc_max_index,i)
   {
      loc_min_index=index_offset;
      loc_max_index=index_offset-1;
      if (type & PATTERN_FORMAT_OFFSET1) {
         #pragma omp for schedule(static) 
         for (i=0;i<n_ptr;++i) {
             if (ptr[i]<ptr[i+1]) {
               qsort(&(index[ptr[i]-1]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
               loc_min_index=MIN(loc_min_index,index[ptr[i]-1]);
               loc_max_index=MAX(loc_max_index,index[ptr[i+1]-2]);
             }
         }
      } else {
         #pragma omp for schedule(static) 
         for (i=0;i<n_ptr;++i) {
             if (ptr[i]<ptr[i+1]) {
               qsort(&(index[ptr[i]]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
               loc_min_index=MIN(loc_min_index,index[ptr[i]]);
               loc_max_index=MAX(loc_max_index,index[ptr[i+1]-1]);
             }
         }
      }
      #pragma omp critical
      {
         min_index=MIN(loc_min_index,min_index);
         max_index=MAX(loc_max_index,max_index);
      }
  }
  if (min_index<index_offset) {
    Paso_setError(TYPE_ERROR,"Matrix pattern index out of range.");
    return NULL;
  }

  out=MEMALLOC(1,Paso_SystemMatrixPattern);
  if (Paso_checkPtr(out)) return NULL;
  out->n_ptr=n_ptr;
  out->n_index=max_index-index_offset+1;
  out->ptr=ptr;
  out->index=index;
  out->len=out->ptr[out->n_ptr]-index_offset;
  out->reference_counter=1;
  out->type=type;
  #ifdef Paso_TRACE
  printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a SystemMatrixPattern: */

void Paso_SystemMatrixPattern_dealloc(Paso_SystemMatrixPattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->ptr);
        MEMFREE(in->index);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}
/* *************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

/*  this routine is used by qsort called in Paso_SystemMatrixPattern_alloc */

int Paso_comparIndex(const void *index1,const void *index2){
   index_t Iindex1,Iindex2;
   Iindex1=*(index_t*)index1;
   Iindex2=*(index_t*)index2;
   if (Iindex1<Iindex2) {
      return -1;
   } else {
      if (Iindex1>Iindex2) {
         return 1;
      } else {
         return 0;
      }
   }
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
