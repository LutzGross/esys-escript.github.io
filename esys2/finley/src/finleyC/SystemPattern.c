/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "Common.h"
#include "SystemPattern.h"

/**************************************************************/

/* allocates a SystemMatrixPattern  */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_alloc(int n_ptr, index_t* ptr,index_t* index) {
  Finley_SystemMatrixPattern*out;
  index_t loc_min_index,loc_max_index,min_index=INDEX_OFFSET,max_index=INDEX_OFFSET-1,k;
  dim_t i;
  Finley_ErrorCode=NO_ERROR;


  #pragma omp parallel private(loc_min_index,loc_max_index,i,k)
  {
     loc_min_index=INDEX_OFFSET;
     loc_max_index=INDEX_OFFSET-1;
     #if PTR_OFFSET>0
        #pragma omp for schedule(static)
        for (i=0;i<n_ptr+1;++i) ptr[i]+=PTR_OFFSET;
     #endif
     #if INDEX_OFFSET>0
        #pragma omp for schedule(static) 
        for (i=0;i<n_ptr;++i) 
             for (k=ptr[i];k<ptr[i+1];++k) index[k]+=INDEX_OFFSET;
     #endif
        
     #pragma omp for schedule(static) 
     for (i=0;i<n_ptr;++i) {
         if (ptr[i]<ptr[i+1]) {
           qsort(&(index[ptr[i]-PTR_OFFSET]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Finley_comparIndex); 
           loc_min_index=MIN(loc_min_index,index[ptr[i]]);
           loc_max_index=MAX(loc_max_index,index[ptr[i+1]-1]);
         }
     }
     #pragma omp critical
     {
        min_index=MIN(loc_min_index,min_index);
        max_index=MAX(loc_max_index,max_index);
     }
  }
  if (min_index<INDEX_OFFSET) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"Matrix pattern index out of range.");
    return NULL;
  }

  out=MEMALLOC(1,Finley_SystemMatrixPattern);
  if (Finley_checkPtr(out)) return NULL;
  out->n_ptr=n_ptr;
  out->n_index=max_index+1-INDEX_OFFSET;
  out->ptr=ptr;
  out->index=index;
  out->len=out->ptr[out->n_ptr];
  out->reference_counter=1;
  #ifdef Finley_TRACE
  printf("Finley_SystemMatrixPattern_dealloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_reference(Finley_SystemMatrixPattern* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}
  
/* deallocates a SystemMatrixPattern: */

void Finley_SystemMatrixPattern_dealloc(Finley_SystemMatrixPattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->ptr);
        MEMFREE(in->index);
        MEMFREE(in);
        #ifdef Finley_TRACE
        printf("Finley_SystemMatrixPattern_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}
/* *************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

/*  this routine is used by qsort called in Finley_SystemMatrixPattern_alloc */

int Finley_comparIndex(const void *index1,const void *index2){
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
 * Revision 1.5  2005/08/23 01:24:30  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-08-23
 *
 * Revision 1.4.2.2  2005/08/19 02:44:09  gross
 * stopping criterion modified to cope with badly balanced equations
 *
 * Revision 1.4.2.1  2005/08/15 12:02:53  gross
 * memory leak fixed. it is still not clear if there is no problem anymore
 *
 * Revision 1.4  2005/07/08 04:07:57  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.2.6  2005/07/01 07:02:13  gross
 * some bug with OPENMP fixed
 *
 * Revision 1.1.2.5  2005/06/30 01:53:56  gross
 * a bug in coloring fixed
 *
 * Revision 1.1.2.4  2005/06/29 02:34:55  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.3  2005/03/04 05:27:07  gross
 * bug in SystemPattern fixed.
 *
 * Revision 1.1.2.2  2004/11/24 01:37:15  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.2.1  2004/11/14 23:49:09  gross
 * the forgotten files
 *
 *
 */
