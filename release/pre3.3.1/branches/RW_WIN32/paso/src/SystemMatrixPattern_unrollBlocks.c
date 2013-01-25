/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern, \
                                           dim_t row_block_size,dim_t col_block_size) {
  Paso_SystemMatrixPattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,iPtr;
  dim_t i,j,k,l;
  Paso_resetError();
  dim_t block_size=row_block_size*col_block_size;
  dim_t new_n_ptr=(pattern->n_ptr)*row_block_size;
  dim_t new_len=(pattern->len)*block_size;

  ptr=MEMALLOC(new_n_ptr+1,index_t);
  index=MEMALLOC(new_len,index_t);


  if (! ( Paso_checkPtr(ptr) || Paso_checkPtr(index) ) )  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<new_n_ptr+1;++i) ptr[i]=0;

        #pragma omp master
        ptr[new_n_ptr]=new_len;

        #pragma omp for private(i,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) 
            for (k=0;k<row_block_size;++k) ptr[i*row_block_size+k]=(pattern->ptr[i]-PTR_OFFSET)*block_size+(pattern->ptr[i+1]-pattern->ptr[i])*col_block_size*k;
          
        #pragma omp for private(i,iPtr) schedule(static) 
        for (i=0;i<new_n_ptr;++i) 
            for (iPtr=ptr[i];iPtr<ptr[i+1];++iPtr) index[iPtr]=0;

        #pragma omp for private(i,j,iPtr,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) {
           for (iPtr=pattern->ptr[i];iPtr<pattern->ptr[i+1];++iPtr)  {
              for (k=0;k<row_block_size;++k) {
                 for (j=0;j<col_block_size;++j) {
                    index[ptr[i*row_block_size+k]+(iPtr-pattern->ptr[i])*col_block_size+j]=(pattern->index[iPtr]-INDEX_OFFSET)*col_block_size+j;
                 }
              }
           }
        }
     }
     /* create return value */
     out=Paso_SystemMatrixPattern_alloc(new_n_ptr,ptr,index);
  }  
  if (! Paso_noError()) {
     MEMFREE(index);
     MEMFREE(ptr);
  }
  return out;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
