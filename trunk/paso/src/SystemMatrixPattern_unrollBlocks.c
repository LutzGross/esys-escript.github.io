/* $Id$ */

/*
********************************************************************************
*               Copyright © 2006 by ACcESS MNRF                                *
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

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern, \
                                           int type, dim_t row_block_size,dim_t col_block_size) {
  index_t index_offset_in=(pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t index_offset_out=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  
  if ((pattern->type & PATTERN_FORMAT_SYM) != (type & PATTERN_FORMAT_SYM)) {
      Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_unrollBlocks: conversion between symmetric and non-symmetric is not implemented yet");
      return NULL;
  }
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
        for (i=0;i<new_n_ptr+1;++i) ptr[i]=index_offset_out;

        #pragma omp master
        ptr[new_n_ptr]=new_len+index_offset_out;

        #pragma omp for private(i,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) 
            for (k=0;k<row_block_size;++k) ptr[i*row_block_size+k]=(pattern->ptr[i]-index_offset_in)*block_size+(pattern->ptr[i+1]-pattern->ptr[i])*col_block_size*k+index_offset_out;
          
        #pragma omp for private(i,iPtr) schedule(static) 
        for (i=0;i<new_n_ptr;++i) 
            for (iPtr=ptr[i]-index_offset_out;iPtr<ptr[i+1]-index_offset_out;++iPtr) index[iPtr]=index_offset_out;

        #pragma omp for private(i,j,iPtr,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) {
           for (iPtr=pattern->ptr[i]-index_offset_in;iPtr<pattern->ptr[i+1]-index_offset_in;++iPtr)  {
              for (k=0;k<row_block_size;++k) {
                 for (j=0;j<col_block_size;++j) {
                    index[ptr[i*row_block_size+k]-index_offset_out+(iPtr-(pattern->ptr[i]-index_offset_in))*col_block_size+j]=(pattern->index[iPtr]-index_offset_in)*col_block_size+j+index_offset_out;
                 }
              }
           }
        }
     }
     /* create return value */
     out=Paso_SystemMatrixPattern_alloc(type,new_n_ptr,ptr,index);
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
