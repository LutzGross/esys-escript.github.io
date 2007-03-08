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

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern, \
                                           int type, dim_t output_block_size,dim_t input_block_size) {
  index_t index_offset_in=(pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t index_offset_out=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  
  if ((pattern->type & PATTERN_FORMAT_SYM) != (type & PATTERN_FORMAT_SYM)) {
      Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_unrollBlocks: conversion between symmetric and non-symmetric is not implemented yet");
      return NULL;
  }
  Paso_MPIInfo* mpi_info=pattern->mpi_info;
  Paso_SystemMatrixPattern*out=NULL;
  Paso_Distribution *input_dist=NULL,*output_dist=NULL;
  index_t *ptr=NULL,*index=NULL,iPtr;
  dim_t i,j,k,l;
  Paso_resetError();
  dim_t block_size=output_block_size*input_block_size;
  dim_t new_myNumOutput=(pattern->myNumOutput)*output_block_size;
  dim_t new_myLen=(pattern->myLen)*block_size;

  ptr=MEMALLOC(new_myNumOutput+1,index_t);
  index=MEMALLOC(new_myLen,index_t);
  if (! ( Paso_checkPtr(ptr) || Paso_checkPtr(index) ) )  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<new_myNumOutput+1;++i) ptr[i]=index_offset_out;

        #pragma omp master
        ptr[new_myNumOutput]=new_myLen+index_offset_out;

        #pragma omp for private(i,k) schedule(static) 
        for (i=0;i<pattern->myNumOutput;++i) 
            for (k=0;k<output_block_size;++k) ptr[i*output_block_size+k]=(pattern->ptr[i]-index_offset_in)*block_size+(pattern->ptr[i+1]-pattern->ptr[i])*input_block_size*k+index_offset_out;
          
        #pragma omp for private(i,iPtr) schedule(static) 
        for (i=0;i<new_myNumOutput;++i) 
            for (iPtr=ptr[i]-index_offset_out;iPtr<ptr[i+1]-index_offset_out;++iPtr) index[iPtr]=index_offset_out;

        #pragma omp for private(i,j,iPtr,k) schedule(static) 
        for (i=0;i<pattern->myNumOutput;++i) {
           for (iPtr=pattern->ptr[i]-index_offset_in;iPtr<pattern->ptr[i+1]-index_offset_in;++iPtr)  {
              for (k=0;k<output_block_size;++k) {
                 for (j=0;j<input_block_size;++j) {
                    index[ptr[i*output_block_size+k]-index_offset_out+(iPtr-(pattern->ptr[i]-index_offset_in))*input_block_size+j]=(pattern->index[iPtr]-index_offset_in)*input_block_size+j+index_offset_out;
                 }
              }
           }
        }
     }
     /* create return value */
     input_dist=Paso_Distribution_alloc(mpi_info,pattern->input_distribution->first_component,
                                        input_block_size, -index_offset_in*input_block_size+index_offset_out);
     output_dist=Paso_Distribution_alloc(mpi_info,pattern->output_distribution->first_component,
                                         output_block_size, -index_offset_in*output_block_size+index_offset_out);
     out=Paso_SystemMatrixPattern_alloc(type,output_dist,input_dist,ptr,index);
  }  
  if (! Paso_noError()) {
     MEMFREE(index);
     MEMFREE(ptr);
  }
  Paso_Distribution_free(input_dist);
  Paso_Distribution_free(output_dist);
  return out;
}
