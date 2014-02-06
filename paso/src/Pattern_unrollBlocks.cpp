
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/* Paso: Pattern_unrollBlocks */

/************************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "Pattern.h"

/************************************************************************************/

/* creates Pattern  */

Paso_Pattern* Paso_Pattern_unrollBlocks(Paso_Pattern* pattern, \
                                        int type, dim_t output_block_size,dim_t input_block_size) {
  Paso_Pattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,iPtr;
  dim_t i,j,k, block_size, new_len, new_numOutput, new_numInput;
  const index_t index_offset_in=(pattern->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  const index_t index_offset_out=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);

  
  Esys_resetError();
  if ( ( output_block_size == 1 ) && (input_block_size == 1) && ((pattern->type & MATRIX_FORMAT_OFFSET1) == (type & MATRIX_FORMAT_OFFSET1) ) ) {
     out = Paso_Pattern_getReference(pattern);
  } else {
     block_size=output_block_size*input_block_size;
     new_len=(pattern->len)*block_size;
     new_numOutput=(pattern->numOutput)*output_block_size;
     new_numInput=(pattern->numInput)*input_block_size;
   
     ptr=new index_t[new_numOutput+1];
     index=new index_t[new_len];
     if (! ( Esys_checkPtr(ptr) || Esys_checkPtr(index) ) )  {
        #pragma omp parallel
        {
           #pragma omp for private(i) schedule(static)
           for (i=0;i<new_numOutput+1;++i) ptr[i]=index_offset_out;
   
           #pragma omp single
           ptr[new_numOutput]=new_len+index_offset_out;
   
           #pragma omp for private(i,k) schedule(static) 
           for (i=0;i<pattern->numOutput;++i) {
               for (k=0;k<output_block_size;++k) {
                    ptr[i*output_block_size+k]=(pattern->ptr[i]-index_offset_in)*block_size+
                                                       (pattern->ptr[i+1]-pattern->ptr[i])*input_block_size*k+index_offset_out;
               }
           }
             
           #pragma omp for private(i,iPtr) schedule(static) 
           for (i=0;i<new_numOutput;++i) {
               #pragma ivdep
               for (iPtr=ptr[i]-index_offset_out;iPtr<ptr[i+1]-index_offset_out;++iPtr) index[iPtr]=index_offset_out;
   	   }
   
           #pragma omp for private(i,j,iPtr,k) schedule(static) 
           for (i=0;i<pattern->numOutput;++i) {
              for (iPtr=pattern->ptr[i]-index_offset_in;iPtr<pattern->ptr[i+1]-index_offset_in;++iPtr)  {
                 for (k=0;k<output_block_size;++k) {
                    #pragma ivdep
                    for (j=0;j<input_block_size;++j) {
                       index[ptr[i*output_block_size+k]-index_offset_out+(iPtr-(pattern->ptr[i]-index_offset_in))*input_block_size+j]=(pattern->index[iPtr]-index_offset_in)*input_block_size+j+index_offset_out;
                    }
                 }
              }
           }
        }
        out=Paso_Pattern_alloc(type,new_numOutput,new_numInput,ptr,index);
     }  
     if (! Esys_noError()) {
        delete[] index;
        delete[] ptr;
     }
  }
  return out;
}
