
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


/**************************************************************/

/* Paso: SparseMatrix :                           */
/*  sets the values of the sparse matrix to a value */

/**************************************************************/

/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

/**************************************************************/

void  Paso_SparseMatrix_setValues(Paso_SparseMatrix* in,double value) {
  index_t index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  dim_t i,j;
  index_t iptr;
  if (! Paso_Pattern_isEmpty(in->pattern)) {
     #pragma omp parallel for private(i,iptr,j) schedule(static)
     for (i=0;i< in->pattern->numOutput;++i) {
        for (iptr=(in->pattern->ptr[i])-index_offset;iptr<(in->pattern->ptr[i+1])-index_offset;++iptr) {
            for (j=0;j<(in->block_size);++j) in->val[iptr*(in->block_size)+j]=value;
        }
     }
 }
}
