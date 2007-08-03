/* $Id:$ */

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

/* Paso: SparseMatrix :                           */
/*  copies a SparseMatrix values into an output array */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005, 2006, 2007 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"


void Paso_SparseMatrix_copy(Paso_SparseMatrix* in,double* array) {
  dim_t i,j;
  index_t iptr;
  index_t index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  #pragma omp parallel for private(i,iptr,j) schedule(static)
  for (i=0;i< in->pattern->numOutput;++i) {
     for (iptr=(in->pattern->ptr[i])-index_offset;iptr<(in->pattern->ptr[i+1])-index_offset; ++iptr) {
         for (j=0;j<in->block_size;j++) array[iptr*(in->block_size)+j]=in->val[iptr*(in->block_size)+j];
     }
  }
}
