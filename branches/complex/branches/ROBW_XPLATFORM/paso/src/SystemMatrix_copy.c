/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrix :                           */
/*  copies a SystemMatrix values into an output array */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"


void Paso_SystemMatrix_copy(Paso_SystemMatrix* in,double* array) {
  dim_t i,j;
  index_t iptr;
  index_t index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  #pragma omp parallel for private(i,iptr,j) schedule(static)
  for (i=0;i< in->pattern->n_ptr;++i) {
     for (iptr=(in->pattern->ptr[i])-index_offset;iptr<(in->pattern->ptr[i+1])-index_offset; ++iptr) {
         for (j=0;j<in->block_size;j++) array[iptr*(in->block_size)+j]=in->val[iptr*(in->block_size)+j];
     }
  }
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
