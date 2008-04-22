
/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/***************************************************************************/

/* Paso: SparseMatrix:  adds the row entries to an array */

/***************************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

void Paso_SparseMatrix_addRow_CSR_OFFSET0(Paso_SparseMatrix* A, double* array) {
   dim_t ir,irow,icb,irb;
   index_t iptr;
   register double fac;
   #pragma omp parallel for private(ir,irb,irow,fac,iptr,icb) schedule(static)
   for (ir=0;ir< A->pattern->numOutput;ir++) {
       for (irb=0;irb< A->row_block_size;irb++) {
	  irow=irb+A->row_block_size*ir;
          fac=0.;
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	      for (icb=0;icb< A->col_block_size;icb++) 
                 fac+=A->val[iptr*A->block_size+irb+A->row_block_size*icb];
          }
          array[irow]+=fac;
        }
   }
}
