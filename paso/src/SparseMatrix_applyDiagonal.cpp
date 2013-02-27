
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/* Paso: SparseMatrix                  

applies diagonal matrices from the left and the right 


************************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"



void Paso_SparseMatrix_applyDiagonal_CSR_OFFSET0(Paso_SparseMatrix* A, const double* left, const double* right) {
  index_t ir,icol,iptr,icb,irb,irow,l;
  const dim_t row_block=A->row_block_size;
  const dim_t col_block=A->col_block_size;
  const dim_t n_block=row_block*col_block;
  
  register double rtmp;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb, rtmp) schedule(static)
  for (ir=0;ir< A->pattern->numOutput;ir++) {
     for (irb=0;irb< row_block;irb++) {
	irow=irb+row_block*ir;
	rtmp=left[irow];
	for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
            #pragma ivdep
	    for (icb=0;icb< A->col_block_size;icb++) {
	       icol=icb+col_block*(A->pattern->index[iptr]);
	       l=iptr*n_block+irb+row_block*icb;
	       A->val[l]*=rtmp*right[icol];
	    }
        }
      }
   }
}
