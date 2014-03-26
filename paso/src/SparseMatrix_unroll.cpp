
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


/****************************************************************************/

/* Paso: SparseMatrix */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

namespace paso {

SparseMatrix* SparseMatrix_unroll(SparseMatrixType type, const SparseMatrix* A) {

  const dim_t row_block_size = A->row_block_size;
  const dim_t col_block_size = A->col_block_size;
  const dim_t block_size=A->block_size;
  const dim_t n = A->numRows;
  const index_t out_type = (type & MATRIX_FORMAT_BLK1) ? type : type + MATRIX_FORMAT_BLK1;
  const index_t A_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  const index_t out_offset= (out_type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t iptr, *start_p, *where_p;
  dim_t i, j, irb, icb, icol, irow, l_col, l_row;
  
  SparseMatrix*out=NULL;
  
  out = SparseMatrix_alloc(out_type, A->pattern, row_block_size, col_block_size, FALSE);
  
  if (Esys_noError()) {
     if (out->type & MATRIX_FORMAT_CSC) {
	
	#pragma omp parallel for private(i, iptr, j, irb, irow, start_p, l_col, icb, icol, where_p) schedule(static)
	for (i=0;i<n;++i) {
	   for (iptr=A->pattern->ptr[i]-A_offset; iptr<A->pattern->ptr[i+1]-A_offset; ++iptr) {
	      j=A->pattern->index[iptr]-A_offset;
	      for (icb=0; icb<col_block_size; ++icb) {
		 icol=j*col_block_size+icb;
		 start_p=&(out->pattern->index[out->pattern->ptr[icol]-out_offset]);
		 l_col=out->pattern->ptr[icol + 1]-out->pattern->ptr[icol];
		 for (irb=0; irb<row_block_size; ++irb) {
		    irow=row_block_size*i+irb+out_offset;
		    where_p=(index_t*)bsearch(&(irow), start_p, l_col,sizeof(index_t), comparIndex);		       
		    if (! (where_p == NULL) ) 
		       out->val[out->pattern->ptr[icol]-out_offset+(index_t)(where_p-start_p)] =
		       A->val[block_size*iptr+irb+row_block_size*icb];
		 }
	      }
	   }
	}
     } else {
	 #pragma omp parallel for private(i, iptr, j, irb, irow, start_p, l_row, icb, icol, where_p) schedule(static)
	 for (i=0;i<n;++i) {
	    for (iptr=A->pattern->ptr[i]-A_offset; iptr<A->pattern->ptr[i+1]-A_offset; ++iptr) {
	       j=A->pattern->index[iptr]-A_offset;
	       for (irb=0; irb<row_block_size; ++irb) {
		  irow=row_block_size*i+irb;
		  start_p=&(out->pattern->index[out->pattern->ptr[irow]-out_offset]);
		  l_row=out->pattern->ptr[irow + 1]-out->pattern->ptr[irow];
		  for (icb=0; icb<col_block_size; ++icb) {
		     icol=j*col_block_size+icb+out_offset;
		     where_p=(index_t*)bsearch(&(icol), start_p, l_row,sizeof(index_t), comparIndex);		       
		     if (! (where_p == NULL) ) 
		         out->val[out->pattern->ptr[irow]-out_offset+(index_t)(where_p-start_p)] =
							 A->val[block_size*iptr+irb+row_block_size*icb];
		  }
	       }
	    }
	 }
    }
  }
 return out;
}

} // namespace paso

