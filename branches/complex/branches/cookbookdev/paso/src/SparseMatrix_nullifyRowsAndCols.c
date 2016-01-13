
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SparseMatrix                                       */

/*  nullify rows and columns in the matrix                    */

/*  the rows and columns are marked by positive values in     */
/*  mask_row and mask_col. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row and      */
/*  mask_col are set to main_diagonal_value                   */


/**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

void Paso_SparseMatrix_nullifyRowsAndCols_CSC_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol;
  #pragma omp parallel for private(irow, iptr,icol) schedule(static)
  for (icol=0;icol< A->pattern->numOutput;icol++) {
     #pragma ivdep
     for (iptr=A->pattern->ptr[icol]-index_offset;iptr<A->pattern->ptr[icol+1]-index_offset; iptr++) {
	  irow=A->pattern->index[iptr]-index_offset;
	  if (mask_col[icol]>0. || mask_row[irow]>0. ) {
            if (irow==icol) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	  }
	}
      }
}
void Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol;
  #pragma omp parallel for private(irow, iptr,icol) schedule(static)
  for (irow=0;irow< A->pattern->numOutput;irow++) {
      /* TODO: test mask_row here amd not inside every loop */
      #pragma ivdep
      for (iptr=A->pattern->ptr[irow]-index_offset;iptr<A->pattern->ptr[irow+1]-index_offset; iptr++) {
        icol=A->pattern->index[iptr]-index_offset;
        if (mask_col[icol]>0. || mask_row[irow]>0. ) {
           if (irow==icol) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	}
     }
  } 
}
void Paso_SparseMatrix_nullifyRowsAndCols_CSC(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t icol,iptr,icb,irb,irow,ic,l;
  #pragma omp parallel for private(l,irow, iptr,icol,ic,irb,icb) schedule(static)
  for (ic=0;ic< A->pattern->numOutput;ic++) {
	for (iptr=A->pattern->ptr[ic]-index_offset;iptr<A->pattern->ptr[ic+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*(A->pattern->index[iptr]-index_offset);
            #pragma ivdep
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*ic;
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
                l=iptr*A->block_size+irb+A->row_block_size*icb;
		if (irow==icol) {
		  A->val[l]=main_diagonal_value;
		} else {
		  A->val[l]=0;
		}
	      }
	    }
	  }
	}
  }
}
void Paso_SparseMatrix_nullifyRowsAndCols_CSR(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t ir,icol,iptr,icb,irb,irow,l;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb) schedule(static)
  for (ir=0;ir< A->pattern->numOutput;ir++) {
	for (iptr=A->pattern->ptr[ir]-index_offset;iptr<A->pattern->ptr[ir+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
            #pragma ivdep
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*(A->pattern->index[iptr]-index_offset);
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
                l=iptr*A->block_size+irb+A->row_block_size*icb;
		if (irow==icol) {
		  A->val[l]=main_diagonal_value;
		} else {
		  A->val[l]=0;
		}
	      }
	    }
	  }
	}
  }
}
void Paso_SparseMatrix_nullifyRows_CSR_BLK1(Paso_SparseMatrix* A, double* mask_row, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol;
  #pragma omp parallel for private(irow, iptr,icol) schedule(static)
  for (irow=0;irow< A->pattern->numOutput;irow++) {
      if (mask_row[irow]>0.) {
         #pragma ivdep
         for (iptr=A->pattern->ptr[irow]-index_offset;iptr<A->pattern->ptr[irow+1]-index_offset; iptr++) {
           icol=A->pattern->index[iptr]-index_offset;
           if (irow==icol) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	 }
     }
  } 
}
void Paso_SparseMatrix_nullifyRows_CSR(Paso_SparseMatrix* A, double* mask_row, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t ir,icol,iptr,icb,irb,irow,l;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb) schedule(static)
  for (ir=0;ir< A->pattern->numOutput;ir++) {
	for (iptr=A->pattern->ptr[ir]-index_offset;iptr<A->pattern->ptr[ir+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
	    if (mask_row[irow]>0. ) {
               #pragma ivdep
	       for (icb=0;icb< A->col_block_size;icb++) {
	           icol=icb+A->col_block_size*(A->pattern->index[iptr]-index_offset);
                   l=iptr*A->block_size+irb+A->row_block_size*icb;
	           if (irow==icol) {
		      A->val[l]=main_diagonal_value;
		    } else {
		      A->val[l]=0;
		    }
	          }
	       }
	    }
	}
  }
}
