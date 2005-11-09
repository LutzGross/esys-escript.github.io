/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrix                                       */

/*  nullify rows and columns in the matrix                    */

/*  the rows and columns are marked by positive values in     */
/*  mask_row and mask_col. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row and      */
/*  mask_col are set to main_diagonal_value                   */


/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {

  index_t ir,icol,iptr,icb,irb,irow,ic,l;

  if (A ->col_block_size==1 && A ->row_block_size ==1) {
    switch(A->type) {
    case CSR:
      #pragma omp parallel for private(irow, iptr,icol) schedule(static)
      for (irow=0;irow< A->pattern->n_ptr;irow++) {
	/* TODO: test mask_row here amd not inside every loop */
	for (iptr=A->pattern->ptr[irow]-PTR_OFFSET;iptr<A->pattern->ptr[irow+1]-PTR_OFFSET; iptr++) {
	  icol=A->pattern->index[iptr]-INDEX_OFFSET;
	  if (mask_col[icol]>0. || mask_row[irow]>0. ) {
            if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	  }
	}
      }
      break;
    case CSC:
      #pragma omp parallel for private(irow, iptr,icol) schedule(static)
      for (icol=0;icol< A->pattern->n_ptr;icol++) {
	for (iptr=A->pattern->ptr[icol]-PTR_OFFSET;iptr<A->pattern->ptr[icol+1]-PTR_OFFSET; iptr++) {
	  irow=A->pattern->index[iptr]-INDEX_OFFSET;
	  if (mask_col[icol]>0. || mask_row[irow]>0. ) {
            if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	  }
	}
      }
      break;
    default:
      Paso_setError(TYPE_ERROR,"Unknown matrix type in row and column elimination.");
    } /* switch A->type */
  } else {
    switch(A->type) {
    case CSR:
      #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb) schedule(static)
      for (ir=0;ir< A->pattern->n_ptr;ir++) {
	for (iptr=A->pattern->ptr[ir]-PTR_OFFSET;iptr<A->pattern->ptr[ir+1]-PTR_OFFSET; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*(A->pattern->index[iptr]-INDEX_OFFSET);
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
                l=iptr*A->block_size+irb+A->row_block_size*icb;
		if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
		  A->val[l]=main_diagonal_value;
		} else {
		  A->val[l]=0;
		}
	      }
	    }
	  }
	}
      }
      break;
    case CSC:
      #pragma omp parallel for private(l,irow, iptr,icol,ic,irb,icb) schedule(static)
      for (ic=0;ic< A->pattern->n_ptr;ic++) {
	for (iptr=A->pattern->ptr[ic]-PTR_OFFSET;iptr<A->pattern->ptr[ic+1]-PTR_OFFSET; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*(A->pattern->index[iptr]-INDEX_OFFSET);
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*ic;
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
                l=iptr*A->block_size+irb+A->row_block_size*icb;
		if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
		  A->val[l]=main_diagonal_value;
		} else {
		  A->val[l]=0;
		}
	      }
	    }
	  }
	}
      }
      break;
    default:
      Paso_setError(TYPE_ERROR,"Unknown matrix type in row and column elimination.");
    } /* switch A->type */
  }
  return;
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
