/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix                                       */

/*  nullify rows and columns in the matrix                    */

/*  the rows and columns are marked by positive values in     */
/*  mask_row and mask_col. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row and      */
/*  mask_col are set to main_diagonal_value                   */


/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "System.h"

/**************************************************************/

void  Finley_SystemMatrix_nullifyRowsAndCols(Finley_SystemMatrix* A, 
                         escriptDataC* mask_row, escriptDataC* mask_col, double main_diagonal_value) {

  /* check structure: */

  if (!isExpanded(mask_row)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"row mask Data object has to be expanded");
  } else if (!isExpanded(mask_col)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"column mask Data object has to be expanded");
  } else if (getLength(mask_col) != A ->col_block_size *A ->num_cols ) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"column mask vector and matrix column dimensions don't match.");
  } else if (getLength(mask_row) != A ->row_block_size * A ->num_rows) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"row mask vector and matrix row dimensions don't match.");
  }

  /*  do the operation: */
  if (Finley_ErrorCode==NO_ERROR) {
     Finley_SystemMatrixNullify(A,getSampleData(mask_row,0),getSampleData(mask_col,0),main_diagonal_value);
  }
 
  return;
}

/****************************************************************/

/* this is the actual ellimination */

void Finley_SystemMatrixNullify(Finley_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {

  maybelong ir,icol,iptr,icb,irb,irow,ic;

  if (A ->col_block_size==1 && A ->row_block_size ==1) {
    switch(A->type) {
    case CSR:
      #pragma omp parallel for private(irow, iptr,icol) schedule(static)
      for (irow=0;irow< A->num_rows;irow++) {
	/* TODO: test mask_row here amd not inside every loop */
	for (iptr=A->ptr[irow]-PTR_OFFSET;iptr<A->ptr[irow+1]-PTR_OFFSET; iptr++) {
	  icol=A->index[iptr]-INDEX_OFFSET;
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
      for (icol=0;icol< A->num_cols;icol++) {
	for (iptr=A->ptr[icol]-PTR_OFFSET;iptr<A->ptr[icol+1]-PTR_OFFSET; iptr++) {
	  irow=A->index[iptr]-INDEX_OFFSET;
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
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Unknown matrix type.");
    } /* switch A->type */
  } else {
    switch(A->type) {
    case CSR:
      #pragma omp parallel for private(irow, iptr,icol,ir,irb,icb) schedule(static)
      for (ir=0;ir< A->num_rows;ir++) {
	for (iptr=A->ptr[ir]-PTR_OFFSET;iptr<A->ptr[ir+1]-PTR_OFFSET; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*(A->index[iptr]-INDEX_OFFSET);
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
		if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
		  A->val[iptr]=main_diagonal_value;
		} else {
		  A->val[iptr]=0;
		}
	      }
	    }
	  }
	}
      }
      break;
    case CSC:
      #pragma omp parallel for private(irow, iptr,icol,ic,irb,icb) schedule(static)
      for (ic=0;ic< A->num_cols;ic++) {
	for (iptr=A->ptr[ic]-PTR_OFFSET;iptr<A->ptr[ic+1]-PTR_OFFSET; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*(A->index[iptr]-INDEX_OFFSET);
	    for (icb=0;icb< A->col_block_size;icb++) {
	      icol=icb+A->col_block_size*ic;
	      if (mask_col[icol]>0. || mask_row[irow]>0. ) {
		if (irow==icol && mask_col[icol]>0. && mask_row[irow]>0.) {
		  A->val[iptr]=main_diagonal_value;
		} else {
		  A->val[iptr]=0;
		}
	      }
	    }
	  }
	}
      }
      break;
    default:
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Unknown matrix type.");
    } /* switch A->type */
  }
  return;
}

/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:46  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
