/* $Id$ */

/**************************************************************/

/* Finley: matrix vector product with sparse matrix           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"

/**************************************************************/

/*  vector update operation: out=out+A*in */

void  Finley_SystemMatrixVector(escriptDataC * out,Finley_SystemMatrix* A, escriptDataC * in) {
  Finley_ScaledSystemMatrixVector(1.0, A, in, 1.0, out);
}

/**************************************************************/

/*  scaled vector update operation: out = alpha * A * in + beta * out */

void  Finley_ScaledSystemMatrixVector(double alpha,
				Finley_SystemMatrix* A,
				escriptDataC * in,
				double beta,
				escriptDataC * out) {
  Finley_ErrorCode=NO_ERROR;
  /* check structure: */

  if (!isExpanded(in)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"matrix vector product : input Data object has to be expanded");
  } else if (!isExpanded(out)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"matrix vector product : output Data object has to be expanded");
  } else if ( getLength(in)  != A ->col_block_size * A ->num_cols) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"matrix vector product : lnegth of input Data object and matrix column dimensions don't match.");
  } else if (getLength(out) != A ->row_block_size * A ->num_rows) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"matrix vector product : lnegth of output Data object and matrix row dimensions don't match.");
  } else if (A ->symmetric) {
    Finley_ErrorCode=SYSTEM_ERROR;
    sprintf(Finley_ErrorMsg,"matrix-vector product for symmetrix matric has not been implemented yet.");
  }

  /* delegate to routine */

  if (Finley_ErrorCode==NO_ERROR) {
     #pragma omp parallel 
     Finley_RawScaledSystemMatrixVector(alpha, A, getSampleData(in,0), beta, getSampleData(out,0));
  }
  return;
}

/**************************************************************************/

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */

/* has to be called within a parallel region                              */
/* barrier synconization is performed to make sure that the input vector available */

void  Finley_RawScaledSystemMatrixVector(double alpha,
    Finley_SystemMatrix* A,
    double* in,
    double beta,
    double* out) {

  maybelong ir,icol,iptr,icb,irb,irow,ic;
  #pragma omp barrier

  if (ABS(beta)>0.) {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->num_rows * A->row_block_size;irow++) 
      out[irow] *= beta;
  } else {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->num_rows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
      
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
      switch(A->type) {
      case CSR:
        #pragma omp for private(irow,iptr) schedule(static)
	for (irow=0;irow< A->num_rows;irow++) {
	  for (iptr=A->ptr[irow]-PTR_OFFSET;iptr<A->ptr[irow+1]-PTR_OFFSET; iptr++) {
	    out[irow] += alpha * A->val[iptr] * in[A->index[iptr]-INDEX_OFFSET];
	  }
	}
	break;
      case CSC:
        /* TODO: parallelize */
        #pragma omp single
	for (icol=0;icol< A->num_cols;icol++) {
	  for (iptr=A->ptr[icol]-PTR_OFFSET;iptr<A->ptr[icol+1]-PTR_OFFSET; iptr++) {
	    out[A->index[iptr]-INDEX_OFFSET]+= alpha * A->val[iptr] * in[icol];
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
        #pragma omp for private(ir,iptr,irb,icb,irow,icol) schedule(static)
	for (ir=0;ir< A->num_rows;ir++) {
	  for (iptr=A->ptr[ir]-PTR_OFFSET;iptr<A->ptr[ir+1]-PTR_OFFSET; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*(A->index[iptr]-INDEX_OFFSET);
		out[irow] += alpha * A->val[iptr*A->row_block_size*A->col_block_size+icb+A->col_block_size*irb] * in[icol];
	      }
	    }
	  }
	}
	break;
      case CSC:
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->num_cols;ic++) {
	  for (iptr=A->ptr[ic]-PTR_OFFSET;iptr<A->ptr[ic+1]-PTR_OFFSET; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*(A->index[iptr]-INDEX_OFFSET);
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*ic;
		out[irow] += alpha * A->val[iptr*A->row_block_size*A->col_block_size+icb+A->col_block_size*irb] * in[icol];
	      }
	    }
	  }
	}
      default:
	Finley_ErrorCode=TYPE_ERROR;
	sprintf(Finley_ErrorMsg,"Unknown matrix type.");
      } /* switch A->type */
    }
  }
  return;
}

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2.2.2  2004/10/26 06:36:39  jgs
 * committing Lutz's changes to branch jgs
 *
 * Revision 1.3  2004/10/21 04:55:54  gross
 * bug in CSC MVM fixed
 *
 * Revision 1.2  2004/08/13 00:12:53  gross
 * Gradtest is getting further now. PDE assemblage has been added but not tested.
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
