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
	for (irow=0;irow< A->pattern->n_ptr;++irow) {
	  for (iptr=(A->pattern->ptr[irow])-PTR_OFFSET;iptr<(A->pattern->ptr[irow+1])-PTR_OFFSET; ++iptr) {
	    out[irow] += alpha * A->val[iptr] * in[A->pattern->index[iptr]-INDEX_OFFSET];
	  }
	}
	break;
      case CSC:
        /* TODO: parallelize */
        #pragma omp single
	for (icol=0;icol< A->pattern->n_ptr;++icol) {
	  for (iptr=A->pattern->ptr[icol]-PTR_OFFSET;iptr<A->pattern->ptr[icol+1]-PTR_OFFSET; ++iptr) {
	    out[A->pattern->index[iptr]-INDEX_OFFSET]+= alpha * A->val[iptr] * in[icol];
	  }
	}
	break;
      default:
	Finley_ErrorCode=TYPE_ERROR;
	sprintf(Finley_ErrorMsg,"Unknown matrix type in MVM.");
      } /* switch A->type */
    } else {
      switch(A->type) {
      case CSR:
        #pragma omp for private(ir,iptr,irb,icb,irow,icol) schedule(static)
	for (ir=0;ir< A->pattern->n_ptr;ir++) {
	  for (iptr=A->pattern->ptr[ir]-PTR_OFFSET;iptr<A->pattern->ptr[ir+1]-PTR_OFFSET; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*(A->pattern->index[iptr]-INDEX_OFFSET);
		out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	    }
	  }
	}
	break;
      case CSC:
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->n_ptr;ic++) {
	  for (iptr=A->pattern->ptr[ic]-PTR_OFFSET;iptr<A->pattern->ptr[ic+1]-PTR_OFFSET; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*(A->pattern->index[iptr]-INDEX_OFFSET);
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*ic;
		out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	    }
	  }
	}
      default:
	Finley_ErrorCode=TYPE_ERROR;
	sprintf(Finley_ErrorMsg,"Unknown matrix type in MVM.");
      } /* switch A->type */
    }
  }
  return;
}

/*
 * $Log$
 * Revision 1.4  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 *
 *
 */
