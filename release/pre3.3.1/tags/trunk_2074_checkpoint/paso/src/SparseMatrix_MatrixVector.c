
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: raw scaled vector update operation: out = alpha * A * in + beta * out */

/**************************************************************/

/* Author: gross@access.edu.au */

/**************************************************************/

#include "SparseMatrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* CSC format with offset 0*/
void  Paso_SparseMatrix_MatrixVector_CSC_OFFSET0(const double alpha,
                                                 const Paso_SparseMatrix* A,
                                                 const double* in,
                                                 const double beta,
                                                 double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;

  if (ABS(beta)>0.) {
    if (beta != 1.) {
        #pragma omp parallel for private(irow) schedule(static)
        for (irow=0;irow < A->numRows * A->row_block_size;irow++) 
          out[irow] *= beta;
    }
  } else {
    #pragma omp parallel for private(irow) schedule(static)
    for (irow=0;irow < A->numRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
  if (Paso_Pattern_isEmpty(A->pattern)) return;
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        /* TODO: parallelize (good luck!) */
	for (icol=0;icol< A->pattern->numOutput;++icol) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[icol];iptr<A->pattern->ptr[icol+1]; ++iptr) {
	    out[A->pattern->index[iptr]]+= alpha * A->val[iptr] * in[icol];
	  }
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic];iptr<A->pattern->ptr[ic+1]; iptr++) {
	       ir=2*(A->pattern->index[iptr]);
	       out[  ir] += alpha * ( A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic] );
	       out[1+ir] += alpha * ( A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic] );
	  }
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic];iptr<A->pattern->ptr[ic+1]; iptr++) {
	      ir=3*(A->pattern->index[iptr]);
              out[  ir] += alpha * ( A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic] );
	      out[1+ir] += alpha * ( A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic] );
	      out[2+ir] += alpha * ( A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic] );
	  }
	}
    } else {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
	  for (iptr=A->pattern->ptr[ic];iptr<A->pattern->ptr[ic+1]; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*(A->pattern->index[iptr]);
              #pragma ivdep
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*ic;
		out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	    }
	  }
	}
    }
  }
  return;
}

/* CSC format with offset 1*/
void  Paso_SparseMatrix_MatrixVector_CSC_OFFSET1(const double alpha,
                                                 const Paso_SparseMatrix* A,
                                                 const double* in,
                                                 const double beta,
                                                 double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;

  if (ABS(beta)>0.) {
    if (beta != 1.) {
        #pragma omp parallel for private(irow) schedule(static)
        for (irow=0;irow < A->numRows * A->row_block_size;irow++) {
          out[irow] *= beta;
        }
    }
  } else {
    #pragma omp parallel for private(irow) schedule(static)
    for (irow=0;irow < A->numRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
      
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        /* TODO: parallelize (good luck!) */
	for (icol=0;icol< A->pattern->numOutput;++icol) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[icol]-1;iptr<A->pattern->ptr[icol+1]-1; ++iptr) {
	    out[A->pattern->index[iptr]-1]+= alpha * A->val[iptr] * in[icol];
	  }
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
	  for (iptr=A->pattern->ptr[ic]-1;iptr<A->pattern->ptr[ic+1]-1; iptr++) {
	       ir=2*(A->pattern->index[iptr]-1);
	       out[  ir] += alpha * ( A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic] );
	       out[1+ir] += alpha * ( A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic] );
	  }
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic]-1;iptr<A->pattern->ptr[ic+1]-1; iptr++) {
	      ir=3*(A->pattern->index[iptr]-1);
              out[  ir] += alpha * ( A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic] );
	      out[1+ir] += alpha * ( A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic] );
	      out[2+ir] += alpha * ( A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic] );
	  }
	}
    } else {
        /* TODO: parallelize */
	for (ic=0;ic< A->pattern->numOutput;ic++) {
	  for (iptr=A->pattern->ptr[ic]-1;iptr<A->pattern->ptr[ic+1]-1; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*(A->pattern->index[iptr]-1);
              #pragma ivdep
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*ic;
		out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	    }
	  }
	}
    }
  }
  return;
}
/* CSR format with offset 1*/
void  Paso_SparseMatrix_MatrixVector_CSR_OFFSET1(const double alpha,
    const Paso_SparseMatrix* A,
    const double* in,
    const double beta,
    double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;
  register double reg,reg1,reg2,reg3;
  if (ABS(beta)>0.) {
    if (beta != 1.) {
        #pragma omp parallel for private(irow) schedule(static)
        for (irow=0;irow < A->numRows * A->row_block_size;irow++) 
          out[irow] *= beta;
    }
  } else {
    #pragma omp parallel for private(irow) schedule(static)
    for (irow=0;irow < A->numRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        #pragma omp parallel for private(irow,iptr,reg) schedule(static)
	for (irow=0;irow< A->pattern->numOutput;++irow) {
          reg=0.;
          #pragma ivdep
	  for (iptr=(A->pattern->ptr[irow])-1;iptr<(A->pattern->ptr[irow+1])-1; ++iptr) {
	      reg += A->val[iptr] * in[A->pattern->index[iptr]-1];
	  }
	  out[irow] += alpha * reg;
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        #pragma omp parallel for private(ir,reg1,reg2,iptr,ic) schedule(static)
	for (ir=0;ir< A->pattern->numOutput;ir++) {
          reg1=0.;
          reg2=0.;
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ir]-1;iptr<A->pattern->ptr[ir+1]-1; iptr++) {
	       ic=2*(A->pattern->index[iptr]-1);
	       reg1 += A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic];
	       reg2 += A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic];
	  }
	  out[  2*ir] += alpha * reg1;
	  out[1+2*ir] += alpha * reg2;
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        #pragma omp parallel for private(ir,reg1,reg2,reg3,iptr,ic) schedule(static)
	for (ir=0;ir< A->pattern->numOutput;ir++) {
          reg1=0.;
          reg2=0.;
          reg3=0.;
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ir]-1;iptr<A->pattern->ptr[ir+1]-1; iptr++) {
	       ic=3*(A->pattern->index[iptr]-1);
	       reg1 += A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic];
	       reg2 += A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic];
	       reg3 += A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic];
	  }
	  out[  3*ir] += alpha * reg1;
	  out[1+3*ir] += alpha * reg2;
	  out[2+3*ir] += alpha * reg3;
	}
    } else {
        #pragma omp parallel for private(ir,iptr,irb,icb,irow,icol,reg) schedule(static)
	for (ir=0;ir< A->pattern->numOutput;ir++) {
	  for (iptr=A->pattern->ptr[ir]-1;iptr<A->pattern->ptr[ir+1]-1; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
              reg=0.;
              #pragma ivdep
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*(A->pattern->index[iptr]-1);
		reg += A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	      out[irow] += alpha * reg;
	    }
	  }
	}
    }
  }
  return;
}
/* CSR format with offset 0*/
void  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(const double alpha,
                                                 const Paso_SparseMatrix* A,
                                                 const double* in,
                                                 const double beta,
                                                 double* out) 
{
/*#define PASO_DYNAMIC_SCHEDULING_MVM */

#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined _OPENMP 
#define USE_DYNAMIC_SCHEDULING
#endif

#ifdef USE_DYNAMIC_SCHEDULING
    char* chksz_chr=NULL;
    dim_t chunk_size=1;
#endif

    dim_t nrow=A->numRows;
    dim_t np=1, len, rest, irow, local_n, p;

#ifdef USE_DYNAMIC_SCHEDULING
    dim_t n_chunks;
#endif

#ifdef _OPENMP
    np=omp_get_max_threads();
#endif
#ifdef USE_DYNAMIC_SCHEDULING
    chksz_chr=getenv("PASO_CHUNK_SIZE_MVM");
    if (chksz_chr!=NULL) sscanf(chksz_chr, "%d",&chunk_size);
    chunk_size=MIN(MAX(1,chunk_size),nrow/np);
    n_chunks=nrow/chunk_size;
    if (n_chunks*chunk_size<nrow) n_chunks+=1;
#else
    len=nrow/np;
    rest=nrow-len*np;
#endif
    
     #pragma omp parallel private(irow, p, local_n)
     {
        #ifdef USE_DYNAMIC_SCHEDULING
          #pragma omp for schedule(dynamic,1)
          for (p=0; p<n_chunks;p++) {
            irow=chunk_size*p;
            local_n=MIN(chunk_size,nrow-chunk_size*p);
        #else
            #pragma omp for schedule(static)
            for (p=0; p<np;p++) {
               irow=len*p+MIN(p,rest);
               local_n=len+(p<rest ? 1 :0 );
        #endif
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(alpha,
                                                          local_n,
                                                          A->row_block_size,
                                                          A->col_block_size,
                                                          &(A->pattern->ptr[irow]),
                                                           A->pattern->index, A->val, in, beta, &out[irow*A->row_block_size]);
        #ifdef USE_DYNAMIC_SCHEDULING
          }
        #else
          }
        #endif
   }
}
/* CSR format with offset 0*/
void  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(const double alpha,
                                                        const dim_t nRows,
                                                        const dim_t row_block_size,
                                                        const dim_t col_block_size,
                                                        const index_t* ptr,
                                                        const index_t* index,
                                                        const double* val,
                                                        const double* in,
                                                        const double beta,
                                                        double* out) 
{
    register index_t ir,icol,iptr,icb,irb,irow,ic,Aiptr;
    register double reg,reg1,reg2,reg3,in1,in2,in3,A00,A10,A20,A01,A11,A21,A02,A12,A22;
    dim_t block_size;
    if (ABS(beta)>0.) {
      if (beta != 1.) {
          for (irow=0;irow < nRows * row_block_size;irow++) 
            out[irow] *= beta;
      }
    } else {
      for (irow=0;irow < nRows * row_block_size;irow++) 
        out[irow] = 0;
    }
    if (ABS(alpha)>0) {
      if (col_block_size==1 && row_block_size ==1) {
  	for (irow=0;irow< nRows;++irow) {
          reg=0.;
          #pragma ivdep
  	  for (iptr=ptr[irow];iptr<ptr[irow+1]; ++iptr) {
  	      reg += val[iptr] * in[index[iptr]];
  	  }
  	  out[irow] += alpha * reg;
  	}
      } else if (col_block_size==2 && row_block_size ==2) {
  	for (ir=0;ir< nRows;ir++) {
          reg1=0.;
          reg2=0.;
          #pragma ivdep
  	  for (iptr=ptr[ir];iptr<ptr[ir+1]; iptr++) {
  	       ic=2*index[iptr];
                 Aiptr=iptr*4;
                 in1=in[ic];
                 in2=in[1+ic];
                 A00=val[Aiptr  ];
                 A10=val[Aiptr+1];
                 A01=val[Aiptr+2];
                 A11=val[Aiptr+3];
  	       reg1 += A00*in1 + A01*in2;
  	       reg2 += A10*in1 + A11*in2;
  	  }
  	  out[  2*ir] += alpha * reg1;
  	  out[1+2*ir] += alpha * reg2;
  	}
      } else if (col_block_size==3 && row_block_size ==3) {
  	for (ir=0;ir< nRows;ir++) {
          reg1=0.;
          reg2=0.;
          reg3=0.;
          #pragma ivdep
  	  for (iptr=ptr[ir];iptr<ptr[ir+1]; iptr++) {
  	       ic=3*index[iptr];
                 Aiptr=iptr*9;
                 in1=in[ic];
                 in2=in[1+ic];
                 in3=in[2+ic];
                 A00=val[Aiptr  ];
                 A10=val[Aiptr+1];
                 A20=val[Aiptr+2];
                 A01=val[Aiptr+3];
                 A11=val[Aiptr+4];
                 A21=val[Aiptr+5];
                 A02=val[Aiptr+6];
                 A12=val[Aiptr+7];
                 A22=val[Aiptr+8];
  	       reg1 += A00*in1 + A01*in2 + A02*in3;
  	       reg2 += A10*in1 + A11*in2 + A12*in3;
  	       reg3 += A20*in1 + A21*in2 + A22*in3;
  	  }
  	  out[  3*ir] += alpha * reg1;
  	  out[1+3*ir] += alpha * reg2;
  	  out[2+3*ir] += alpha * reg3;
  	}
      } else {
        block_size=col_block_size*row_block_size;
  	for (ir=0;ir< nRows;ir++) {
  	  for (iptr=ptr[ir];iptr<ptr[ir+1]; iptr++) {
  	    for (irb=0;irb< row_block_size;irb++) {
  	      irow=irb+row_block_size*ir;
              reg=0.;
              #pragma ivdep
  	      for (icb=0;icb< col_block_size;icb++) {
  		icol=icb+col_block_size*index[iptr];
  		reg += val[iptr*block_size+irb+row_block_size*icb] * in[icol];
  	      }
  	      out[irow] += alpha * reg;
  	    }
  	  }
  	}
      }
    }
    return;
}
