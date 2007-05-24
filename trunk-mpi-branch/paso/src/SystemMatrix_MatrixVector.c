/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: matrix vector product with sparse matrix           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/**************************************************************************/

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */

/* has to be called within a parallel region                              */
/* barrier synconization is performed to make sure that the input vector available */

void  Paso_SystemMatrix_MatrixVector(double alpha,
    Paso_SystemMatrix* A,
    double* in,
    double beta,
    double* out) {

  double *buffer0=NULL, *buffer1=NULL;
  dim_t N;
  Paso_MPIInfo *mpi_info=A->mpi_info;

  if (A->type & MATRIX_FORMAT_CSC) {
     if ( mpi_info->size>1) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSC is not supported by MPI.");
           return;
     } else {
       if (A->type & MATRIX_FORMAT_OFFSET1) {
         Paso_SystemMatrix_MatrixVector_CSC_OFFSET1(alpha,A,in,beta,out);
       } else {
         Paso_SystemMatrix_MatrixVector_CSC_OFFSET0(alpha,A,in,beta,out);
       }
     }
  } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: TRILINOS is not supported with MPI.");
           return;
  } else {
     if (A->type & MATRIX_FORMAT_OFFSET1) {
           if ( mpi_info->size>1) {
              Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSR with index 1 is not supported by MPI.");
              return;
           } else {
              Paso_SystemMatrix_MatrixVector_CSR_OFFSET1(alpha,A,in,beta,out);
           }
     } else {
          if ( mpi_info->size>1) {
               N=A->maxNumCols * A->col_block_size;
               buffer0=TMPMEMALLOC(N,double);
               buffer1=TMPMEMALLOC(N,double);
               if (Finley_checkPtr(buffer0) ||  Finley_checkPtr(buffer1) ) {
                   TMPMEMFREE(buffer0);
                   TMPMEMFREE(buffer1);
                   return;
               }
          } else{
               buffer0=NULL;
               buffer1=NULL;
          }
          Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(alpha,A,in,beta,out,buffer0,buffer1);
          if ( mpi_info->size>1) {
               TMPMEMFREE(buffer0);
               TMPMEMFREE(buffer1);
          }
     }
  }
  return;
}

void  Paso_SystemMatrix_MatrixVector_CSC_OFFSET0(double alpha,
    Paso_SystemMatrix* A,
    double* in,
    double beta,
    double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;
  register double reg,reg1,reg2,reg3;
  #pragma omp barrier

  if (ABS(beta)>0.) {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] *= beta;
  } else {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
      
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        /* TODO: parallelize (good luck!) */
        #pragma omp single
	for (icol=0;icol< A->pattern->myNumOutput;++icol) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[icol];iptr<A->pattern->ptr[icol+1]; ++iptr) {
	    out[A->pattern->index[iptr]]+= alpha * A->val[iptr] * in[icol];
	  }
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic];iptr<A->pattern->ptr[ic+1]; iptr++) {
	       ic=2*(A->pattern->index[iptr]);
	       out[  2*ir] += alpha * ( A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic] );
	       out[1+2*ir] += alpha * ( A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic] );
	  }
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic];iptr<A->pattern->ptr[ic+1]; iptr++) {
	      ir=3*(A->pattern->index[iptr]);
              out[  3*ir] += alpha * ( A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic] );
	      out[1+3*ir] += alpha * ( A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic] );
	      out[2+3*ir] += alpha * ( A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic] );
	  }
	}
    } else {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
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

void  Paso_SystemMatrix_MatrixVector_CSC_OFFSET1(double alpha,
    Paso_SystemMatrix* A,
    double* in,
    double beta,
    double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;
  register double reg,reg1,reg2,reg3;
  #pragma omp barrier

  if (ABS(beta)>0.) {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] *= beta;
  } else {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
      
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        /* TODO: parallelize (good luck!) */
        #pragma omp single
	for (icol=0;icol< A->pattern->myNumOutput;++icol) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[icol]-1;iptr<A->pattern->ptr[icol+1]-1; ++iptr) {
	    out[A->pattern->index[iptr]-1]+= alpha * A->val[iptr] * in[icol];
	  }
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
	  for (iptr=A->pattern->ptr[ic]-1;iptr<A->pattern->ptr[ic+1]-1; iptr++) {
	       ic=2*(A->pattern->index[iptr]-1);
	       out[  2*ir] += alpha * ( A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic] );
	       out[1+2*ir] += alpha * ( A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic] );
	  }
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ic]-1;iptr<A->pattern->ptr[ic+1]-1; iptr++) {
	      ir=3*(A->pattern->index[iptr]-1);
              out[  3*ir] += alpha * ( A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic] );
	      out[1+3*ir] += alpha * ( A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic] );
	      out[2+3*ir] += alpha * ( A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic] );
	  }
	}
    } else {
        /* TODO: parallelize */
        #pragma omp single
	for (ic=0;ic< A->pattern->myNumOutput;ic++) {
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

void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET1(double alpha,
    Paso_SystemMatrix* A,
    double* in,
    double beta,
    double* out) {

  register index_t ir,icol,iptr,icb,irb,irow,ic;
  register double reg,reg1,reg2,reg3;
  #pragma omp barrier

  if (ABS(beta)>0.) {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] *= beta;
  } else {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;irow++) 
      out[irow] = 0;
  }
  /*  do the operation: */
  if (ABS(alpha)>0) {
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        #pragma omp for private(irow,iptr,reg) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;++irow) {
          reg=0.;
          #pragma ivdep
	  for (iptr=(A->pattern->ptr[irow])-1;iptr<(A->pattern->ptr[irow+1])-1; ++iptr) {
	      reg += A->val[iptr] * in[A->pattern->index[iptr]-1];
	  }
	  out[irow] += alpha * reg;
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        #pragma omp for private(ir,reg1,reg2,iptr,ic) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
          reg1=0.;
          reg2=0.;
	  for (iptr=A->pattern->ptr[ir]-1;iptr<A->pattern->ptr[ir+1]-1; iptr++) {
	       ic=2*(A->pattern->index[iptr]-1);
	       reg1 += A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic];
	       reg2 += A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic];
	  }
	  out[  2*ir] += alpha * reg1;
	  out[1+2*ir] += alpha * reg2;
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        #pragma omp for private(ir,reg1,reg2,reg3,iptr,ic) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
          reg1=0.;
          reg2=0.;
          reg3=0.;
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
        #pragma omp for private(ir,iptr,irb,icb,irow,icol,reg) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
	  for (iptr=A->pattern->ptr[ir]-1;iptr<A->pattern->ptr[ir+1]-1; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
              reg=0.;
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

void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha,
    Paso_SystemMatrix* A,
    double* in,
    double beta,
    double* out,
    double* buffer0,
    double* buffer1) {

  register index_t irow;
  Paso_MPIInfo *mpi_info=A->mpi_info;
  #pragma omp barrier
  if (ABS(beta)>0.) {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;++irow)
      out[irow] *= beta;
  } else {
    #pragma omp for private(irow) schedule(static)
    for (irow=0;irow < A->myNumRows * A->row_block_size;++irow) 
      out[irow] = 0;
  }
  /*  do the operation: */
  if (ABS(alpha)>0) {
      if ( mpi_info->size>1) {
         #ifdef PASO_MPI
         {
             MPI_Status status;
             MPI_Request snd_rqst, rcv_rqst;
             double *snd_buf=in, *rcv_buf=buffer0, *swap;
             dim_t len_snd_buf= A->myNumCols * A->col_block_size;
             dim_t len_rcv_buf=A->maxNumCols * A->col_block_size;
             dim_t h, fromRank, toRank;
             dim_t myRank=mpi_info->rank;
             dim_t numProc=mpi_info->size;
             index_t rank_of_snd_buf=mpi_info->rank;

             for (h=0; h<A->pattern->numHops; ++h) {

                 /* start asynchronos communication */
                 if (h<A->pattern->numHops-1) {
                   fromRank=PASO_MPI_mod(myRank-(A->pattern->hop[h]),numProc);
                   toRank=PASO_MPI_mod(myRank+(A->pattern->hop[h]),numProc);
                   mpi_info->msg_tag_counter++;

                   #pragma omp master
                   {
                      MPI_Issend(snd_buf,len_snd_buf,MPI_DOUBLE,toRank,
                                 mpi_info->msg_tag_counter,mpi_info->comm,&rcv_rqst);
                      MPI_Irecv(rcv_buf,len_rcv_buf,MPI_DOUBLE,fromRank, 
                                mpi_info->msg_tag_counter,mpi_info->comm,&snd_rqst);
                   }
                 }
                 /* MVM with send buffer as input */
                 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_P(alpha,
                                                    A,
                                                    snd_buf,
                                                    A->pattern->input_distribution->first_component[rank_of_snd_buf],
                                                    A->pattern->input_distribution->first_component[rank_of_snd_buf+1],
                                                    out);
                 /* wait communication to be finished */
                 if (h<A->pattern->numHops-1) {
                    #pragma omp master
                    {
                       MPI_Wait(&rcv_rqst,&status);
                       MPI_Wait(&snd_rqst,&status);
                    }
                 }
                 /*  swap recieve and send buffer  */
                 if (h==0) {
                     snd_buf = rcv_buf;
                     len_snd_buf = len_rcv_buf;
                     rcv_buf = buffer1;
                 } else {
                     swap = snd_buf;
                     snd_buf = rcv_buf;
                     rcv_buf = swap;
                 }
                 rank_of_snd_buf=PASO_MPI_mod(rank_of_snd_buf-(A->pattern->hop[h]),numProc);
                 
           }
        }
        #endif
      } else {
           Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_S(alpha,A,in,out);
      }

  }
}

void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_S(
           double alpha,
           Paso_SystemMatrix* A,
           double* in,
           double* out) 
{
    register index_t ir,icol,iptr,icb,irb,irow,ic,Aiptr;
    register double reg,reg1,reg2,reg3,in1,in2,in3,A00,A10,A20,A01,A11,A21,A02,A12,A22;
    if (A ->col_block_size==1 && A->row_block_size ==1) {
        #pragma omp for private(irow,iptr,reg) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;++irow) {
          reg=0.;
	  for (iptr=(A->pattern->ptr[irow]);iptr<(A->pattern->ptr[irow+1]); ++iptr) {
	      reg += A->val[iptr] * in[A->pattern->index[iptr]];
	  }
	  out[irow] += alpha * reg;
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
        #pragma omp for private(ir,reg1,reg2,iptr,ic,Aiptr,in1,in2,A00,A10,A01,A11) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
          reg1=0.;
          reg2=0.;
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	       ic=2*(A->pattern->index[iptr]);
               Aiptr=iptr*4;
               in1=in[ic];
               in2=in[1+ic];
               A00=A->val[Aiptr  ];
               A10=A->val[Aiptr+1];
               A01=A->val[Aiptr+2];
               A11=A->val[Aiptr+3];
	       reg1 += A00*in1 + A01*in2;
	       reg2 += A10*in1 + A11*in2;
	  }
	  out[  2*ir] += alpha * reg1;
	  out[1+2*ir] += alpha * reg2;
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
        #pragma omp for private(ir,reg1,reg2,reg3,iptr,ic,Aiptr,in1,in2,in3,A00,A10,A20,A01,A11,A21,A02,A12,A22) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
          reg1=0.;
          reg2=0.;
          reg3=0.;
          #pragma ivdep
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	       ic=3*(A->pattern->index[iptr]);
               Aiptr=iptr*9;
               in1=in[ic];
               in2=in[1+ic];
               in3=in[2+ic];
               A00=A->val[Aiptr  ];
               A10=A->val[Aiptr+1];
               A20=A->val[Aiptr+2];
               A01=A->val[Aiptr+3];
               A11=A->val[Aiptr+4];
               A21=A->val[Aiptr+5];
               A02=A->val[Aiptr+6];
               A12=A->val[Aiptr+7];
               A22=A->val[Aiptr+8];
	       reg1 += A00*in1 + A01*in2 + A02*in3;
	       reg2 += A10*in1 + A11*in2 + A12*in3;
	       reg3 += A20*in1 + A21*in2 + A22*in3;
	  }
	  out[  3*ir] += alpha * reg1;
	  out[1+3*ir] += alpha * reg2;
	  out[2+3*ir] += alpha * reg3;
	}
    } else {
        #pragma omp for private(ir,iptr,irb,icb,irow,icol,reg) schedule(static)
	for (ir=0;ir< A->pattern->myNumOutput;ir++) {
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
              reg=0.;
              #pragma ivdep
	      for (icb=0;icb< A->col_block_size;icb++) {
		icol=icb+A->col_block_size*(A->pattern->index[iptr]);
		reg += A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
	      }
	      out[irow] += alpha * reg;
	    }
	  }
	}
    }
    return;
}
void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_P(
           double alpha,
           Paso_SystemMatrix* A,
           double* in,
           index_t min_index, index_t max_index,
           double* out) 
{
    register index_t iptr,icb,irb,irow,ic,Aiptr,itmp,irtmp,ictmp;
    register double reg,reg1,reg2,reg3,in1,in2,in3,A00,A10,A20,A01,A11,A21,A02,A12,A22;
    if (A ->col_block_size==1 && A->row_block_size ==1) {
	int iptr0;
        #pragma omp for private(irow,iptr,reg,itmp,iptr0) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;++irow) {
          reg=0.;
          #pragma ivdep
	  /* index[] is an increasing sequence, can find min first */
	  for (iptr0=(A->pattern->ptr[irow]);iptr0<(A->pattern->ptr[irow+1]); ++iptr0) {
              itmp=A->pattern->index[iptr0];
              if ( min_index<= itmp) { break; }
	  }
          #pragma ivdep
	  for (iptr=iptr0;iptr<(A->pattern->ptr[irow+1]); ++iptr) {
              itmp=A->pattern->index[iptr];
	      /* index[] is an increasing sequence, can quit loop if at max */
	      if (itmp >= max_index) { break; }
	      reg += A->val[iptr] * in[itmp-min_index];
	  }
	  out[irow] += alpha * reg;
	}
    } else if (A ->col_block_size==2 && A->row_block_size ==2) {
	int iptr0;
        #pragma omp for private(irow,reg1,reg2,iptr,ic,Aiptr,in1,in2,A00,A10,A01,A11,itmp,iptr0) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;++irow) {
          reg1=0.;
          reg2=0.;
          #pragma ivdep
	  /* index[] is an increasing sequence, can find min first */
	  for (iptr0=(A->pattern->ptr[irow]);iptr0<(A->pattern->ptr[irow+1]); ++iptr0) {
              itmp=A->pattern->index[iptr0];
              if ( min_index<= itmp) { break; }
	  }
          #pragma ivdep
	  for (iptr=iptr0;iptr<(A->pattern->ptr[irow+1]); ++iptr) {
              itmp=A->pattern->index[iptr];
	      /* index[] is an increasing sequence, can quit loop if at max */
	      if (itmp >= max_index) { break; }
	      ic=2*(itmp-min_index);
              Aiptr=iptr*4;
              in1=in[ic];
              in2=in[1+ic];
              A00=A->val[Aiptr  ];
              A10=A->val[Aiptr+1];
              A01=A->val[Aiptr+2];
              A11=A->val[Aiptr+3];
	      reg1 += A00*in1 + A01*in2;
	      reg2 += A10*in1 + A11*in2;
	  }
	  out[  2*irow] += alpha * reg1;
	  out[1+2*irow] += alpha * reg2;
	}
    } else if (A ->col_block_size==3 && A->row_block_size ==3) {
	int iptr0;
        #pragma omp for private(irow,reg1,reg2,reg3,iptr,ic,Aiptr,in1,in2,in3,A00,A10,A20,A01,A11,A21,A02,A12,A22,itmp,iptr0) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;irow++) {
          reg1=0.;
          reg2=0.;
          reg3=0.;
          #pragma ivdep
	  /* index[] is an increasing sequence, can find min first */
	  for (iptr0=(A->pattern->ptr[irow]);iptr0<(A->pattern->ptr[irow+1]); ++iptr0) {
              itmp=A->pattern->index[iptr0];
              if ( min_index<= itmp) { break; }
	  }
          #pragma ivdep
	  for (iptr=iptr0;iptr<(A->pattern->ptr[irow+1]); ++iptr) {
              itmp=A->pattern->index[iptr];
	      /* index[] is an increasing sequence, can quit loop if at max */
	      if (itmp >= max_index) { break; }
	      ic=3*(itmp-min_index);
              Aiptr=iptr*9;
              in1=in[ic];
              in2=in[1+ic];
              in3=in[2+ic];
              A00=A->val[Aiptr  ];
              A10=A->val[Aiptr+1];
              A20=A->val[Aiptr+2];
              A01=A->val[Aiptr+3];
              A11=A->val[Aiptr+4];
              A21=A->val[Aiptr+5];
              A02=A->val[Aiptr+6];
              A12=A->val[Aiptr+7];
              A22=A->val[Aiptr+8];
	      reg1 += A00*in1 + A01*in2 + A02*in3;
	      reg2 += A10*in1 + A11*in2 + A12*in3;
	      reg3 += A20*in1 + A21*in2 + A22*in3;
	  }
	  out[  3*irow] += alpha * reg1;
	  out[1+3*irow] += alpha * reg2;
	  out[2+3*irow] += alpha * reg3;
	}
    } else {
	int iptr0;
        #pragma omp for private(irow,iptr,irb,icb,irtmp,ictmp,reg,itmp,iptr0) schedule(static)
	for (irow=0;irow< A->pattern->myNumOutput;irow++) {
          #pragma ivdep
	  /* index[] is an increasing sequence, can find min first */
	  for (iptr0=(A->pattern->ptr[irow]);iptr0<(A->pattern->ptr[irow+1]); ++iptr0) {
              itmp=A->pattern->index[iptr0];
              if ( min_index<= itmp) { break; }
	  }
          #pragma ivdep
	  for (iptr=iptr0;iptr<(A->pattern->ptr[irow+1]); ++iptr) {
              itmp=A->pattern->index[iptr];
	      /* index[] is an increasing sequence, can quit loop if at max */
	      if (itmp >= max_index) { break; }
	      for (irb=0;irb< A->row_block_size;irb++) {
	        irtmp=irb+A->row_block_size*irow;
                reg=0.;
	        for (icb=0;icb< A->col_block_size;icb++) {
		  ictmp=icb+A->col_block_size*(itmp-min_index);
		  reg += A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[ictmp];
	        }
	        out[irtmp] += alpha * reg;
	      }
          }
        }
    }
    return;
}
