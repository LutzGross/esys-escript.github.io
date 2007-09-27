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

  dim_t N=A->maxNumCols * A->col_block_size;
  Paso_MPIInfo *mpi_info=A->mpi_info;
  double* buffer0, *buffer1;
  if (mpi_info->size>1) {
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC with block size 1 is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
          Paso_SystemMatrix_nullifyRows_CSR_BLK1(A,mask_row,main_diagonal_value);
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
           Paso_SystemMatrix_nullifyRows_CSR(A,mask_row,main_diagonal_value);
       }
     }
     #ifdef PASO_MPI
     {
         buffer0=TMPMEMALLOC(N,double);
         buffer1=TMPMEMALLOC(N,double);
         if (Finley_checkPtr(buffer0) ||  Finley_checkPtr(buffer1) ) {
             TMPMEMFREE(buffer0);
             TMPMEMFREE(buffer1);
             return;
         }
         MPI_Status status;
         MPI_Request snd_rqst, rcv_rqst;
         double* snd_buf=mask_col, *rcv_buf=buffer0, *swap;
         dim_t len_snd_buf=N;
         dim_t len_rcv_buf=N;
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
                      MPI_Irecv(rcv_buf,len_rcv_buf,MPI_DOUBLE,fromRank,
                                mpi_info->msg_tag_counter,mpi_info->comm,&snd_rqst);
                      MPI_Issend(snd_buf,len_snd_buf,MPI_DOUBLE,toRank,
                                 mpi_info->msg_tag_counter,mpi_info->comm,&rcv_rqst);
                }
             }
             /* annulate colums as input */
             if (A ->col_block_size==1 && A ->row_block_size ==1) {
                 Paso_SystemMatrix_nullifyCols_CSR_BLK1(A,snd_buf,main_diagonal_value,
                                                    A->pattern->input_distribution->first_component[rank_of_snd_buf],
                                                    A->pattern->input_distribution->first_component[rank_of_snd_buf+1]);
             } else {
                 Paso_SystemMatrix_nullifyCols_CSR(A,snd_buf,main_diagonal_value,
                                                   A->pattern->input_distribution->first_component[rank_of_snd_buf],
                                                   A->pattern->input_distribution->first_component[rank_of_snd_buf+1]);
             }
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
             rank_of_snd_buf=fromRank;
         }
         TMPMEMFREE(buffer0);
         TMPMEMFREE(buffer1);
     }
     #endif
  } else { 
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_SystemMatrix_nullifyRowsAndCols_CSC_BLK1(A,mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported.");
       } else {
           Paso_SystemMatrix_nullifyRowsAndCols_CSR_BLK1(A,mask_row,mask_col,main_diagonal_value);
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_SystemMatrix_nullifyRowsAndCols_CSC(A,mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
       } else {
           Paso_SystemMatrix_nullifyRowsAndCols_CSR(A,mask_row,mask_col,main_diagonal_value);
       }
     }
  }
  return;
}

void Paso_SystemMatrix_nullifyRowsAndCols_CSC_BLK1(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol;
  #pragma omp parallel for private(irow, iptr,icol) schedule(static)
  for (icol=0;icol< A->pattern->myNumOutput;icol++) {
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
void Paso_SystemMatrix_nullifyRowsAndCols_CSR_BLK1(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol;
  #pragma omp parallel for private(irow, iptr,icol) schedule(static)
  for (irow=0;irow< A->pattern->myNumOutput;irow++) {
      /* TODO: test mask_row here amd not inside every loop */
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
void Paso_SystemMatrix_nullifyRowsAndCols_CSC(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t ir,icol,iptr,icb,irb,irow,ic,l;
  #pragma omp parallel for private(l,irow, iptr,icol,ic,irb,icb) schedule(static)
  for (ic=0;ic< A->pattern->myNumOutput;ic++) {
	for (iptr=A->pattern->ptr[ic]-index_offset;iptr<A->pattern->ptr[ic+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*(A->pattern->index[iptr]-index_offset);
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
void Paso_SystemMatrix_nullifyRowsAndCols_CSR(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t ir,icol,iptr,icb,irb,irow,ic,l;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb) schedule(static)
  for (ir=0;ir< A->pattern->myNumOutput;ir++) {
	for (iptr=A->pattern->ptr[ir]-index_offset;iptr<A->pattern->ptr[ir+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
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
void Paso_SystemMatrix_nullifyRows_CSR_BLK1(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t irow, iptr, icol_global, irow_global;
  index_t myFirstRow=A->myFirstRow;
  #pragma omp parallel for private(irow, iptr,icol_global, irow_global) schedule(static)
  for (irow=0;irow< A->pattern->myNumOutput;irow++) {
      if (mask_row[irow]>0.) {
         irow_global=myFirstRow+irow;
         for (iptr=A->pattern->ptr[irow]-index_offset;iptr<A->pattern->ptr[irow+1]-index_offset; iptr++) {
           icol_global=A->pattern->index[iptr]-index_offset;
           if (irow_global==icol_global) {
	      A->val[iptr]=main_diagonal_value;
            } else {
	      A->val[iptr]=0;
            }
	 }
     }
  } 
}
void Paso_SystemMatrix_nullifyRows_CSR(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t myFirstRow=A->myFirstRow;
  index_t ir,icol,iptr,icb,irb,irow,ic,l, irow_global, icol_global;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb, irow_global, icol_global) schedule(static)
  for (ir=0;ir< A->pattern->myNumOutput;ir++) {
	for (iptr=A->pattern->ptr[ir]-index_offset;iptr<A->pattern->ptr[ir+1]-index_offset; iptr++) {
	  for (irb=0;irb< A->row_block_size;irb++) {
	    irow=irb+A->row_block_size*ir;
	    if (mask_row[irow]>0. ) {
               irow_global=irow+myFirstRow*A->row_block_size;
	       for (icb=0;icb< A->col_block_size;icb++) {
	           icol_global=icb+A->col_block_size*(A->pattern->index[iptr]-index_offset);
                   l=iptr*A->block_size+irb+A->row_block_size*icb;
	           if (irow_global==icol_global) {
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
void Paso_SystemMatrix_nullifyCols_CSR_BLK1(Paso_SystemMatrix* A,
                                              double* mask_col, 
                                              double main_diagonal_value,
                                              index_t min_index, index_t max_index) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t myFirstRow=A->myFirstRow;
  index_t irow, iptr, icol, icol_global, irow_global;
  #pragma omp parallel for private(irow, iptr,icol, icol_global, irow_global) schedule(static)
  for (irow=0;irow< A->pattern->myNumOutput;irow++) {
      irow_global=irow+myFirstRow;
      for (iptr=A->pattern->ptr[irow]-index_offset;iptr<A->pattern->ptr[irow+1]-index_offset; iptr++) {
         icol_global=A->pattern->index[iptr]-index_offset;
         if (min_index <= icol_global && icol_global < max_index) {
            icol=icol_global-min_index;
            if (mask_col[icol]>0.) {
               if (irow_global==icol_global) {
	          A->val[iptr]=main_diagonal_value;
                } else {
	          A->val[iptr]=0;
                }
	    }
         }
      } 
  }
}
void Paso_SystemMatrix_nullifyCols_CSR(Paso_SystemMatrix* A,
                                              double* mask_col, 
                                              double main_diagonal_value,
                                              index_t min_index, index_t max_index) {
  index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t myFirstRow=A->myFirstRow;
  index_t ir,icol,iptr,icb,irb,irow,ic,l, icol_global, irow_global, ic_global;
  #pragma omp parallel for private(l,irow, iptr,icol,ir,irb,icb, icol_global, irow_global, ic_global) schedule(static)
  for (ir=0;ir< A->pattern->myNumOutput;ir++) {
      for (iptr=A->pattern->ptr[ir]-index_offset;iptr<A->pattern->ptr[ir+1]-index_offset; iptr++) {
	  ic_global=A->pattern->index[iptr]-index_offset;
          if (min_index <= ic_global && ic_global < max_index) {
	     for (icb=0;icb< A->col_block_size;icb++) {
	         icol_global=icb+A->col_block_size*ic_global;
	         icol=icb+A->col_block_size*(ic_global-min_index);
	         if (mask_col[icol]>0.) {
	            for (irb=0;irb< A->row_block_size;irb++) {
	               irow_global=irb+A->row_block_size*(ir+myFirstRow);
                       l=iptr*A->block_size+irb+A->row_block_size*icb;
		       if (irow_global==icol_global) {
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
}
