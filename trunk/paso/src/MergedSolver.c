
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: AMG preconditioner  (local version)                  */

/**************************************************************/

/* Author: lgao@uq.edu.au, l.gross@uq.edu.au                 */

/**************************************************************/


#include "Paso.h"
#include "MergedSolver.h"
#include "Preconditioner.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "MKL.h"
#include<stdio.h>

Paso_MergedSolver* Paso_MergedSolver_alloc(Paso_SystemMatrix *A, Paso_Options* options)
{
   const index_t rank = A->mpi_info->rank;
   const index_t size = A->mpi_info->size;
   const dim_t global_n = Paso_SystemMatrix_getGlobalNumRows(A);
   const dim_t n_block = A->mainBlock->row_block_size;
   const dim_t* dist = A->pattern->input_distribution->first_component;
   dim_t i;
   Paso_SparseMatrix *A_temp;

   Paso_MergedSolver* out = NULL;
   A_temp = Paso_MergedSolver_mergeSystemMatrix(A); 
   out = MEMALLOC(1, Paso_MergedSolver);
   Esys_checkPtr(out);
   if (Esys_noError()) {
   
       /* merge matrix */
       out->A = NULL;
       out->mpi_info=Esys_MPIInfo_getReference(A->mpi_info);
       out->reordering = options->reordering;
       out->refinements = options->coarse_matrix_refinements;
       /*out->verbose = options->verbose; */
       out->verbose = FALSE;
       out->sweeps = options->pre_sweeps+options->post_sweeps;

       /* First, gather x and b into rank 0 */
       out->b = TMPMEMALLOC(global_n*n_block, double);
       out->x = TMPMEMALLOC(global_n*n_block, double);
       out->counts = TMPMEMALLOC(size, index_t);
       out->offset = TMPMEMALLOC(size, index_t);
       
       if (! (Esys_checkPtr(out->b) || Esys_checkPtr(out->x) || Esys_checkPtr(out->counts) || Esys_checkPtr(out->offset) ) ){
           #pragma omp parallel for private(i)
           for (i=0; i<size; i++) {
             const dim_t p = dist[i];
             out->counts[i] = (dist[i+1] - p)*n_block;
             out->offset[i] = p*n_block;
           }

          if (rank == 0) {
             /* solve locally */
             #ifdef MKL
               out->A = Paso_SparseMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, A_temp);
               out->A->solver_package = PASO_MKL;
             #else 
               #ifdef UMFPACK
	         out->A  = Paso_SparseMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_CSC, A_temp);
	         out->A->solver_package = PASO_UMFPACK;
               #else
	         out->A->solver_p = Paso_Preconditioner_LocalSmoother_alloc(out->A, (options->smoother == PASO_JACOBI), out->verbose);
	         out->A->solver_package = PASO_SMOOTHER;
               #endif
             #endif
           }
       }
   }

   Paso_SparseMatrix_free(A_temp);
   if ( Esys_noError()) {
      return out;
   } else {
     Paso_MergedSolver_free(out); 
      return NULL;
   }
}

void Paso_MergedSolver_free(Paso_MergedSolver* in) {
     if (in!=NULL) {
        Paso_SparseMatrix_free(in->A);
	MEMFREE(in->x);
	MEMFREE(in->b);
	MEMFREE(in->counts);
	MEMFREE(in->offset);
	MEMFREE(in);
     }
}

/* Merge the system matrix which is distributed on ranks into a complete 
   matrix on rank 0, then solve this matrix on rank 0 only */
Paso_SparseMatrix* Paso_MergedSolver_mergeSystemMatrix(Paso_SystemMatrix* A) {
  index_t i, iptr, j, n, remote_n, global_n, len, offset, tag;
  index_t row_block_size, col_block_size, block_size;
  index_t size=A->mpi_info->size;
  index_t rank=A->mpi_info->rank;
  index_t *ptr=NULL, *idx=NULL, *ptr_global=NULL, *idx_global=NULL;
  index_t *temp_n=NULL, *temp_len=NULL;
  double  *val=NULL;
  Paso_Pattern *pattern=NULL;
  Paso_SparseMatrix *out=NULL;
  #ifdef ESYS_MPI
    MPI_Request* mpi_requests=NULL;
    MPI_Status* mpi_stati=NULL;
  #else
    int *mpi_requests=NULL, *mpi_stati=NULL;
  #endif

  if (size == 1) {
    n = A->mainBlock->numRows;
    ptr = TMPMEMALLOC(n, index_t); 
    #pragma omp parallel for private(i)
    for (i=0; i<n; i++) ptr[i] = i;
    out = Paso_SparseMatrix_getSubmatrix(A->mainBlock, n, n, ptr, ptr);
    TMPMEMFREE(ptr);
    return out;
  }

  n = A->mainBlock->numRows;
  block_size = A->block_size;

  /* Merge MainBlock and CoupleBlock to get a complete column entries
     for each row allocated to current rank. Output (ptr, idx, val) 
     contains all info needed from current rank to merge a system 
     matrix  */
  Paso_SystemMatrix_mergeMainAndCouple(A, &ptr, &idx, &val);

  #ifdef ESYS_MPI
    mpi_requests=TMPMEMALLOC(size*2,MPI_Request);
    mpi_stati=TMPMEMALLOC(size*2,MPI_Status);
  #else
    mpi_requests=TMPMEMALLOC(size*2,int);
    mpi_stati=TMPMEMALLOC(size*2,int);
  #endif

  /* Now, pass all info to rank 0 and merge them into one sparse 
     matrix */
  if (rank == 0) {
    /* First, copy local ptr values into ptr_global */
    global_n=Paso_SystemMatrix_getGlobalNumRows(A);
    ptr_global = MEMALLOC(global_n+1, index_t);
    memcpy(ptr_global, ptr, (n+1) * sizeof(index_t));
    iptr = n+1;
    MEMFREE(ptr);
    temp_n = TMPMEMALLOC(size, index_t);
    temp_len = TMPMEMALLOC(size, index_t);
    temp_n[0] = iptr;
    
    /* Second, receive ptr values from other ranks */
    for (i=1; i<size; i++) {
      remote_n = A->row_distribution->first_component[i+1] -
		 A->row_distribution->first_component[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(ptr_global[iptr]), remote_n, MPI_INT, i, 
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      #endif
      temp_n[i] = remote_n;
      iptr += remote_n;
    }
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;

    /* Then, prepare to receive idx and val from other ranks */
    len = 0;
    offset = -1;
    for (i=0; i<size; i++) {
      if (temp_n[i] > 0) {
	offset += temp_n[i];
	len += ptr_global[offset];
	temp_len[i] = ptr_global[offset];
      }else 
	temp_len[i] = 0;
    }

    idx_global = MEMALLOC(len, index_t);
    iptr = temp_len[0];
    offset = n+1;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(idx_global[iptr]), len, MPI_INT, i,
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      #endif
      remote_n = temp_n[i];
      for (j=0; j<remote_n; j++) {
	ptr_global[j+offset] = ptr_global[j+offset] + iptr;
      }
      offset += remote_n;
      iptr += len;
    }
    memcpy(idx_global, idx, temp_len[0] * sizeof(index_t));
    MEMFREE(idx);
    row_block_size = A->mainBlock->row_block_size;
    col_block_size = A->mainBlock->col_block_size;
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_n);

    /* Then generate the sparse matrix */
    pattern = Paso_Pattern_alloc(A->mainBlock->pattern->type, global_n,
			global_n, ptr_global, idx_global);
    out = Paso_SparseMatrix_alloc(A->mainBlock->type, pattern, 
			row_block_size, col_block_size, FALSE);
    Paso_Pattern_free(pattern);

    /* Finally, receive and copy the value */
    iptr = temp_len[0] * block_size;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(out->val[iptr]), len * block_size, MPI_DOUBLE, i,
                        A->mpi_info->msg_tag_counter+i,
                        A->mpi_info->comm,
                        &mpi_requests[i]);
      #endif
      iptr += (len * block_size);
    }
    memcpy(out->val, val, temp_len[0] * sizeof(double) * block_size);
    MEMFREE(val);
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_len);
  } else { /* it's not rank 0 */

    /* First, send out the local ptr */
    tag = A->mpi_info->msg_tag_counter+rank;
    #ifdef ESYS_MPI
    MPI_Issend(&(ptr[1]), n, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[0]);
    #endif

    /* Next, send out the local idx */
    len = ptr[n];
    tag += size;
    #ifdef ESYS_MPI
    MPI_Issend(idx, len, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[1]);
    #endif

    /* At last, send out the local val */
    len *= block_size;
    tag += size;
    #ifdef ESYS_MPI
    MPI_Issend(val, len, MPI_DOUBLE, 0, tag, A->mpi_info->comm, 
			&mpi_requests[2]);

    MPI_Waitall(3, mpi_requests, mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter = tag + size - rank;
    MEMFREE(ptr);
    MEMFREE(idx);
    MEMFREE(val);

    out = NULL;
  }

  TMPMEMFREE(mpi_requests);
  TMPMEMFREE(mpi_stati);
  return out;
}


void Paso_MergedSolver_solve(Paso_MergedSolver* ms, double* local_x, double* local_b) {

  const index_t rank = ms->mpi_info->rank;
  const dim_t count = ms->counts[rank];


  #ifdef ESYS_MPI
  {
     MPI_Gatherv(local_b, count, MPI_DOUBLE, ms->b, ms->counts, ms->offset, MPI_DOUBLE, 0, ms->mpi_info->comm);
  }
  #else
  {
      dim_t i;
      #pragma omp parallel for private(i)
      for (i=0; i<count; i++) {
             ms->b[i]=local_b[i];
             ms->x[i]=local_x[i];
      }
  }
  #endif
  if (rank == 0) {
            switch (ms->A->solver_package) {
               case (PASO_MKL):
                  Paso_MKL(ms->A, ms->x,ms->b, ms->reordering, ms->refinements, ms->verbose);
                  break;
               case (PASO_UMFPACK):
                  Paso_UMFPACK(ms->A, ms->x,ms->b, ms->refinements, ms->verbose);
                  break;
               case (PASO_SMOOTHER):
                  Paso_Preconditioner_LocalSmoother_solve(ms->A, ms->A->solver_p,ms->x,ms->b,ms->sweeps, FALSE);
                  break;
            }
  }
  #ifdef ESYS_MPI
  /* now we need to distribute the solution to all ranks */
  {
     MPI_Scatterv(ms->x, ms->counts, ms->offset, MPI_DOUBLE, local_x, count, MPI_DOUBLE, 0, ms->mpi_info->comm);
  }
  #else
  {
      dim_t i;
      #pragma omp parallel for private(i)
      for (i=0; i<count; i++) local_x[i]=ms->x[i];
  }
  #endif
}
