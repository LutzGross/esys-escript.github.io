
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix                                       */

/*  Merge the MainBlock and CoupleBlock in the matrix         */


/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_mergeMainAndCouple(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val){
  if (A->type & MATRIX_FORMAT_CSC) {
      if (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) {
          Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1(A, p_ptr, p_idx, p_val);
      } else {
          Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: CSC with index 0 or block size larger than 1 is not supported. ");
      }
  } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
      Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: TRILINOS is not supported. ");
  } else {
      Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: CRS is not supported. ");
  }
  return;
}

void Paso_SystemMatrix_Coupler_collectGlobalID(Paso_Coupler* coupler, index_t **global, index_t offset) {
  Paso_MPIInfo *mpi_info = coupler->mpi_info;
  index_t *recv_buffer=NULL;
  index_t *send_buffer=NULL;
  dim_t i;
  
  if (mpi_info->size <= 1) {
      *global = NULL;
      return;
  }

  /* initialization */
  send_buffer = MEMALLOC(coupler->connector->send->numSharedComponents, index_t);
  recv_buffer = MEMALLOC(coupler->connector->recv->numSharedComponents, index_t);

  /* start receiving input */
  for (i=0; i< coupler->connector->recv->numNeighbors; ++i) {
      #ifdef PASO_MPI
      MPI_Irecv(&(recv_buffer[coupler->connector->recv->offsetInShared[i]]),
                coupler->connector->recv->offsetInShared[i+1]- coupler->connector->recv->offsetInShared[i],
                MPI_INT,
                coupler->connector->recv->neighbor[i],
                mpi_info->msg_tag_counter+coupler->connector->recv->neighbor[i],
                mpi_info->comm,
                &(coupler->mpi_requests[i]));
      #endif
  }

  /* collect values into buffer */
  #pragma omp parallel for private(i)
  for (i=0; i < coupler->connector->send->numSharedComponents;++i) 
      send_buffer[i]=coupler->connector->send->shared[i] + offset;

  /* send buffer out */
  for (i=0; i< coupler->connector->send->numNeighbors; ++i) {
      #ifdef PASO_MPI
      MPI_Issend(&(send_buffer[coupler->connector->send->offsetInShared[i]]),
                 coupler->connector->send->offsetInShared[i+1]- coupler->connector->send->offsetInShared[i],
                 MPI_INT,
                 coupler->connector->send->neighbor[i],
                 mpi_info->msg_tag_counter+mpi_info->rank,
                 mpi_info->comm,
                 &(coupler->mpi_requests[i+ coupler->connector->recv->numNeighbors]));
      #endif
  }

  mpi_info->msg_tag_counter+=mpi_info->size;

  /* wait for receive */
  #ifdef PASO_MPI
  MPI_Waitall(coupler->connector->recv->numNeighbors+coupler->connector->send->numNeighbors,
              coupler->mpi_requests,
              coupler->mpi_stati);
  #endif
  
  /* finalization */
  *global = recv_buffer;
  MEMFREE(send_buffer);

  return;
}

void Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val) {

  index_t i, j, i_ub, j_ub, col, num_vals, main_num_vals;
  index_t couple_num_vals, idx, rank, main_offset;
  index_t main_num_cols=A->mainBlock->pattern->numOutput;
  index_t couple_num_cols=A->col_coupleBlock->pattern->numOutput;
  index_t *main_ptr=A->mainBlock->pattern->ptr;
  index_t *main_idx=A->mainBlock->pattern->index;
  double  *main_val=A->mainBlock->val;
  index_t *couple_ptr=A->col_coupleBlock->pattern->ptr;
  index_t *couple_idx=A->col_coupleBlock->pattern->index;
  double  *couple_val=A->col_coupleBlock->val;
  index_t *couple_global=NULL;
  Paso_Coupler* coupler=A->col_coupler;

fprintf(stderr, "CHP 0\n");


  if (A->mainBlock->col_block_size!=1 || 
      A->mainBlock->row_block_size!=1 ||
      A->col_coupleBlock->col_block_size!=1 ||
      A->col_coupleBlock->row_block_size!=1) {
      Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1: requires format with block size 1.");
      return;
  }
  
  if (main_num_cols != couple_num_cols) {
      Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1: number of collums do not match.");
      return;
  }

  /* allocate arrays "ptr", "index" and "val" */
  main_num_vals = main_ptr[main_num_cols]-1;
  couple_num_vals = couple_ptr[couple_num_cols]-1;
  num_vals = main_num_vals + couple_num_vals;
  *p_ptr = MEMALLOC(main_num_cols+1, index_t);
  *p_idx = MEMALLOC(num_vals, index_t);
  *p_val = MEMALLOC(num_vals, double);

  /* initialize before merge */
  (*p_ptr)[0] = 1;
  rank = A->mpi_info->rank;
  main_offset = A->col_distribution->first_component[rank];
  i=0;
  j=0;

  /* gather global index for entries in couple block */
  Paso_SystemMatrix_Coupler_collectGlobalID(coupler, &couple_global, main_offset+1);

  /* merge mainBlock and col_coupleBlock */
  for (col=1; col<=main_num_cols; col++) {
      i_ub = main_ptr[col]-1;
      j_ub = couple_ptr[col]-1;
      while (i < i_ub || j < j_ub) {
          if (j < j_ub) {
              /* switch from coupleBlock index to global row index of matrix */
              if (A->mpi_info->size == 1) {
                 Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1: requires more than 1 MPI rank when coupleBlock exists.");
              }

              idx = couple_global[couple_idx[j]-1];
          }
          if (j == j_ub || (i < i_ub && (main_idx[i] + main_offset) < idx)){
              (*p_idx)[i+j] = main_idx[i] + main_offset;
              (*p_val)[i+j] = main_val[i];
              i++;
          } else {
              (*p_idx)[i+j] = idx;
              (*p_val)[i+j] = couple_val[j];
              j++;
          }
      }
      (*p_ptr)[col] = i+j+1;
  }

  MEMFREE(couple_global);
  return;
}

void Paso_SystemMatrix_copyMain_CSC_OFFSET1(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val) {

  index_t i, i_ub, col, num_vals, idx;
  index_t main_num_cols=A->mainBlock->pattern->numOutput;
  index_t *main_ptr=A->mainBlock->pattern->ptr;
  index_t *main_idx=A->mainBlock->pattern->index;
  double  *main_val=A->mainBlock->val;

  if (A->mainBlock->col_block_size!=1 || A->mainBlock->row_block_size!=1) {
      Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1: requires format with block size 1.");
      return;
  }

  /* allocate arrays "ptr", "index" and "val" */
  num_vals = main_ptr[main_num_cols]-1;
  *p_ptr = MEMALLOC(main_num_cols+1, index_t);
  *p_idx = MEMALLOC(num_vals, index_t);
  *p_val = MEMALLOC(num_vals, double);

  /* copy from mainBlock */
  (*p_ptr)[0] = 1;
  i=0;
  for (col=1; col<=main_num_cols; col++) {
      i_ub = main_ptr[col]-1;
      while (i < i_ub) {
          (*p_idx)[i] = main_idx[i];
          (*p_val)[i] = main_val[i];
          i++;
      }
      (*p_ptr)[col] = i+1;
  }
  return;
}

