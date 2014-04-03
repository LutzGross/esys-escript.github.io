
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/
/* Paso: SystemMatrix							*/
/*									*/
/*  Copy mainBlock and col_coupleBlock in other ranks			*/
/*  into remote_coupleBlock						*/
/*									*/
/*  WARNING: function uses mpi_requests of the coupler attached to A.	*/
/*									*/
/****************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"

void Paso_SystemMatrix_copyRemoteCoupleBlock(Paso_SystemMatrix* A, bool recreatePattern)
{
  paso::Coupler_ptr coupler;
  paso::SharedComponents_ptr send, recv;
  double *cols=NULL, *send_buf=NULL;
  index_t *global_id=NULL, *cols_array=NULL, *ptr_ptr=NULL, *ptr_idx=NULL;
  index_t *send_idx=NULL, *send_offset=NULL, *recv_buf=NULL, *recv_offset=NULL;
  index_t i, j, k, l, m, n, p, q, j_ub, k_lb, k_ub, l_lb, l_ub, i0;
  index_t offset, len, overlapped_n, base, block_size, block_size_size;
  index_t num_main_cols, num_couple_cols, num_neighbors, row, neighbor;
  index_t *recv_degree=NULL, *send_degree=NULL;
  dim_t rank=A->mpi_info->rank, mpi_size=A->mpi_info->size;

  if (mpi_size == 1)
      return;

  if (recreatePattern)
      A->remote_coupleBlock.reset();

  if (A->remote_coupleBlock)
      return;

  /* sending/receiving unknown's global ID */
  num_main_cols = A->mainBlock->numCols;
  cols = new double[num_main_cols];
  offset = A->col_distribution->first_component[rank];
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<num_main_cols; ++i) cols[i] = offset + i;
  if (A->global_id == NULL) {
    coupler.reset(new paso::Coupler(A->col_coupler->connector, 1));
    coupler->startCollect(cols);
  }

  recv_buf = new index_t[mpi_size];
  recv_degree = new index_t[mpi_size];
  recv_offset = new index_t[mpi_size+1];
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<mpi_size; i++){
    recv_buf[i] = 0;
    recv_degree[i] = 1;
    recv_offset[i] = i;
  }

  num_couple_cols = A->col_coupleBlock->numCols;
  overlapped_n = A->row_coupleBlock->numRows;
  send = A->row_coupler->connector->send;
  recv = A->row_coupler->connector->recv;
  num_neighbors = send->numNeighbors;
  block_size = A->block_size;
  block_size_size = block_size * sizeof(double);

  /* waiting for receiving unknown's global ID */
  if (A->global_id == NULL) {
        coupler->finishCollect();
        global_id = new index_t[num_couple_cols+1];
#pragma omp parallel for private(i) schedule(static)
        for (i=0; i<num_couple_cols; ++i) 
	    global_id[i] = coupler->recv_buffer[i];
        A->global_id = global_id;
        coupler.reset();
  } else 
  global_id = A->global_id;

  /* distribute the number of cols in current col_coupleBlock to all ranks */
  #ifdef ESYS_MPI
  MPI_Allgatherv(&num_couple_cols, 1, MPI_INT, recv_buf, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm); 
  #endif

  /* distribute global_ids of cols to be considered to all ranks*/
  len = 0;
  for (i=0; i<mpi_size; i++){
    recv_degree[i] = recv_buf[i];
    recv_offset[i] = len;
    len += recv_buf[i];
  }
  recv_offset[mpi_size] = len;
  cols_array = new index_t[len];

  #ifdef ESYS_MPI
  MPI_Allgatherv(global_id, num_couple_cols, MPI_INT, cols_array, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm);
  #endif

  /* first, prepare the ptr_ptr to be received */
  ptr_ptr = new index_t[overlapped_n+1];
  for (p=0; p<recv->numNeighbors; p++) {
    row = recv->offsetInShared[p];
    i = recv->offsetInShared[p+1];
    #ifdef ESYS_MPI
    MPI_Irecv(&(ptr_ptr[row]), i-row, MPI_INT, recv->neighbor[p], 
		A->mpi_info->msg_tag_counter+recv->neighbor[p],
		A->mpi_info->comm,
		&(A->row_coupler->mpi_requests[p]));
    #endif
  }

  /* now prepare the rows to be sent (the degree, the offset and the data) */
  p = send->offsetInShared[num_neighbors];
  len = 0;
  for (i=0; i<num_neighbors; i++) {
    /* #cols per row X #rows */  
    len += recv_buf[send->neighbor[i]] * 
		(send->offsetInShared[i+1] - send->offsetInShared[i]);
  }
  send_buf = new double[len*block_size];
  send_idx = new index_t[len];
  send_offset = new index_t[p+1];
  send_degree = new index_t[num_neighbors];

  len = 0;
  base = 0;
  i0 = 0;
  for (p=0; p<num_neighbors; p++) {
    i = i0;
    neighbor = send->neighbor[p];
    l_ub = recv_offset[neighbor+1];
    l_lb = recv_offset[neighbor];
    j_ub = send->offsetInShared[p + 1];
    for (j=send->offsetInShared[p]; j<j_ub; j++) {
	row = send->shared[j];

        /* check col_coupleBlock for data to be passed @row */
        l = l_lb;
        k_ub = A->col_coupleBlock->pattern->ptr[row+1];
        k = A->col_coupleBlock->pattern->ptr[row];
	q = A->mainBlock->pattern->index[A->mainBlock->pattern->ptr[row]] + offset;
        while (k<k_ub && l<l_ub) {
          m = global_id[A->col_coupleBlock->pattern->index[k]];
	  if (m > q) break;
          n = cols_array[l];
          if (m == n) {
            send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->col_coupleBlock->val[block_size*k]), block_size_size);
            len++;
            l++;
            k++;
          } else if (m < n) {
            k++;
          } else {
            l++;
          }
        }
	k_lb = k;

	/* check mainBlock for data to be passed @row */
	k_ub = A->mainBlock->pattern->ptr[row+1];
	k=A->mainBlock->pattern->ptr[row]; 
	while (k<k_ub && l<l_ub) {
	  m = A->mainBlock->pattern->index[k] + offset;
	  n = cols_array[l];
	  if (m == n) {
	    send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->mainBlock->val[block_size*k]), block_size_size);
	    len++;
	    l++; 
	    k++;
	  } else if (m < n) {
	    k++;
	  } else {
	    l++;
	  }
	} 

	/* check col_coupleBlock for data to be passed @row */
	k_ub = A->col_coupleBlock->pattern->ptr[row+1];
	k=k_lb;
	while (k<k_ub && l<l_ub) {
	  m = global_id[A->col_coupleBlock->pattern->index[k]];
	  n = cols_array[l];
	  if (m == n) {
	    send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->col_coupleBlock->val[block_size*k]), block_size_size);
	    len++;
	    l++;
	    k++;
	  } else if (m < n) {
	    k++;
	  } else {
	    l++;
	  }
	}

	send_offset[i] = len - base;
	base = len;
	i++;
    }

    /* sending */
    #ifdef ESYS_MPI
    MPI_Issend(&(send_offset[i0]), i-i0, MPI_INT, send->neighbor[p],
		A->mpi_info->msg_tag_counter+rank,
		A->mpi_info->comm,
		&(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    #endif
    send_degree[p] = len;
    i0 = i;
  }

  #ifdef ESYS_MPI
  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
		A->row_coupler->mpi_requests,
		A->row_coupler->mpi_stati);
  #endif
  ESYS_MPI_INC_COUNTER(*(A->mpi_info), mpi_size)

  len = 0;
  for (i=0; i<overlapped_n; i++) {
    p = ptr_ptr[i];
    ptr_ptr[i] = len;
    len += p;
  }
  ptr_ptr[overlapped_n] = len;
  ptr_idx = new index_t[len];

  /* send/receive index array */
  j=0;
  for (p=0; p<recv->numNeighbors; p++) {
    i = ptr_ptr[recv->offsetInShared[p+1]] - ptr_ptr[recv->offsetInShared[p]];
    #ifdef ESYS_MPI
    if (i > 0) 
	MPI_Irecv(&(ptr_idx[j]), i, MPI_INT, recv->neighbor[p],
                A->mpi_info->msg_tag_counter+recv->neighbor[p],
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p]));
    #endif
    j += i;
  }

  j=0;
  for (p=0; p<num_neighbors; p++) {
    i = send_degree[p] - j;
    #ifdef ESYS_MPI
    if (i > 0) 
	MPI_Issend(&(send_idx[j]), i, MPI_INT, send->neighbor[p],
                A->mpi_info->msg_tag_counter+rank,
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    #endif
    j = send_degree[p];
  }

  #ifdef ESYS_MPI
  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
                A->row_coupler->mpi_requests,
                A->row_coupler->mpi_stati);
  #endif
  ESYS_MPI_INC_COUNTER(*(A->mpi_info), mpi_size)

  /* allocate pattern and sparsematrix for remote_coupleBlock */
  paso::Pattern_ptr pattern(new paso::Pattern(A->row_coupleBlock->pattern->type,
                overlapped_n, num_couple_cols, ptr_ptr, ptr_idx));
  A->remote_coupleBlock.reset(new paso::SparseMatrix(A->row_coupleBlock->type,
                pattern, A->row_block_size, A->col_block_size, false));

  /* send/receive value array */
  j=0;
  for (p=0; p<recv->numNeighbors; p++) {
    i = ptr_ptr[recv->offsetInShared[p+1]] - ptr_ptr[recv->offsetInShared[p]];
    #ifdef ESYS_MPI
    if (i > 0) 
	MPI_Irecv(&(A->remote_coupleBlock->val[j]), i * block_size, 
		MPI_DOUBLE, recv->neighbor[p],
                A->mpi_info->msg_tag_counter+recv->neighbor[p],
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p]));
    #endif
    j += (i * block_size);
  }

  j=0;
  for (p=0; p<num_neighbors; p++) {
    i = send_degree[p] - j;
    #ifdef ESYS_MPI
    if (i > 0)
	MPI_Issend(&(send_buf[j*block_size]), i*block_size, MPI_DOUBLE, send->neighbor[p],
                A->mpi_info->msg_tag_counter+rank,
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    #endif
    j = send_degree[p];
  }

  #ifdef ESYS_MPI
  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
                A->row_coupler->mpi_requests,
                A->row_coupler->mpi_stati);
  #endif
  ESYS_MPI_INC_COUNTER(*(A->mpi_info), mpi_size)

  /* release all temp memory allocation */
  delete[] cols;
  delete[] cols_array;
  delete[] recv_offset;
  delete[] recv_degree;
  delete[] recv_buf;
  delete[] send_buf;
  delete[] send_offset;
  delete[] send_degree;
  delete[] send_idx;
}

