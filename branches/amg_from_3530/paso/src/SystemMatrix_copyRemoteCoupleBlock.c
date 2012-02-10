
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


/************************************************************************/
/* Paso: SystemMatrix							*/
/*									*/
/*  Copy mainBlock and col_coupleBlock in other ranks			*/
/*  into remote_coupleBlock						*/
/*									*/
/*  WARNING: function uses mpi_reqests of the coupler attached to A.	*/
/*									*/
/************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MY_DEBUG 0
#define MY_DEBUG1 0

void Paso_SystemMatrix_copyRemoteCoupleBlock(Paso_SystemMatrix* A, const bool_t recreatePattern){
  Paso_Pattern *pattern=NULL;
  Paso_Coupler *coupler=NULL;
  Paso_SharedComponents *send=NULL, *recv=NULL;
  double *cols=NULL, *send_buf=NULL;
  index_t *global_id=NULL, *cols_array=NULL, *ptr_ptr=NULL, *ptr_idx=NULL;
  index_t *send_idx=NULL, *send_offset=NULL, *recv_buf=NULL, *recv_offset=NULL;
  index_t i, j, k, l, m, n, p, q, j_ub, k_lb, k_ub, l_lb, l_ub, i0;
  index_t offset, len, overlapped_n, base, block_size, block_size_size;
  index_t num_main_cols, num_couple_cols, num_neighbors, row, neighbor;
  index_t *recv_degree=NULL, *send_degree=NULL;
  dim_t rank=A->mpi_info->rank, mpi_size=A->mpi_info->size;

  if (mpi_size == 1) return;

  if (recreatePattern) Paso_SparseMatrix_free(A->remote_coupleBlock);

  if (A->remote_coupleBlock) return;
if (MY_DEBUG) fprintf(stderr, "rank%d XCP1\n", rank);

  /* sending/receiving unknown's global ID */
  num_main_cols = A->mainBlock->numCols;
  cols = TMPMEMALLOC(num_main_cols, double);
  offset = A->col_distribution->first_component[rank];
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<num_main_cols; ++i) cols[i] = offset + i;
  if (A->global_id == NULL) {
    coupler=Paso_Coupler_alloc(A->col_coupler->connector, 1);
    Paso_Coupler_startCollect(coupler, cols);
  }

  recv_buf = TMPMEMALLOC(mpi_size, index_t);
  recv_degree = TMPMEMALLOC(mpi_size, index_t);
  recv_offset = TMPMEMALLOC(mpi_size+1, index_t);
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
    Paso_Coupler_finishCollect(coupler);
    global_id = MEMALLOC(num_couple_cols+1, index_t);
//fprintf(stderr, "rank %d global_id just allocated %d. index_t %d\n", rank, num_couple_cols, sizeof(index_t));
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_couple_cols; ++i) 
	global_id[i] = coupler->recv_buffer[i];
    A->global_id = global_id;
    Paso_Coupler_free(coupler);
  } else 
//fprintf(stderr, "rank %d global_id already allocated. index_t %d\n", rank, sizeof(index_t));
  global_id = A->global_id;
if (MY_DEBUG) fprintf(stderr, "rank%d XCP2\n", rank);
if (MY_DEBUG && rank == 1){
  int p, sum;
  char *str1, *str2;
  sum = num_couple_cols+1;
  str1 = TMPMEMALLOC(sum*100+100, char);
  str2 = TMPMEMALLOC(100, char);
  sprintf(str1, "global_id[%d]= ", sum);
  for (p=0; p<sum; ++p) {
    sprintf(str2, "%d ", global_id[p]);
    strcat(str1, str2);
    if ( (p+1) == ((p+1)/300)*300) sprintf(str1, "%s\n", str1);
  }
  fprintf(stderr, "%s\n", str1);
  TMPMEMFREE(str1);
  TMPMEMFREE(str2);
}

  /* distribute the number of cols in current col_coupleBlock to all ranks */
  MPI_Allgatherv(&num_couple_cols, 1, MPI_INT, recv_buf, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm); 

  /* distribute global_ids of cols to be considered to all ranks*/
  len = 0;
  for (i=0; i<mpi_size; i++){
    recv_degree[i] = recv_buf[i];
    recv_offset[i] = len;
    len += recv_buf[i];
  }
  recv_offset[mpi_size] = len;
  cols_array = TMPMEMALLOC(len, index_t);
  if (Esys_checkPtr(cols_array)) fprintf(stderr, "rank %d MALLOC has trouble\n", rank);
if (MY_DEBUG)
fprintf(stderr, "rank %d num_couple_cols %d len %d (%d %d %d %d)\n", rank, num_couple_cols, len, recv_degree[0], recv_degree[1], recv_offset[0], recv_offset[1]);
/*{
char *str1, *str2;
int my_i, sum = num_couple_cols;
str1 = TMPMEMALLOC(sum*30+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d global_id[%d]=(", rank, sum);
for (my_i=0;my_i<sum;my_i++) {
  sprintf(str2, "%d-%d ", my_i, global_id[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}*/
  MPI_Allgatherv(global_id, num_couple_cols, MPI_INT, cols_array, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm);
/*{
char *str1, *str2;
int my_i, sum = len;
str1 = TMPMEMALLOC(sum*30+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d cols_array[%d]=(", rank, sum);
for (my_i=0;my_i<sum;my_i++) {
  sprintf(str2, "%d-%d ", my_i, cols_array[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}*/

if (MY_DEBUG) fprintf(stderr, "rank%d XCP3\n", rank);

  /* first, prepare the ptr_ptr to be received */
  ptr_ptr = MEMALLOC(overlapped_n+1, index_t);
  for (p=0; p<recv->numNeighbors; p++) {
    row = recv->offsetInShared[p];
    i = recv->offsetInShared[p+1];
    MPI_Irecv(&(ptr_ptr[row]), i-row, MPI_INT, recv->neighbor[p], 
		A->mpi_info->msg_tag_counter+recv->neighbor[p],
		A->mpi_info->comm,
		&(A->row_coupler->mpi_requests[p]));
  }

  /* now prepare the rows to be sent (the degree, the offset and the data) */
  p = send->offsetInShared[num_neighbors];
  len = 0;
  for (i=0; i<num_neighbors; i++) {
    /* #cols per row X #rows */  
    len += recv_buf[send->neighbor[i]] * 
		(send->offsetInShared[i+1] - send->offsetInShared[i]);
  }
  send_buf = TMPMEMALLOC(len*block_size, double);
  send_idx = TMPMEMALLOC(len, index_t);
  send_offset = TMPMEMALLOC(p+1, index_t);
  send_degree = TMPMEMALLOC(num_neighbors, index_t);

  len = 0;
  base = 0;
  i0 = 0;
  for (p=0; p<num_neighbors; p++) {
    i = i0;
    neighbor = send->neighbor[p];
    l_ub = recv_offset[neighbor+1];
    l_lb = recv_offset[neighbor];
    j_ub = send->offsetInShared[p + 1];
if (MY_DEBUG) fprintf(stderr, "rank%d: p%d neighbor%d l_ub%d l_lb%d j_ub%d %d %d\n", rank, p, neighbor, l_ub, l_lb, j_ub, send->offsetInShared[0], send->offsetInShared[1]);
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
//if (MY_DEBUG) fprintf(stderr, "rank%d: found in couple(1) %d\n", rank, m);
            send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->col_coupleBlock->val[block_size*k]), block_size_size);
//            send_buf[len] = A->col_coupleBlock->val[k]; /*??blocksize>1??*/
            len++;
            l++;
            k++;
          } else if (m < n) {
//if (MY_DEBUG) fprintf(stderr, "rank%d: checked in couple(1) %d %d\n", rank, m, n);
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
//if (MY_DEBUG) fprintf(stderr, "rank%d: found in main %d\n", rank, m);
	    send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->mainBlock->val[block_size*k]), block_size_size);
//	    send_buf[len] = A->mainBlock->val[k]; /*??blocksize>1??*/
	    len++;
	    l++; 
	    k++;
	  } else if (m < n) {
//if (MY_DEBUG) fprintf(stderr, "rank%d: checked in main %d %d\n", rank, m, n);
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
//if (MY_DEBUG) fprintf(stderr, "rank%d: found in couple(2) %d\n", rank, m);
	    send_idx[len] = l - l_lb;
	    memcpy(&(send_buf[len*block_size]), &(A->col_coupleBlock->val[block_size*k]), block_size_size);
//	    send_buf[len] = A->col_coupleBlock->val[k]; /*??blocksize>1??*/
	    len++;
	    l++;
	    k++;
	  } else if (m < n) {
//if (MY_DEBUG) fprintf(stderr, "rank%d: checked in couple(2) %d %d\n", rank, m, n);
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
    MPI_Issend(&(send_offset[i0]), i-i0, MPI_INT, send->neighbor[p],
		A->mpi_info->msg_tag_counter+rank,
		A->mpi_info->comm,
		&(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    send_degree[p] = len;
    i0 = i;
  }

  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
		A->row_coupler->mpi_requests,
		A->row_coupler->mpi_stati);
  A->mpi_info->msg_tag_counter += mpi_size;
if (MY_DEBUG) fprintf(stderr, "rank%d XCP4\n", rank);

  len = 0;
  for (i=0; i<overlapped_n; i++) {
    p = ptr_ptr[i];
    ptr_ptr[i] = len;
    len += p;
  }
if (MY_DEBUG) fprintf(stderr, "rank%d XCP4 cp1\n", rank);
  ptr_ptr[overlapped_n] = len;
if (MY_DEBUG) fprintf(stderr, "rank%d XCP4 cp2\n", rank);
  ptr_idx = MEMALLOC(len, index_t);
if (MY_DEBUG) fprintf(stderr, "rank%d XCP4_1 %d %d\n", rank, recv->numNeighbors, num_neighbors);

  /* send/receive index array */
  j=0;
  for (p=0; p<recv->numNeighbors; p++) {
    i = ptr_ptr[recv->offsetInShared[p+1]] - ptr_ptr[recv->offsetInShared[p]];
//if (i==0)
if (MY_DEBUG) fprintf(stderr, "rank%d p%d offset%d neighbor%d\n", rank, p, recv->offsetInShared[p+1], recv->numNeighbors);
    if (i > 0) 
	MPI_Irecv(&(ptr_idx[j]), i, MPI_INT, recv->neighbor[p],
                A->mpi_info->msg_tag_counter+recv->neighbor[p],
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p]));
    j += i;
  }
if (MY_DEBUG) fprintf(stderr, "rank%d XCP4_2\n", rank);

  j=0;
  for (p=0; p<num_neighbors; p++) {
    i = send_degree[p] - j;
//if (i==0)
if (MY_DEBUG) fprintf(stderr, "rank%d p%d send_degree%d neighbor%d\n", rank, p, send_degree[p], num_neighbors);
    if (i > 0) 
	MPI_Issend(&(send_idx[j]), i, MPI_INT, send->neighbor[p],
                A->mpi_info->msg_tag_counter+rank,
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    j = send_degree[p];
  }

if (MY_DEBUG) fprintf(stderr, "rank%d XCP4_3\n", rank);
  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
                A->row_coupler->mpi_requests,
                A->row_coupler->mpi_stati);
  A->mpi_info->msg_tag_counter += mpi_size;
if (MY_DEBUG) fprintf(stderr, "rank%d XCP5 %d %d\n", rank, overlapped_n, num_couple_cols);

if (MY_DEBUG1 && rank ==1) {
int my_i, sum1 = overlapped_n;
int sum = ptr_ptr[sum1];
char *str1, *str2;
fprintf(stderr, "rank %d %d\n", rank, sum);
str1 = TMPMEMALLOC(sum*10+100, char);
fprintf(stderr, "rank %d %d\n", rank, sum);
str2 = TMPMEMALLOC(30, char);
fprintf(stderr, "rank %d %d\n", rank, ptr_ptr[sum1]);
sprintf(str1, "rank %d ptr_ptr[%d] = (", rank, sum1);
for (my_i=0; my_i<sum1; my_i++){
  sprintf(str2, "%d ", ptr_idx[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
sprintf(str1, "rank %d ptr_idx[%d] = (", rank, sum);
for (my_i=0; my_i<sum; my_i++){
  sprintf(str2, "%d ", ptr_idx[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

  /* allocate pattern and sparsematrix for remote_coupleBlock */
  pattern = Paso_Pattern_alloc(A->row_coupleBlock->pattern->type,
                overlapped_n, num_couple_cols, ptr_ptr, ptr_idx);
  A->remote_coupleBlock = Paso_SparseMatrix_alloc(A->row_coupleBlock->type,
                pattern, A->row_block_size, A->col_block_size, 
                FALSE);
  Paso_Pattern_free(pattern);
if (MY_DEBUG) fprintf(stderr, "rank%d XCP6\n", rank);

  /* send/receive value array */
  j=0;
  for (p=0; p<recv->numNeighbors; p++) {
    i = ptr_ptr[recv->offsetInShared[p+1]] - ptr_ptr[recv->offsetInShared[p]];
    if (i > 0) 
	MPI_Irecv(&(A->remote_coupleBlock->val[j]), i * block_size, 
		MPI_DOUBLE, recv->neighbor[p],
                A->mpi_info->msg_tag_counter+recv->neighbor[p],
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p]));
    j += (i * block_size);
  }

  j=0;
  for (p=0; p<num_neighbors; p++) {
    i = send_degree[p] - j;
    if (i > 0)
	MPI_Issend(&(send_buf[j*block_size]), i*block_size, MPI_DOUBLE, send->neighbor[p],
                A->mpi_info->msg_tag_counter+rank,
                A->mpi_info->comm,
                &(A->row_coupler->mpi_requests[p+recv->numNeighbors]));
    j = send_degree[p];
  }

  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
                A->row_coupler->mpi_requests,
                A->row_coupler->mpi_stati);
  A->mpi_info->msg_tag_counter += mpi_size;

  /* release all temp memory allocation */
  TMPMEMFREE(cols);
  TMPMEMFREE(cols_array);
  TMPMEMFREE(recv_offset);
  TMPMEMFREE(recv_degree);
  TMPMEMFREE(recv_buf);
  TMPMEMFREE(send_buf);
  TMPMEMFREE(send_offset);
  TMPMEMFREE(send_degree);
  TMPMEMFREE(send_idx);

if (MY_DEBUG1) {
    char *str1, *str2;
    str1 = TMPMEMALLOC((overlapped_n + num_couple_cols)*15+50, char);
    str2 = TMPMEMALLOC(15, char);
    sprintf(str1, "Rank%d: Ptr=", rank);
    for (i=1; i<=overlapped_n; i++) {
	sprintf(str2, "(row%d)%d,", i-1, A->remote_coupleBlock->pattern->ptr[i]);
	strcat(str1, str2);
    }
    fprintf(stderr, "%s\n", str1); 

    for (i=1; i<=overlapped_n; i++) {
	sprintf(str1, "Rank%d: idx row %d=", rank, i-1);
	for (j= A->remote_coupleBlock->pattern->ptr[i-1]; 
		j< A->remote_coupleBlock->pattern->ptr[i]; j++) {
	  sprintf(str2, "%d,", A->remote_coupleBlock->pattern->index[j]);
	  strcat(str1, str2);
	}
	fprintf(stderr, "%s\n", str1);
	sprintf(str1, "Rank%d: val row %d=", rank, i-1);
        for (j= A->remote_coupleBlock->pattern->ptr[i-1];
                j< A->remote_coupleBlock->pattern->ptr[i]; j++) {
          sprintf(str2, "%g,", A->remote_coupleBlock->val[j]);
          strcat(str1, str2);
        }
        fprintf(stderr, "%s\n", str1);
    }
    TMPMEMFREE(str1);
    TMPMEMFREE(str2);
  }
}



