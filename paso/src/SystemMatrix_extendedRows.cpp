
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


/**********************************************************************************/
/* Paso: SystemMatrix                                       */
/*							    */
/*  Extend the ST sets of rows in row_coupleBlock 	    */
/*  Input: SystemMatrix A and ST info			    */
/*  Output: 						    */
/*	degree_ST:					    */
/*	offset_ST:					    */
/*	ST: 					 	    */
/**********************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/**********************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void Paso_SystemMatrix_extendedRowsForST(Paso_SystemMatrix* A, dim_t* degree_ST, index_t* offset_ST, index_t* ST)
{
    paso::Coupler_ptr coupler;
    double *cols=NULL, *rows=NULL;
    index_t *global_id=NULL, *B=NULL, *send_buf=NULL;
    index_t *send_ST=NULL, *recv_ST=NULL, *recv_offset_ST=NULL;
    index_t i, j, k, p, z, z0, z1, size, rank, offset, my_n, len, overlapped_n;
    index_t num_main_cols, num_couple_cols;
    bool flag;
    dim_t *recv_degree_ST=NULL;

    if (A->mpi_info->size == 1) return;

  /* sending/receiving unknown's global ID */
  num_main_cols = A->mainBlock->numCols;
  cols = new double[num_main_cols];
  rank = A->mpi_info->rank;
  offset = A->col_distribution->first_component[rank];
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<num_main_cols; ++i) cols[i] = offset + i;
  if (A->global_id == NULL) {
    coupler.reset(new paso::Coupler(A->col_coupler->connector, 1));
    coupler->startCollect(cols);
  }

  my_n = A->mainBlock->numRows;
  rows = new double[my_n];
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<my_n; i++) rows[i] = degree_ST[i];

  num_couple_cols = A->col_coupleBlock->numCols;
  size = num_main_cols + num_couple_cols;
  overlapped_n = A->row_coupleBlock->numRows;
  recv_offset_ST = new index_t[overlapped_n+1];
  recv_degree_ST = new dim_t[overlapped_n];
  send_ST = new index_t[offset_ST[my_n]];
  len = A->row_coupler->connector->send->offsetInShared[A->row_coupler->connector->send->numNeighbors] * size;
  send_buf = new index_t[len];


  /* waiting for receiving unknown's global ID */
  if (A->global_id == NULL) {
      coupler->finishCollect();
    global_id = new index_t[num_couple_cols];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_couple_cols; ++i) 
	global_id[i] = coupler->recv_buffer[i];
    A->global_id = global_id;
  }
  global_id = A->global_id;

  /* sending/receiving the degree_ST */
  coupler.reset(new paso::Coupler(A->row_coupler->connector, 1));
  coupler->startCollect(rows);

  /* prepare ST with global ID */
  B = new index_t[size*2];
  /* find the point z in array of global_id that 
     forall i < z, global_id[i] < offset; and
     forall i >= z, global_id[i] > offset + my_n */
  for (i=0; i<num_couple_cols; i++) 
    if (global_id[i] > offset) break;
  z = i;
  for (i=0; i<num_main_cols; i++) {
    p = offset_ST[i+1];
    z0 = 0;
    z1 = offset_ST[i];
    if (degree_ST[i] > 0 && (ST[p-1] < my_n || z == 0)) flag = 1;
    else flag = 0;
    for (j=offset_ST[i]; j<p; j++) {
	k = ST[j];
	if (flag == 0 && k < my_n) {
	  send_buf[z0] = k + offset;
	  z0++;
	} else if (flag == 0 && k - my_n >= z) {
	  if (z0 >0) 
	    memcpy(&(send_ST[z1]), &(send_buf[0]), z0*sizeof(index_t));
	  z1 += z0;
	  flag = 1;
	  send_ST[z1] = global_id[k - my_n];
	  z1++;
	} else if (flag == 0 && j == p-1) {
	  send_ST[z1] = global_id[k - my_n];
          z1++;
	  if (z0 >0)
	    memcpy(&(send_ST[z1]), &(send_buf[0]), z0*sizeof(index_t));
	  flag = 1;
          z1 += z0;
	} else if (k < my_n) {
	  send_ST[z1] = k + offset;
	  z1++;
	} else {
	  send_ST[z1] = global_id[k - my_n];
	  z1++;
	}
    }
  } 

  /* waiting for receiving the degree_ST */
  coupler->finishCollect();
  delete[] rows;

  /* preparing degree_ST and offset_ST for the to-be-received extended rows */
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<overlapped_n; i++) recv_degree_ST[i] = coupler->recv_buffer[i];
  recv_offset_ST[0] = 0;
  for (i=0; i<overlapped_n; i++) {
    recv_offset_ST[i+1] = recv_offset_ST[i] + coupler->recv_buffer[i];
  }
  recv_ST = new index_t[recv_offset_ST[overlapped_n]];
  coupler.reset();

  /* receiving ST for the extended rows */
  z = 0;
  for (p=0; p<A->row_coupler->connector->recv->numNeighbors; p++) {
    const index_t j_min = A->row_coupler->connector->recv->offsetInShared[p];
    const index_t j_max = A->row_coupler->connector->recv->offsetInShared[p+1];
    j = recv_offset_ST[j_max] - recv_offset_ST[j_min];
    #ifdef ESYS_MPI
    MPI_Irecv(&(recv_ST[z]), j,  MPI_INT, 
		A->row_coupler->connector->recv->neighbor[p],
		A->mpi_info->msg_tag_counter+A->row_coupler->connector->recv->neighbor[p],
		A->mpi_info->comm,
		&(A->row_coupler->mpi_requests[p]) );
    #endif
    z += j;
  }

  /* sending ST for the extended rows */
  z0 = 0;
  for (p=0; p<A->row_coupler->connector->send->numNeighbors; p++) {
    const index_t j_min = A->row_coupler->connector->send->offsetInShared[p];
    const index_t j_max = A->row_coupler->connector->send->offsetInShared[p+1];
    z = z0;
    for (j=j_min; j<j_max; j++) {
	const index_t row=A->row_coupler->connector->send->shared[j];
	if (degree_ST[row] > 0) {
	  memcpy(&(send_buf[z]), &(send_ST[offset_ST[row]]), degree_ST[row] * sizeof(index_t));
	  z += degree_ST[row];
	}
    }
    #ifdef ESYS_MPI
    MPI_Issend(&(send_buf[z0]),z-z0, MPI_INT, 
		 A->row_coupler->connector->send->neighbor[p],
		 A->mpi_info->msg_tag_counter+A->mpi_info->rank,
		 A->mpi_info->comm,
		 &(A->row_coupler->mpi_requests[p+A->row_coupler->connector->recv->numNeighbors]));
    #endif
    z0 = z;
  }

  /* first merge "cols" and "global_id" into array "B" */
  i = 0;
  j = 0;
  z = 0;
  while (i+j < size) {
    if (i < num_main_cols && j < num_couple_cols) {
        if (cols[i] < global_id[j]) {
          B[2*z] = cols[i];
	  B[2*z+1] = i;
          i++;
        } else {
          B[2*z] = global_id[j];
	  B[2*z+1] = j + num_main_cols;
          j++;
        }
        z++;
    } else if (i >= num_main_cols) {
        for (; j<num_couple_cols; j++, z++) {
	  B[2*z] = global_id[j];
	  B[2*z+1] = j + num_main_cols;
	}
    } else {/* j >= num_couple_cols */
        for (; i<num_main_cols; i++, z++) {
	  B[2*z] = cols[i];
	  B[2*z+1] = i;
	}
    }
  }

  /* wait til everything is done */
  #ifdef ESYS_MPI
  MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
	      A->row_coupler->mpi_requests,
	      A->row_coupler->mpi_stati);
  #endif
  ESYS_MPI_INC_COUNTER(*(A->mpi_info), A->mpi_info->size)

  /* filter the received ST (for extended rows) with cols in mainBlock as
     well as cols in col_coupleBlock, their global ids are listed in "B" */
  len = offset_ST[my_n];
  size = 2 * size;
  for (i=0; i<overlapped_n; i++) {
    p = recv_offset_ST[i+1];
    z = 0;
    for (j=recv_offset_ST[i]; j<p; j++) {
	if (recv_ST[j] == B[z]) {
	  ST[len] = B[z+1];
	  len++;
	  z+=2;
	  if (z >= size) break;
	} else if (recv_ST[j] > B[z]) {
	  while (recv_ST[j] > B[z] && z < size) z+=2;
	  if (z >= size) break;
	  if (recv_ST[j] == B[z]) {
	    ST[len] = B[z+1];
	    len++;
	    z+=2;
	    if (z >= size) break;
	  }
	}
    }
    j = my_n + i;
    degree_ST[j] = len - offset_ST[j];
    offset_ST[j+1] = len;
    qsort(&(ST[offset_ST[j]]), (size_t) degree_ST[j], sizeof(index_t), paso::comparIndex);
  }

  /* release all temp memory allocation */
  delete[] cols;
  delete[] B;
  delete[] send_buf;
  delete[] send_ST;
  delete[] recv_ST;
  delete[] recv_offset_ST;
  delete[] recv_degree_ST;
}


