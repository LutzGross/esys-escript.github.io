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


/************************************************************************************/

/* Paso: defines AMG Restriction Operator  */

/************************************************************************************/

/* Author: Lin Gao, lgao@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

/************************************************************************************

    Methods necessary for AMG preconditioner

    construct n_C x n the Restriction matrix R from A_p.
    
    R->mainBlock is the transpose of P->mainBlock, but we need
    to recover R->col_coupleBlock from P's data in other ranks.
 
*************************************************************************************/

Paso_SystemMatrix* Paso_Preconditioner_AMG_getRestriction(Paso_SystemMatrix* P)
{ 
   Esys_MPIInfo *mpi_info=Esys_MPIInfo_getReference(P->mpi_info);
   paso::SparseMatrix *main_block=NULL, *couple_block=NULL;
   Paso_SystemMatrix *out=NULL;
   paso::SystemMatrixPattern *pattern=NULL;
   Paso_Distribution *input_dist=NULL, *output_dist=NULL;
   Paso_SharedComponents *send =NULL, *recv=NULL;
   paso::Connector *col_connector=NULL;
   paso::Pattern *couple_pattern=NULL;
   const dim_t row_block_size=P->row_block_size;
   const dim_t col_block_size=P->col_block_size;
   const dim_t n=P->mainBlock->numRows;
   const dim_t n_C=P->mainBlock->numCols;
   index_t size=mpi_info->size, rank=mpi_info->rank, *dist=NULL;
   index_t *ptr=NULL, *idx=NULL, *degree_set=NULL, *offset_set=NULL;
   index_t *send_ptr=NULL, *recv_ptr=NULL, *recv_idx=NULL;
   index_t *temp=NULL, *where_p=NULL;
   index_t num_Pcouple_cols, num_Rcouple_cols, numNeighbors;
   index_t i, j, j_ub, k, p, iptr, iptr_ub, icb, irb;
   index_t block_size, copy_block_size, sum, offset, len, msgs;
   double  *val=NULL, *data_set=NULL, *recv_val=NULL;
   index_t *shared=NULL, *offsetInShared=NULL;
   Esys_MPI_rank *neighbor=NULL;
   #ifdef ESYS_MPI
     MPI_Request* mpi_requests=NULL;
     MPI_Status* mpi_stati=NULL;
   #else
     int *mpi_requests=NULL, *mpi_stati=NULL;
   #endif

   /* get main_block of R from the transpose of P->mainBlock */
   main_block = paso::SparseMatrix_getTranspose(P->mainBlock);

   /* prepare "ptr" for the col_coupleBlock of R, start with get info about
      the degree_set (associated with "ptr"), offset_set (associated with
      "idx" and data_set (associated with "val") to be sent to other ranks */
   couple_block = P->col_coupleBlock;
   num_Pcouple_cols = couple_block->numCols;
   block_size = P->block_size;
   copy_block_size = block_size * sizeof(double);
   degree_set = new index_t[num_Pcouple_cols];
   send_ptr = new index_t[num_Pcouple_cols+1];
   memset(degree_set, 0, sizeof(index_t) * num_Pcouple_cols);
   for (i=0; i<n; i++) {
     iptr_ub = couple_block->pattern->ptr[i+1];
     for (iptr=couple_block->pattern->ptr[i]; iptr<iptr_ub; iptr++) {
	j = couple_block->pattern->index[iptr];
	degree_set[j] ++;
     }
   }

   send_ptr[0] = 0;
   for (i=0; i<num_Pcouple_cols; i++) {
     send_ptr[i+1] = send_ptr[i] + degree_set[i];
   }

   memset(degree_set, 0, sizeof(index_t) * num_Pcouple_cols);
   sum = couple_block->pattern->ptr[n];
   offset_set = new index_t[sum];
   data_set = new double[sum * block_size];
   offset = P->pattern->output_distribution->first_component[rank];

   if (P->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
     for (i=0; i<n; i++) {
       iptr_ub = couple_block->pattern->ptr[i+1];
       for (iptr=couple_block->pattern->ptr[i]; iptr<iptr_ub; iptr++) {
	j = couple_block->pattern->index[iptr];
	k = send_ptr[j] + degree_set[j];
	offset_set[k] = i + offset;   /* now we have the global id for row i,
                                        which will be used as col index of R */
	memcpy(&(data_set[k*block_size]), &(couple_block->val[iptr*block_size]), copy_block_size);
	degree_set[j] ++;
       }
     }
   } else {
     for (i=0; i<n; i++) {
       iptr_ub = couple_block->pattern->ptr[i+1];
       for (iptr=couple_block->pattern->ptr[i]; iptr<iptr_ub; iptr++) {
        j = couple_block->pattern->index[iptr];
	k = send_ptr[j] + degree_set[j];
	offset_set[k] = i + offset;   /* now we have the global id for row i,
					which will be used as col index of R */
	for (irb=0 ; irb < row_block_size; irb++) {
	  for (icb =0 ; icb < col_block_size; icb++) {
	    data_set[k*block_size+icb+col_block_size*irb] = couple_block->val[iptr*block_size+irb+row_block_size*icb];
	  }
	}
        degree_set[j] ++;
       }
     }
   }


   #ifdef ESYS_MPI
     mpi_requests=new MPI_Request[size*4];
     mpi_stati=new MPI_Status[size*4];
   #else
     mpi_requests=new int[size*4];
     mpi_stati=new int[size*4];
   #endif

   /* send/receive degree_set to build the "ptr" for R->col_coupleBlock */
   msgs = 0;
   send = P->col_coupler->connector->send;
   recv = P->col_coupler->connector->recv;
   recv_ptr = new index_t[send->offsetInShared[send->numNeighbors]];
   for (p=0; p<send->numNeighbors; p++) {
     i = send->offsetInShared[p];
     j = send->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
	#ifdef ESYS_MPI
	MPI_Irecv(&(recv_ptr[i]), k, MPI_INT, send->neighbor[p],
		mpi_info->msg_tag_counter+send->neighbor[p],
		mpi_info->comm, &mpi_requests[msgs]); 
	#endif
	msgs++;
     }
   }

   for (p=0; p<recv->numNeighbors; p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
	#ifdef ESYS_MPI
	MPI_Issend(&(degree_set[i]), k, MPI_INT, recv->neighbor[p],
		mpi_info->msg_tag_counter+rank, mpi_info->comm, 
                &mpi_requests[msgs]);
	#endif
	msgs++;
     }
   }

   #ifdef ESYS_MPI
   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   #endif
   ESYS_MPI_INC_COUNTER(*mpi_info, size)

   delete[] degree_set;
   degree_set = new index_t[send->numNeighbors];
   memset(degree_set, 0, sizeof(index_t)*send->numNeighbors);
   for (p=0, sum=0; p<send->numNeighbors; p++) {
     iptr_ub = send->offsetInShared[p+1];
     for (iptr = send->offsetInShared[p]; iptr < iptr_ub; iptr++){
	degree_set[p] += recv_ptr[iptr];
     }
     sum += degree_set[p];
   }

   /* send/receive offset_set and data_set to build the "idx" and "val"
      for R->col_coupleBlock */
   msgs = 0;
   recv_idx = new index_t[sum];
   recv_val = new double[sum * block_size];
   for (p=0, offset=0; p<send->numNeighbors; p++) {
     if (degree_set[p]) {
	#ifdef ESYS_MPI
	MPI_Irecv(&(recv_idx[offset]), degree_set[p], MPI_INT,
		send->neighbor[p], mpi_info->msg_tag_counter+send->neighbor[p],
		mpi_info->comm, &mpi_requests[msgs]);
	msgs++;
	MPI_Irecv(&(recv_val[offset*block_size]), degree_set[p] * block_size, 
		MPI_DOUBLE, send->neighbor[p], 
		mpi_info->msg_tag_counter+send->neighbor[p]+size,
		mpi_info->comm, &mpi_requests[msgs]);
	offset += degree_set[p];
	#endif
	msgs++;
     }
   }

   for (p=0; p<recv->numNeighbors; p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = send_ptr[j] - send_ptr[i];
     if (k > 0) {
	#ifdef ESYS_MPI
	MPI_Issend(&(offset_set[send_ptr[i]]), k, MPI_INT,
		recv->neighbor[p], mpi_info->msg_tag_counter+rank,
		mpi_info->comm, &mpi_requests[msgs]);
	msgs++;
	MPI_Issend(&(data_set[send_ptr[i]*block_size]), k*block_size, MPI_DOUBLE,
		recv->neighbor[p], mpi_info->msg_tag_counter+rank+size,
		mpi_info->comm, &mpi_requests[msgs]);
	#endif
	msgs++;
     }
   }

   len = send->offsetInShared[send->numNeighbors];
   temp = new index_t[len];
   memset(temp, 0, sizeof(index_t)*len);
   for (p=1; p<len; p++) {
     temp[p] = temp[p-1] + recv_ptr[p-1];
   }

   #ifdef ESYS_MPI
   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   #endif
   ESYS_MPI_INC_COUNTER(*mpi_info, 2*size)
   delete[] degree_set;
   delete[] offset_set;
   delete[] data_set;
   delete[] send_ptr;
   delete[] mpi_requests;
   delete[] mpi_stati;

   /* construct "ptr", "idx" and "val" for R->col_coupleBlock */
   ptr = new  index_t[n_C + 1];
   idx = new  index_t[sum];
   val = new  double[sum*block_size];
   ptr[0] = 0;
   for (i=0; i<n_C; i++) {
     icb = 0;
     for (p=0; p<send->numNeighbors; p++) {
        k = send->offsetInShared[p+1];
        for (j=send->offsetInShared[p]; j<k; j++) {
          if (send->shared[j] == i) {
	    offset = ptr[i] + icb;
	    len = recv_ptr[j];
	    memcpy(&idx[offset], &recv_idx[temp[j]], sizeof(index_t)*len);
	    memcpy(&val[offset*block_size], &recv_val[temp[j]*block_size], sizeof(double)*len*block_size); 
	    icb += len;
            break;
          }
        }
     }
     ptr[i+1] = ptr[i] + icb;
   }
   sum = ptr[n_C];
   delete[] temp;
   delete[] recv_ptr;
   delete[] recv_val;

   /* count the number of cols (num_Rcouple_cols) in R->col_coupleBlock, 
      and convert the global id in "idx" into local id */
   num_Rcouple_cols = 0;
   if (sum) {
     #ifdef USE_QSORTG
       qsortG(recv_idx, (size_t)sum, sizeof(index_t), paso::comparIndex);
     #else
       qsort(recv_idx, (size_t)sum, sizeof(index_t), paso::comparIndex);
     #endif
     num_Rcouple_cols = 1;
     i = recv_idx[0];
     for (j=1; j<sum; j++) {
	if (recv_idx[j] > i) {
	  i = recv_idx[j];
	  recv_idx[num_Rcouple_cols] = i;
	  num_Rcouple_cols++;
	}
     }
     #pragma omp parallel for private(i,where_p) schedule(static)
     for (i=0; i<sum; i++) {
	where_p = (index_t *)bsearch(&(idx[i]), recv_idx, num_Rcouple_cols,
				sizeof(index_t), paso::comparIndex);
	idx[i] = (index_t)(where_p - recv_idx);
     }
   }

   /* prepare the receiver for the col_connector */
   dist = P->pattern->output_distribution->first_component;
   offsetInShared = new index_t[size+1];
   shared = new index_t[num_Rcouple_cols];
   numNeighbors = send->numNeighbors;
   neighbor = send->neighbor;
   memset(offsetInShared, 0, sizeof(index_t) * (size+1));
   if (num_Rcouple_cols > 0) offset = dist[neighbor[0] + 1];
   for (i=0, p=0; i<num_Rcouple_cols; i++) {
     /* cols i is received from rank neighbor[p] when it's still smaller
	than "offset", otherwise, it is received from rank neighbor[p+1] */
     while (recv_idx[i] >= offset) { 
	p++;
	offsetInShared[p] = i;
	offset = dist[neighbor[p] + 1];
     }
     shared[i] = i + n;  /* n is the number of cols in R->mainBlock */
   }
   #pragma omp parallel for private(i) schedule(static)
   for (i=p; i<numNeighbors; i++) {
     offsetInShared[i+1] = num_Rcouple_cols;
   }
   recv = Paso_SharedComponents_alloc(n, numNeighbors,
		neighbor, shared, offsetInShared, 1, 0, mpi_info);
   delete[] recv_idx;

   /* prepare the sender for the col_connector */
   delete[] shared;
   numNeighbors = P->col_coupler->connector->recv->numNeighbors;
   neighbor = P->col_coupler->connector->recv->neighbor;
   shared = new index_t[n * numNeighbors];
   couple_pattern = P->col_coupleBlock->pattern;
   sum=0;
   memset(offsetInShared, 0, sizeof(index_t) * (size+1));
   for (p=0; p<numNeighbors; p++) {
     j = P->col_coupler->connector->recv->offsetInShared[p];
     j_ub = P->col_coupler->connector->recv->offsetInShared[p+1];
     for (i=0; i<n; i++) {
	iptr = couple_pattern->ptr[i];
	iptr_ub = couple_pattern->ptr[i+1];
	for (; iptr<iptr_ub; iptr++) {
	  k = couple_pattern->index[iptr];
	  if (k >= j && k < j_ub) {
	    shared[sum] = i;
	    sum++;
	    break;
	  }
	}
     }
     offsetInShared[p+1] = sum;
   }
   send = Paso_SharedComponents_alloc(n, numNeighbors,
                neighbor, shared, offsetInShared, 1, 0, mpi_info);

   /* build the col_connector based on sender and receiver */
   col_connector = paso::Connector_alloc(send, recv);
   Paso_SharedComponents_free(recv);
   Paso_SharedComponents_free(send);
   delete[] offsetInShared;
   delete[] shared;   

   couple_pattern = paso::Pattern_alloc(MATRIX_FORMAT_DEFAULT, n_C,
                        num_Rcouple_cols, ptr, idx);

   input_dist = Paso_Distribution_alloc(mpi_info, dist, 1, 0);
   dist = P->pattern->input_distribution->first_component;
   output_dist = Paso_Distribution_alloc(mpi_info, dist, 1, 0);

   /* now we need to create the System Matrix 
      TO BE FIXED: at this stage, we only construction col_couple_pattern
      and col_connector for Restriction matrix R. To be completed, 
      row_couple_pattern and row_connector need to be constructed as well */
   if (Esys_noError()) {
     pattern = new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT, 
		output_dist, input_dist, main_block->pattern, couple_pattern, 
		couple_pattern, col_connector, col_connector);
     out = Paso_SystemMatrix_alloc(MATRIX_FORMAT_DIAGONAL_BLOCK, pattern,
		row_block_size, col_block_size, FALSE);
   } 

   /* now fill in the matrix */
   memcpy(out->mainBlock->val, main_block->val, 
		main_block->len * sizeof(double));
   memcpy(out->col_coupleBlock->val, val,
		out->col_coupleBlock->len * sizeof(double));
   delete[] val;

   /* clean up */ 
   paso::SparseMatrix_free(main_block);
   paso::SystemMatrixPattern_free(pattern);
   paso::Pattern_free(couple_pattern);
   paso::Connector_free(col_connector);
   Paso_Distribution_free(output_dist);
   Paso_Distribution_free(input_dist);

   if (Esys_noError()) {
      return out;
   } else {
      Paso_SystemMatrix_free(out);
      return NULL;
   }
}

