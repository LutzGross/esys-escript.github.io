//
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: defines AMG Restriction Operator  */

/**************************************************************/

/* Author: Lin Gao, lgao@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

/**************************************************************

    Methods nessecary for AMG preconditioner

    construct n_C x n the Restriction matrix R from A_p.
    
    R is the transpose of P. 
*/

#define MY_DEBUG 0
#define MY_DEBUG1 0
 
Paso_SystemMatrix* Paso_Preconditioner_AMG_getRestriction(Paso_SystemMatrix* P)
{ 
   Esys_MPIInfo *mpi_info=Esys_MPIInfo_getReference(P->mpi_info);
   Paso_SparseMatrix *main_block=NULL, *couple_block=NULL;
   Paso_SystemMatrix *out=NULL;
   Paso_SystemMatrixPattern *pattern=NULL;
   Paso_Distribution *input_dist=NULL, *output_dist=NULL;
   Paso_SharedComponents *send =NULL, *recv=NULL;
   Paso_Connector *col_connector=NULL, *row_connector=NULL;
   Paso_Coupler *coupler=NULL;
   Paso_Pattern *couple_pattern=NULL;
   const dim_t row_block_size=P->row_block_size;
   const dim_t col_block_size=P->col_block_size;
   const dim_t n=P->mainBlock->numRows;
   const dim_t n_C=P->mainBlock->numCols;
   const dim_t num_threads=omp_get_max_threads();
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

if (MY_DEBUG1)
fprintf(stderr, "Into Restriction %d %d\n", mpi_info->rank, rank);

   /* get main_block of R from the transpose of P->mainBlock */
   main_block = Paso_SparseMatrix_getTranspose(P->mainBlock);
if (MY_DEBUG1)
fprintf(stderr, "Rank %d CP1\n", rank);

   /* prepare "ptr" for the col_coupleBlock of R, start with get info about
      the degree_set (associated with "ptr"), offset_set (associated with
      "idx" and data_set (associated with "val") to be sent to other ranks */
   couple_block = P->col_coupleBlock;
   num_Pcouple_cols = couple_block->numCols;
   block_size = P->block_size;
   copy_block_size = block_size * sizeof(double);
   degree_set = TMPMEMALLOC(num_Pcouple_cols, index_t);
   send_ptr = TMPMEMALLOC(num_Pcouple_cols+1, index_t);
   memset(degree_set, 0, sizeof(index_t) * num_Pcouple_cols);
   for (i=0; i<n; i++) {
     iptr_ub = couple_block->pattern->ptr[i+1];
     for (iptr=couple_block->pattern->ptr[i]; iptr<iptr_ub; iptr++) {
	j = couple_block->pattern->index[iptr];
	degree_set[j] ++;
     }
   }

if (MY_DEBUG1)
fprintf(stderr, "rank %d Pcouple_cols %d block %d col_block %d\n", rank, num_Pcouple_cols, block_size, couple_block->block_size);
   send_ptr[0] = 0;
   for (i=0; i<num_Pcouple_cols; i++) {
     send_ptr[i+1] = send_ptr[i] + degree_set[i];
   }

   memset(degree_set, 0, sizeof(index_t) * num_Pcouple_cols);
   sum = couple_block->pattern->ptr[n];
   offset_set = TMPMEMALLOC(sum, index_t);
   data_set = TMPMEMALLOC(sum * block_size, double);
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
if (MY_DEBUG && rank == 3 && k <= 5) fprintf(stderr, "offset_set k %d i %d offset %d j %d send_ptr[j] %d degree_set[j] %d\n", k, i, offset, j, send_ptr[j], degree_set[j]);
	offset_set[k] = i + offset;   /* now we have the global id for row i,
					which will be used as col index of R */
//        memcpy(&(data_set[k*block_size]), &(couple_block->val[iptr*block_size]), copy_block_size);
	for (irb=0 ; irb < row_block_size; irb++) {
	  for (icb =0 ; icb < col_block_size; icb++) {
if (MY_DEBUG1 && iptr*block_size+irb+row_block_size*icb >= couple_block->len)
fprintf(stderr, "rank %d ACCESS OVERFLOW %d(iptr %d block %d irb %d row %d icb%d) >= %d (%d, %d)\n", rank, iptr*block_size+irb+row_block_size*icb, iptr, block_size, irb, row_block_size, icb, couple_block->len, couple_block->pattern->len, couple_block->pattern->ptr[n]); 
	    data_set[k*block_size+icb+col_block_size*irb] = couple_block->val[iptr*block_size+irb+row_block_size*icb];
	  }
	}
        degree_set[j] ++;
       }
     }
   }


   #ifdef ESYS_MPI
     mpi_requests=TMPMEMALLOC(size*4,MPI_Request);
     mpi_stati=TMPMEMALLOC(size*4,MPI_Status);
   #else
     mpi_requests=TMPMEMALLOC(size*4,int);
     mpi_stati=TMPMEMALLOC(size*4,int);
   #endif

   /* send/receive degree_set to build the "ptr" for R->col_coupleBlock */
   msgs = 0;
   send = P->col_coupler->connector->send;
if (MY_DEBUG) {
int my_sum = send->offsetInShared[1];
char *str1, *str2;
str1 = TMPMEMALLOC(my_sum*20+100, char);
str2 = TMPMEMALLOC(20, char);
sprintf(str1, "rank %d send->shared[%d] = \n", rank, my_sum);
for (i=0; i<my_sum; i++) {
  if (send->shared[i] < 0 || send->shared[i] >=n_C) {
  sprintf(str2, "(%d %d), ",i, send->shared[i]);
  strcat(str1, str2);
  }
}
fprintf(stderr, "%s\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
   recv = P->col_coupler->connector->recv;
   recv_ptr = TMPMEMALLOC(send->offsetInShared[send->numNeighbors], index_t);
   for (p=0; p<send->numNeighbors; p++) {
     i = send->offsetInShared[p];
     j = send->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
	MPI_Irecv(&(recv_ptr[i]), k, MPI_INT, send->neighbor[p],
		mpi_info->msg_tag_counter+send->neighbor[p],
		mpi_info->comm, &mpi_requests[msgs]); 
	msgs++;
     }
   }

   for (p=0; p<recv->numNeighbors; p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
if (MY_DEBUG && rank == 3 && recv->neighbor[p] == 4)
fprintf(stderr, "3 Send degree_set %d to 4 (%d) offset %d p %d\n", k, degree_set[i], i, p);
	MPI_Issend(&(degree_set[i]), k, MPI_INT, recv->neighbor[p],
		mpi_info->msg_tag_counter+rank, mpi_info->comm, 
                &mpi_requests[msgs]);
	msgs++;
     }
   }

   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   mpi_info->msg_tag_counter += size;
if (MY_DEBUG1)
fprintf(stderr, "rank %d Waitall(1)\n", rank);

   TMPMEMFREE(degree_set);
   degree_set = TMPMEMALLOC(send->numNeighbors, index_t);
   memset(degree_set, 0, sizeof(index_t)*send->numNeighbors);
   for (p=0, sum=0; p<send->numNeighbors; p++) {
     iptr_ub = send->offsetInShared[p+1];
     for (iptr = send->offsetInShared[p]; iptr < iptr_ub; iptr++){
	degree_set[p] += recv_ptr[iptr];
     }
     sum += degree_set[p];
   }
if (MY_DEBUG && rank == 4)
fprintf(stderr, "rank %d degree_set %d (%d %d)\n", rank, sum, degree_set[0], degree_set[1]);

   /* send/receive offset_set and data_set to build the "idx" and "val"
      for R->col_coupleBlock */
   msgs = 0;
   recv_idx = TMPMEMALLOC(sum, index_t);
   recv_val = TMPMEMALLOC(sum * block_size, double);
   for (p=0, offset=0; p<send->numNeighbors; p++) {
     if (degree_set[p]) {
if (MY_DEBUG && rank == 4)
fprintf(stderr, "4 Recv %d from %d\n", degree_set[p], send->neighbor[p]);
	MPI_Irecv(&(recv_idx[offset]), degree_set[p], MPI_INT,
		send->neighbor[p], mpi_info->msg_tag_counter+send->neighbor[p],
		mpi_info->comm, &mpi_requests[msgs]);
	msgs++;
	MPI_Irecv(&(recv_val[offset*block_size]), degree_set[p] * block_size, 
		MPI_DOUBLE, send->neighbor[p], 
		mpi_info->msg_tag_counter+send->neighbor[p]+size,
		mpi_info->comm, &mpi_requests[msgs]);
	offset += degree_set[p];
	msgs++;
     }
   }

   for (p=0; p<recv->numNeighbors; p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = send_ptr[j] - send_ptr[i];
     if (k > 0) {
if (MY_DEBUG && rank == 3 && recv->neighbor[p] == 4)
fprintf(stderr, "3 Send %d to %d (%d %d) offset %d\n", k, recv->neighbor[p], offset_set[i], offset_set[i+1], i);
	MPI_Issend(&(offset_set[send_ptr[i]]), k, MPI_INT,
		recv->neighbor[p], mpi_info->msg_tag_counter+rank,
		mpi_info->comm, &mpi_requests[msgs]);
	msgs++;
	MPI_Issend(&(data_set[send_ptr[i]*block_size]), k*block_size, MPI_DOUBLE,
		recv->neighbor[p], mpi_info->msg_tag_counter+rank+size,
		mpi_info->comm, &mpi_requests[msgs]);
	msgs++;
     }
   }

   len = send->offsetInShared[send->numNeighbors];
   temp = TMPMEMALLOC(len, index_t);
   memset(temp, 0, sizeof(index_t)*len);
   for (p=1; p<len; p++) {
     temp[p] = temp[p-1] + recv_ptr[p-1];
   }

   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   mpi_info->msg_tag_counter += 2*size;
   TMPMEMFREE(degree_set);
   TMPMEMFREE(offset_set);
   TMPMEMFREE(data_set);
   TMPMEMFREE(send_ptr);
   TMPMEMFREE(mpi_requests);
   TMPMEMFREE(mpi_stati);
if (MY_DEBUG && rank == 4){
char *str1, *str2;
int my_sum, my_i;
my_sum = 2;
str1 = TMPMEMALLOC(my_sum*100+100, char);
str2 = TMPMEMALLOC(100, char);
sprintf(str1, "recv_idx[%d] = (", my_sum);
for (my_i=0; my_i<my_sum; my_i++) {
  sprintf(str2, "%d ", recv_idx[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}


if (MY_DEBUG1)
fprintf(stderr, "rank %d Construct col_coupleBlock for R: %d %d \n", rank, n_C, sum);

   /* construct "ptr", "idx" and "val" for R->col_coupleBlock */
   ptr = MEMALLOC(n_C + 1, index_t);
   idx = MEMALLOC(sum, index_t);
   val = MEMALLOC(sum*block_size, double);
   ptr[0] = 0;
   for (i=0; i<n_C; i++) {
     icb = 0;
     for (p=0; p<send->numNeighbors; p++) {
        k = send->offsetInShared[p+1];
        for (j=send->offsetInShared[p]; j<k; j++) {
if (MY_DEBUG)
fprintf(stderr, "rank %d send %d, i%d len%d\n", rank, send->shared[j], i, recv_ptr[j]);
//if (send->shared[j] < 0 || send->shared[j] >= n_C) 
//fprintf(stderr, "rank %d shared ele out of range %d %d %d\n", rank, j, send->shared[j], n_C);
          if (send->shared[j] == i) {
	    offset = ptr[i] + icb;
	    len = recv_ptr[j];
	    memcpy(&idx[offset], &recv_idx[temp[j]], sizeof(index_t)*len);
if (temp[j] + len > sum) 
fprintf(stderr, "rank %d ACCESS to recv_idx overflowed i%d j%d offset%d t_j%d len %d\n", rank, i, j, offset, temp[j], len);
	    memcpy(&val[offset*block_size], &recv_val[temp[j]*block_size], sizeof(double)*len*block_size); 
	    icb += len;
if (MY_DEBUG)
fprintf(stderr, "rank %d send %d, i%d len%d\n", rank, send->shared[j], i, recv_ptr[j]);
if (MY_DEBUG) fprintf(stderr, "rank %d offset %d len %d block %d\n", rank, offset, len, block_size);
            break;
          }
        }
     }
     ptr[i+1] = ptr[i] + icb;
   }
   sum = ptr[n_C];
   TMPMEMFREE(temp);
   TMPMEMFREE(recv_ptr);
   TMPMEMFREE(recv_val);

if (MY_DEBUG1){
fprintf(stderr, "rank %d col_coupleBlock non-zero: %d\n", rank, sum);
}

   /* count the number of cols (num_Rcouple_cols) in R->col_coupleBlock, 
      and convert the global id in "idx" into local id */
   num_Rcouple_cols = 0;
   if (sum) {
     #ifdef USE_QSORTG
       qsortG(recv_idx, (size_t)sum, sizeof(index_t), Paso_comparIndex);
     #else
       qsort(recv_idx, (size_t)sum, sizeof(index_t), Paso_comparIndex);
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
     for (i=0; i<sum; i++) {
	where_p = (index_t *)bsearch(&(idx[i]), recv_idx, num_Rcouple_cols,
				sizeof(index_t), Paso_comparIndex);
if (where_p == NULL)
fprintf(stderr, "rank %d index out of range, idx[%d]=%d\n", rank, i, idx[i]);
	idx[i] = (index_t)(where_p - recv_idx);
if (idx[i] >= num_Rcouple_cols) 
fprintf(stderr, "rank %d index out of range (%d %d %d)\n", rank, i, idx[i], num_Rcouple_cols);
     }
   }

if (MY_DEBUG1)
fprintf(stderr, "rank %d Count num_Rcouple_cols %d\n", rank, num_Rcouple_cols);

   /* prepare the receiver for the col_connector */
   dist = P->pattern->output_distribution->first_component;
   offsetInShared = TMPMEMALLOC(size+1, index_t);
   shared = TMPMEMALLOC(num_Rcouple_cols, index_t);
   numNeighbors = send->numNeighbors;
   neighbor = send->neighbor;
   memset(offsetInShared, 0, sizeof(index_t) * (size+1));
if (MY_DEBUG && rank == 0){
char *str1, *str2;
int sum, my_i;
sum = num_Rcouple_cols;
str1 = TMPMEMALLOC((sum+size)*100+100, char);
str2 = TMPMEMALLOC(100, char);
sprintf(str1, "rank %d distribution[%d] = (", rank, size);
for (my_i=0; my_i<size; my_i++) {
  sprintf(str2, "%d ", dist[my_i+1]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
sprintf(str1, "rank %d recv_idx[%d] = (", rank, sum);
for (my_i=0; my_i<sum; my_i++) {
  sprintf(str2, "%d ", recv_idx[my_i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
if (MY_DEBUG1) fprintf(stderr, "rank %d CXXP1\n", rank);
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
if (MY_DEBUG1) fprintf(stderr, "rank %d CXXP2\n", rank);
   for (i=p; i<numNeighbors; i++) {
     offsetInShared[i+1] = num_Rcouple_cols;
   }
//   offsetInShared[numNeighbors] = num_Rcouple_cols;
   recv = Paso_SharedComponents_alloc(n, numNeighbors,
		neighbor, shared, offsetInShared, 1, 0, mpi_info);
   TMPMEMFREE(recv_idx);

if (MY_DEBUG1)
fprintf(stderr, "rank %d Receiver!! %d\n", rank, numNeighbors);

   /* prepare the sender for the col_connector */
   TMPMEMFREE(shared);
   numNeighbors = P->col_coupler->connector->recv->numNeighbors;
   neighbor = P->col_coupler->connector->recv->neighbor;
   shared = TMPMEMALLOC(n * numNeighbors, index_t);
   couple_pattern = P->col_coupleBlock->pattern;
   sum=0;
   memset(offsetInShared, 0, sizeof(index_t) * (size+1));
   for (p=0; p<numNeighbors; p++) {
     j = P->col_coupler->connector->recv->offsetInShared[p];
     j_ub = P->col_coupler->connector->recv->offsetInShared[p+1];
if (MY_DEBUG && rank ==3 && P->col_coupler->connector->recv->neighbor[p] == 4) {
fprintf(stderr, "3sendto4 with P=%d\n", p);
}
     for (i=0; i<n; i++) {
	iptr = couple_pattern->ptr[i];
	iptr_ub = couple_pattern->ptr[i+1];
	for (; iptr<iptr_ub; iptr++) {
	  k = couple_pattern->index[iptr];
	  if (k >= j && k < j_ub) {
if (MY_DEBUG && rank ==3 && p == 1 && P->col_coupler->connector->recv->neighbor[p] == 2) {
fprintf(stderr, "3 send to 1: [%d]=k%d i%d n %d (%d, %d)\n", sum, k, i, n, j, j_ub);
}
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
if (MY_DEBUG1)
fprintf(stderr, "rank %d Sender!! %d %d\n", rank, n_C, num_Rcouple_cols);

   /* build the col_connector based on sender and receiver */
   col_connector = Paso_Connector_alloc(send, recv);
   Paso_SharedComponents_free(recv);
   Paso_SharedComponents_free(send);
   TMPMEMFREE(offsetInShared);
   TMPMEMFREE(shared);   

   couple_pattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, n_C,
                        num_Rcouple_cols, ptr, idx);

   input_dist = Paso_Distribution_alloc(mpi_info, dist, 1, 0);
   dist = P->pattern->input_distribution->first_component;
   output_dist = Paso_Distribution_alloc(mpi_info, dist, 1, 0);

if (MY_DEBUG)
fprintf(stderr, "rank %d Before alloc System Pattern for Restriction %d %s\n", rank, Esys_getErrorType(), Esys_getErrorMessage());

   /* now we need to create the System Matrix 
      TO BE FIXED: at this stage, we only construction col_couple_pattern
      and col_connector for Restriction matrix R. To be completed, 
      row_couple_pattern and row_connector need to be constructed as well */
   if (Esys_noError()) {
     pattern = Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT, 
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
   MEMFREE(val);

if (MY_DEBUG)
fprintf(stderr, "rank %d After alloc System Pattern for Restriction\n", rank);


   /* clean up */ 
   Paso_SparseMatrix_free(main_block);
   Paso_SystemMatrixPattern_free(pattern);
   Paso_Pattern_free(couple_pattern);
   Paso_Connector_free(col_connector);
   Paso_Distribution_free(output_dist);
   Paso_Distribution_free(input_dist);

   if (Esys_noError()) {
      return out;
   } else {
      Paso_SystemMatrix_free(out);
      return NULL;
   }
}

