/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: defines AMG Restriction Operator  */

/****************************************************************************/

/* Author: Lin Gao, lgao@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "SparseMatrix.h"

#include <cstring> // memcpy

namespace paso {

/****************************************************************************

    Methods necessary for AMG preconditioner

    construct n_C x n the Restriction matrix R from A_p.

    R->mainBlock is the transpose of P->mainBlock, but we need
    to recover R->col_coupleBlock from P's data in other ranks.

*****************************************************************************/

SystemMatrix_ptr Preconditioner_AMG_getRestriction(SystemMatrix_ptr P)
{
   escript::JMPI mpi_info(P->mpi_info);
   escript::Distribution_ptr input_dist, output_dist;
   Connector_ptr col_connector;
   const dim_t row_block_size=P->row_block_size;
   const dim_t col_block_size=P->col_block_size;
   const dim_t n=P->mainBlock->numRows;
   const dim_t n_C=P->mainBlock->numCols;
   index_t size=mpi_info->size, rank=mpi_info->rank;
   index_t *ptr=NULL, *idx=NULL, *degree_set=NULL, *offset_set=NULL;
   index_t *send_ptr=NULL, *recv_ptr=NULL, *recv_idx=NULL;
   index_t *temp=NULL, *where_p=NULL;
   index_t num_Pcouple_cols, num_Rcouple_cols, numNeighbors;
   index_t i, j, j_ub, k, p, iptr, iptr_ub, icb, irb;
   index_t block_size, copy_block_size, sum, offset, len, msgs;
   double  *val=NULL, *data_set=NULL, *recv_val=NULL;
   index_t *shared=NULL;
   #ifdef ESYS_MPI
     MPI_Request* mpi_requests=NULL;
     MPI_Status* mpi_stati=NULL;
   #else
     int *mpi_requests=NULL, *mpi_stati=NULL;
   #endif

   /* get main_block of R from the transpose of P->mainBlock */
   SparseMatrix_ptr main_block(P->mainBlock->getTranspose());

   /* prepare "ptr" for the col_coupleBlock of R, start with get info about
      the degree_set (associated with "ptr"), offset_set (associated with
      "idx" and data_set (associated with "val") to be sent to other ranks */
   SparseMatrix_ptr couple_block(P->col_coupleBlock);
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
   SharedComponents_ptr send(P->col_coupler->connector->send);
   SharedComponents_ptr recv(P->col_coupler->connector->recv);
   recv_ptr = new index_t[send->numSharedComponents];
   for (p=0; p<send->neighbour.size(); p++) {
     i = send->offsetInShared[p];
     j = send->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
#ifdef ESYS_MPI
        MPI_Irecv(&(recv_ptr[i]), k, MPI_INT, send->neighbour[p],
                mpi_info->counter()+send->neighbour[p],
                mpi_info->comm, &mpi_requests[msgs]);
#endif
        msgs++;
     }
   }

   for (p=0; p<recv->neighbour.size(); p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = j - i;
     if (k > 0) {
#ifdef ESYS_MPI
        MPI_Issend(&degree_set[i], k, MPI_INT, recv->neighbour[p],
                mpi_info->counter()+rank, mpi_info->comm,
                &mpi_requests[msgs]);
#endif
        msgs++;
     }
   }

#ifdef ESYS_MPI
   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   mpi_info->incCounter(size);
#endif

   delete[] degree_set;
   degree_set = new index_t[send->neighbour.size()];
   memset(degree_set, 0, sizeof(index_t)*send->neighbour.size());
   for (p=0, sum=0; p<send->neighbour.size(); p++) {
     iptr_ub = send->offsetInShared[p+1];
     for (iptr = send->offsetInShared[p]; iptr < iptr_ub; iptr++) {
        degree_set[p] += recv_ptr[iptr];
     }
     sum += degree_set[p];
   }

   /* send/receive offset_set and data_set to build the "idx" and "val"
      for R->col_coupleBlock */
   msgs = 0;
   recv_idx = new index_t[sum];
   recv_val = new double[sum * block_size];
   for (p=0, offset=0; p<send->neighbour.size(); p++) {
     if (degree_set[p]) {
#ifdef ESYS_MPI
        MPI_Irecv(&recv_idx[offset], degree_set[p], MPI_INT,
                send->neighbour[p], mpi_info->counter()+send->neighbour[p],
                mpi_info->comm, &mpi_requests[msgs]);
        msgs++;
        MPI_Irecv(&recv_val[offset*block_size], degree_set[p] * block_size,
                MPI_DOUBLE, send->neighbour[p],
                mpi_info->counter()+send->neighbour[p]+size,
                mpi_info->comm, &mpi_requests[msgs]);
        offset += degree_set[p];
#endif
        msgs++;
     }
   }

   for (p=0; p<recv->neighbour.size(); p++) {
     i = recv->offsetInShared[p];
     j = recv->offsetInShared[p+1];
     k = send_ptr[j] - send_ptr[i];
     if (k > 0) {
        #ifdef ESYS_MPI
        MPI_Issend(&offset_set[send_ptr[i]], k, MPI_INT,
                recv->neighbour[p], mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[msgs]);
        msgs++;
        MPI_Issend(&data_set[send_ptr[i]*block_size], k*block_size, MPI_DOUBLE,
                recv->neighbour[p], mpi_info->counter()+rank+size,
                mpi_info->comm, &mpi_requests[msgs]);
        #endif
        msgs++;
     }
   }

   len = send->numSharedComponents;
   temp = new index_t[len];
   memset(temp, 0, sizeof(index_t)*len);
   for (p=1; p<len; p++) {
     temp[p] = temp[p-1] + recv_ptr[p-1];
   }

#ifdef ESYS_MPI
   MPI_Waitall(msgs, mpi_requests, mpi_stati);
   mpi_info->incCounter(2*size);
#endif
   delete[] degree_set;
   delete[] offset_set;
   delete[] data_set;
   delete[] send_ptr;
   delete[] mpi_requests;
   delete[] mpi_stati;

   /* construct "ptr", "idx" and "val" for R->col_coupleBlock */
   ptr = new index_t[n_C + 1];
   idx = new index_t[sum];
   val = new double[sum*block_size];
   ptr[0] = 0;
   for (i=0; i<n_C; i++) {
     icb = 0;
     for (p=0; p<send->neighbour.size(); p++) {
        k = send->offsetInShared[p+1];
        for (j = send->offsetInShared[p]; j<k; j++) {
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
     qsort(recv_idx, (size_t)sum, sizeof(index_t), util::comparIndex);
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
                                sizeof(index_t), util::comparIndex);
        idx[i] = (index_t)(where_p - recv_idx);
     }
   }

   /* prepare the receiver for the col_connector */
   const std::vector<index_t> dist(P->pattern->output_distribution->first_component);
   std::vector<index_t> offsetInShared(size+1);
   shared = new index_t[num_Rcouple_cols];
   numNeighbors = send->neighbour.size();
   std::vector<int> neighbour = send->neighbour;
   if (num_Rcouple_cols > 0) offset = dist[neighbour[0] + 1];
   for (i=0, p=0; i<num_Rcouple_cols; i++) {
     /* cols i is received from rank neighbor[p] when it's still smaller
        than "offset", otherwise, it is received from rank neighbor[p+1] */
     while (recv_idx[i] >= offset) {
        p++;
        offsetInShared[p] = i;
        offset = dist[neighbour[p] + 1];
     }
     shared[i] = i + n;  /* n is the number of cols in R->mainBlock */
   }
   #pragma omp parallel for private(i) schedule(static)
   for (i=p; i<numNeighbors; i++) {
     offsetInShared[i+1] = num_Rcouple_cols;
   }
   recv.reset(new SharedComponents(n, neighbour, shared, offsetInShared));
   delete[] recv_idx;

   /* prepare the sender for the col_connector */
   delete[] shared;
   numNeighbors = P->col_coupler->connector->recv->neighbour.size();
   neighbour = P->col_coupler->connector->recv->neighbour;
   shared = new index_t[n * numNeighbors];
   Pattern_ptr couple_pattern(P->col_coupleBlock->pattern);
   sum=0;
   offsetInShared.assign(size+1, 0);
   for (p = 0; p < numNeighbors; p++) {
     j = P->col_coupler->connector->recv->offsetInShared[p];
     j_ub = P->col_coupler->connector->recv->offsetInShared[p+1];
     for (i = 0; i < n; i++) {
        iptr = couple_pattern->ptr[i];
        iptr_ub = couple_pattern->ptr[i+1];
        for (; iptr < iptr_ub; iptr++) {
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
   send.reset(new SharedComponents(n, neighbour, shared, offsetInShared));

   // build the col_connector based on sender and receiver
   col_connector.reset(new Connector(send, recv));
   delete[] shared;

   couple_pattern.reset(new Pattern(MATRIX_FORMAT_DEFAULT, n_C,
                        num_Rcouple_cols, ptr, idx));

   input_dist.reset(new escript::Distribution(mpi_info, dist));
   output_dist.reset(new escript::Distribution(mpi_info, P->pattern->input_distribution->first_component));

    /* now we need to create the System Matrix
       TO BE FIXED: at this stage, we only construction col_couple_pattern
       and col_connector for Restriction matrix R. To be completed,
       row_couple_pattern and row_connector need to be constructed as well */
    SystemMatrix_ptr out;
    SystemMatrixPattern_ptr pattern;
    pattern.reset(new SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
              output_dist, input_dist, main_block->pattern, couple_pattern,
              couple_pattern, col_connector, col_connector));
    out.reset(new SystemMatrix(MATRIX_FORMAT_DIAGONAL_BLOCK, pattern,
              row_block_size, col_block_size, false,
              P->getRowFunctionSpace(), P->getColumnFunctionSpace()));

    /* now fill in the matrix */
    memcpy(out->mainBlock->val, main_block->val,
                main_block->len * sizeof(double));
    memcpy(out->col_coupleBlock->val, val,
                out->col_coupleBlock->len * sizeof(double));
    delete[] val;
    return out;
}

} // namespace paso

