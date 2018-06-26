/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

/* Paso: defines AMG Interpolation                            */

/****************************************************************************/

/* Author: l.gao@uq.edu.au                                    */

/****************************************************************************/

#include "Paso.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "SparseMatrix.h"

#include <cstring> // memcpy

namespace paso {

/****************************************************************************

    Methods necessary for AMG preconditioner

    construct n_C x n_C interpolation operator A_C from matrix A
    and prolongation matrix P.

    The coarsening operator A_C is defined by RAP where R=P^T.

*****************************************************************************/

/* Extend system matrix B with extra two sparse matrices:
        B_ext_main and B_ext_couple
   The combination of this two sparse matrices represents
   the portion of B that is stored on neighbour procs and
   needed locally for triple matrix product (B^T) A B.

   FOR NOW, we assume that the system matrix B has a NULL
   row_coupleBlock and a NULL remote_coupleBlock. Based on
   this assumption, we use link row_coupleBlock for sparse
   matrix B_ext_main, and link remote_coupleBlock for sparse
   matrix B_ext_couple.

   To be specific, B_ext (on processor k) are group of rows
   in B, where the list of rows from processor q is given by
        A->col_coupler->connector->send->shared[rPtr...]
        rPtr=A->col_coupler->connector->send->OffsetInShared[k]
   on q.
*/
void Preconditioner_AMG_extendB(SystemMatrix_ptr A, SystemMatrix_ptr B)
{
    if (A->mpi_info->size == 1) return;

    if (B->remote_coupleBlock.get()) {
        throw PasoException("Preconditioner_AMG_extendB: the link to "
                            "remote_coupleBlock has already been set.");
    }
#ifdef ESYS_MPI
    B->row_coupleBlock.reset();

    Pattern_ptr pattern_main, pattern_couple;
    Coupler_ptr<real_t> coupler;
    double *ptr_val=NULL;
    index_t *global_id=NULL, *cols_array=NULL, *ptr_ptr=NULL, *ptr_idx=NULL;
    index_t *ptr_main=NULL, *ptr_couple=NULL, *idx_main=NULL, *idx_couple=NULL;
    index_t *idx_m=NULL, *idx_c=NULL;
    index_t i, j, k, m, q, j_ub, k_lb, k_ub, m_lb, m_ub, l_m, l_c, i0;
    index_t max_num_cols;
    index_t row, neighbor;
    const int rank = A->mpi_info->rank;
    const int size = A->mpi_info->size;

    // sending/receiving unknown's global ID
    dim_t num_main_cols = B->mainBlock->numCols;
    double* cols = new double[num_main_cols];
    index_t offset = B->col_distribution->first_component[rank];

#pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_main_cols; ++i)
        cols[i] = offset + i;

    if (B->global_id == NULL) {
        coupler.reset(new Coupler<real_t>(B->col_coupler->connector, 1, A->mpi_info));
        coupler->startCollect(cols);
    }

    index_t* recv_buf = new index_t[size];
    int* recv_degree = new int[size];
    int* recv_offset = new int[size+1];

#pragma omp parallel for private(i) schedule(static)
    for (i=0; i<size; i++) {
        recv_buf[i] = 0;
        recv_degree[i] = 1;
        recv_offset[i] = i;
    }

    const index_t block_size = B->block_size;
    const size_t block_size_size = block_size * sizeof(double);
    dim_t num_couple_cols = B->col_coupleBlock->numCols;
    SharedComponents_ptr send(A->col_coupler->connector->send);
    SharedComponents_ptr recv(A->col_coupler->connector->recv);
    const int num_neighbors = send->neighbour.size();
    index_t p = send->offsetInShared[num_neighbors];
    index_t len = p * B->col_distribution->first_component[size];
    double* send_buf = new double[len * block_size];
    index_t* send_idx = new index_t[len];
    int* send_offset = new int[(p+1)*2];
    int* send_degree = new int[num_neighbors];
    i = num_main_cols + num_couple_cols;
    double* send_m = new double[i * block_size];
    double* send_c = new double[i * block_size];
    idx_m = new index_t[i];
    idx_c = new index_t[i];

    /* waiting for receiving unknown's global ID */
    if (B->global_id == NULL) {
        coupler->finishCollect();
        global_id = new index_t[num_couple_cols];
#pragma omp parallel for private(i) schedule(static)
        for (i=0; i<num_couple_cols; ++i)
            global_id[i] = coupler->recv_buffer[i];
        B->global_id = global_id;
    } else
        global_id = B->global_id;

    /* distribute the number of cols in current col_coupleBlock to all ranks */
    MPI_Allgatherv(&num_couple_cols, 1, MPI_INT, recv_buf, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm);

    /* distribute global_ids of cols to be considered to all ranks*/
    len = 0;
    max_num_cols = 0;
    for (i=0; i<size; i++){
        recv_degree[i] = recv_buf[i];
        recv_offset[i] = len;
        len += recv_buf[i];
        if (max_num_cols < recv_buf[i])
            max_num_cols = recv_buf[i];
    }
    recv_offset[size] = len;
    cols_array = new index_t[len];
    MPI_Allgatherv(global_id, num_couple_cols, MPI_INT, cols_array, recv_degree, recv_offset, MPI_INT, A->mpi_info->comm);

    // first, prepare the ptr_ptr to be received
    q = recv->neighbour.size();
    len = recv->offsetInShared[q];
    ptr_ptr = new index_t[(len+1) * 2];
    for (p=0; p<q; p++) {
        row = recv->offsetInShared[p];
        m = recv->offsetInShared[p + 1];
        MPI_Irecv(&(ptr_ptr[2*row]), 2 * (m-row), MPI_INT, recv->neighbour[p],
                A->mpi_info->counter()+recv->neighbour[p],
                A->mpi_info->comm,
                &(A->col_coupler->mpi_requests[p]));
    }

    // now prepare the rows to be sent (the degree, the offset and the data)
    len = 0;
    i0 = 0;
    for (p=0; p<num_neighbors; p++) {
        i = i0;
        neighbor = send->neighbour[p];
        m_lb = B->col_distribution->first_component[neighbor];
        m_ub = B->col_distribution->first_component[neighbor + 1];
        j_ub = send->offsetInShared[p + 1];
        for (j=send->offsetInShared[p]; j<j_ub; j++) {
            row = send->shared[j];
            l_m = 0;
            l_c = 0;
            k_ub = B->col_coupleBlock->pattern->ptr[row + 1];
            k_lb = B->col_coupleBlock->pattern->ptr[row];

            /* check part of col_coupleBlock for data to be passed @row */
            for (k=k_lb; k<k_ub; k++) {
                m = global_id[B->col_coupleBlock->pattern->index[k]];
                if (m > offset) break;
                if (m>= m_lb && m < m_ub) {
                    /* data to be passed to sparse matrix B_ext_main */
                    idx_m[l_m] = m - m_lb;
                    memcpy(&(send_m[l_m*block_size]), &(B->col_coupleBlock->val[block_size*k]), block_size_size);
                    l_m++;
                } else {
                    /* data to be passed to sparse matrix B_ext_couple */
                    idx_c[l_c] = m;
                    memcpy(&(send_c[l_c*block_size]), &(B->col_coupleBlock->val[block_size*k]), block_size_size);
                    l_c++;
                }
            }
            k_lb = k;

            /* check mainBlock for data to be passed @row to sparse
            matrix B_ext_couple */
            k_ub = B->mainBlock->pattern->ptr[row + 1];
            k = B->mainBlock->pattern->ptr[row];
            memcpy(&(send_c[l_c*block_size]), &(B->mainBlock->val[block_size*k]), block_size_size * (k_ub-k));
            for (; k<k_ub; k++) {
                m = B->mainBlock->pattern->index[k] + offset;
                idx_c[l_c] = m;
                l_c++;
            }

            /* check the rest part of col_coupleBlock for data to
            be passed @row to sparse matrix B_ext_couple */
            k = k_lb;
            k_ub = B->col_coupleBlock->pattern->ptr[row + 1];
            for (k=k_lb; k<k_ub; k++) {
                m = global_id[B->col_coupleBlock->pattern->index[k]];
                if (m>= m_lb && m < m_ub) {
                    /* data to be passed to sparse matrix B_ext_main */
                    idx_m[l_m] = m - m_lb;
                    memcpy(&(send_m[l_m*block_size]), &(B->col_coupleBlock->val[block_size*k]), block_size_size);
                    l_m++;
                } else {
                    /* data to be passed to sparse matrix B_ext_couple */
                    idx_c[l_c] = m;
                    memcpy(&(send_c[l_c*block_size]), &(B->col_coupleBlock->val[block_size*k]), block_size_size);
                    l_c++;
                }
            }

            memcpy(&(send_buf[len*block_size]), send_m, block_size_size*l_m);
            memcpy(&(send_idx[len]), idx_m, l_m * sizeof(index_t));
            send_offset[2*i] = l_m;
            len += l_m;
            memcpy(&(send_buf[len*block_size]), send_c, block_size_size*l_c);
            memcpy(&(send_idx[len]), idx_c, l_c * sizeof(index_t));
            send_offset[2*i+1] = l_c;
            len += l_c;
            i++;
        }

        /* sending */
        MPI_Issend(&send_offset[2*i0], 2*(i-i0), MPI_INT, send->neighbour[p],
                A->mpi_info->counter()+rank,
                A->mpi_info->comm,
                &A->col_coupler->mpi_requests[p+recv->neighbour.size()]);
        send_degree[p] = len;
        i0 = i;
    }
    delete[] send_m;
    delete[] send_c;
    delete[] idx_m;
    delete[] idx_c;

    q = recv->neighbour.size();
    len = recv->offsetInShared[q];
    ptr_main = new index_t[(len+1)];
    ptr_couple = new index_t[(len+1)];

    MPI_Waitall(A->col_coupler->connector->send->neighbour.size() +
                    A->col_coupler->connector->recv->neighbour.size(),
                A->col_coupler->mpi_requests, A->col_coupler->mpi_stati);
    A->mpi_info->incCounter(size);

    j = 0;
    k = 0;
    ptr_main[0] = 0;
    ptr_couple[0] = 0;
    for (i=0; i<len; i++) {
        j += ptr_ptr[2*i];
        k += ptr_ptr[2*i+1];
        ptr_main[i+1] = j;
        ptr_couple[i+1] = k;
    }

    delete[] ptr_ptr;
    idx_main = new index_t[j];
    idx_couple = new index_t[k];
    ptr_idx = new index_t[j+k];
    ptr_val = new double[(j+k) * block_size];

    /* send/receive index array */
    j=0;
    k_ub = 0;
    for (p=0; p<recv->neighbour.size(); p++) {
        k = recv->offsetInShared[p];
        m = recv->offsetInShared[p+1];
        i = ptr_main[m] - ptr_main[k] + ptr_couple[m] - ptr_couple[k];
        if (i > 0) {
            k_ub ++;
            MPI_Irecv(&(ptr_idx[j]), i, MPI_INT, recv->neighbour[p],
                A->mpi_info->counter()+recv->neighbour[p],
                A->mpi_info->comm,
                &(A->col_coupler->mpi_requests[p]));
        }
        j += i;
    }

    j=0;
    k_ub = 0;
    for (p=0; p<num_neighbors; p++) {
        i = send_degree[p] - j;
        if (i > 0){
            k_ub ++;
            MPI_Issend(&(send_idx[j]), i, MPI_INT, send->neighbour[p],
                A->mpi_info->counter()+rank,
                A->mpi_info->comm,
                &(A->col_coupler->mpi_requests[p+recv->neighbour.size()]));
        }
        j = send_degree[p];
    }

    MPI_Waitall(A->col_coupler->connector->send->neighbour.size() +
                    A->col_coupler->connector->recv->neighbour.size(),
                A->col_coupler->mpi_requests,
                A->col_coupler->mpi_stati);
    A->mpi_info->incCounter(size);

#pragma omp parallel for private(i,j,k,m,p) schedule(static)
    for (i=0; i<len; i++) {
        j = ptr_main[i];
        k = ptr_main[i+1];
        m = ptr_couple[i];
        for (p=j; p<k; p++) {
            idx_main[p] = ptr_idx[m+p];
        }
        j = ptr_couple[i+1];
        for (p=m; p<j; p++) {
            idx_couple[p] = ptr_idx[k+p];
        }
    }
    delete[] ptr_idx;

    /* allocate pattern and sparsematrix for B_ext_main */
    pattern_main.reset(new Pattern(B->col_coupleBlock->pattern->type,
                len, num_main_cols, ptr_main, idx_main));
    B->row_coupleBlock.reset(new SparseMatrix(B->col_coupleBlock->type,
                pattern_main, B->row_block_size, B->col_block_size,
                false));

    /* allocate pattern and sparsematrix for B_ext_couple */
    pattern_couple.reset(new Pattern(B->col_coupleBlock->pattern->type,
                len, B->col_distribution->first_component[size],
                ptr_couple, idx_couple));
    B->remote_coupleBlock.reset(new SparseMatrix(B->col_coupleBlock->type,
                pattern_couple, B->row_block_size, B->col_block_size,
                false));

    /* send/receive value array */
    j=0;
    for (p=0; p<recv->neighbour.size(); p++) {
        k = recv->offsetInShared[p];
        m = recv->offsetInShared[p+1];
        i = ptr_main[m] - ptr_main[k] + ptr_couple[m] - ptr_couple[k];
        if (i > 0)
            MPI_Irecv(&(ptr_val[j]), i * block_size,
                MPI_DOUBLE, recv->neighbour[p],
                A->mpi_info->counter()+recv->neighbour[p],
                A->mpi_info->comm,
                &(A->col_coupler->mpi_requests[p]));
        j += (i * block_size);
    }

    j=0;
    for (p=0; p<num_neighbors; p++) {
        i = send_degree[p] - j;
        if (i > 0)
            MPI_Issend(&send_buf[j*block_size], i*block_size, MPI_DOUBLE,
                    send->neighbour[p], A->mpi_info->counter()+rank,
                    A->mpi_info->comm,
                    &A->col_coupler->mpi_requests[p+recv->neighbour.size()]);
        j = send_degree[p];
    }

    MPI_Waitall(A->col_coupler->connector->send->neighbour.size() +
                A->col_coupler->connector->recv->neighbour.size(),
                A->col_coupler->mpi_requests, A->col_coupler->mpi_stati);
    A->mpi_info->incCounter(size);

#pragma omp parallel for private(i,j,k,m,p) schedule(static)
    for (i=0; i<len; i++) {
        j = ptr_main[i];
        k = ptr_main[i+1];
        m = ptr_couple[i];
        for (p=j; p<k; p++) {
            memcpy(&(B->row_coupleBlock->val[p*block_size]), &(ptr_val[(m+p)*block_size]), block_size_size);
        }
        j = ptr_couple[i+1];
        for (p=m; p<j; p++) {
            memcpy(&(B->remote_coupleBlock->val[p*block_size]), &(ptr_val[(k+p)*block_size]), block_size_size);
        }
    }

    delete[] ptr_val;
    delete[] cols;
    delete[] cols_array;
    delete[] recv_offset;
    delete[] recv_degree;
    delete[] recv_buf;
    delete[] send_buf;
    delete[] send_offset;
    delete[] send_degree;
    delete[] send_idx;
#endif // ESYS_MPI
}

/* As defined, sparse matrix (let's called it T) defined by T(ptr, idx, val)
   has the same number of rows as P->col_coupleBlock->numCols. Now, we need
   to copy block of data in T to neighbour processors, defined by
        P->col_coupler->connector->recv->neighbour[k] where k is in
        [0, P->col_coupler->connector->recv->numNeighbours).
   Rows to be copied to neighbour processor k is in the list defined by
        P->col_coupler->connector->recv->offsetInShared[k] ...
        P->col_coupler->connector->recv->offsetInShared[k+1]  */
void Preconditioner_AMG_CopyRemoteData(SystemMatrix_ptr P,
        index_t **p_ptr, index_t **p_idx, double **p_val,
        index_t *global_id, index_t block_size)
{
    index_t i, j, p, m, n;
    index_t *ptr=*p_ptr, *idx=*p_idx;
    double  *val=*p_val;
#ifdef ESYS_MPI
    int rank = P->mpi_info->rank;
    int size = P->mpi_info->size;
#endif

    SharedComponents_ptr send(P->col_coupler->connector->recv);
    SharedComponents_ptr recv(P->col_coupler->connector->send);
    int send_neighbors = send->neighbour.size();
    int recv_neighbors = recv->neighbour.size();
    dim_t send_rows = P->col_coupleBlock->numCols;
    dim_t recv_rows = recv->offsetInShared[recv_neighbors];

    index_t* send_degree = new index_t[send_rows];
    index_t* recv_ptr = new index_t[recv_rows + 1];
#pragma omp for schedule(static) private(i)
    for (i=0; i<send_rows; i++)
        send_degree[i] = ptr[i+1] - ptr[i];

    // First, send/receive the degree
    for (p = 0; p < recv_neighbors; p++) { // Receiving
        m = recv->offsetInShared[p];
        n = recv->offsetInShared[p+1];
#ifdef ESYS_MPI
        MPI_Irecv(&recv_ptr[m], n-m, MPI_INT, recv->neighbour[p],
                  P->mpi_info->counter() + recv->neighbour[p],
                  P->mpi_info->comm, &P->col_coupler->mpi_requests[p]);
#endif
    }
    for (p = 0; p < send_neighbors; p++) { // Sending
        m = send->offsetInShared[p];
        n = send->offsetInShared[p+1];
#ifdef ESYS_MPI
        MPI_Issend(&send_degree[m], n-m, MPI_INT, send->neighbour[p],
                   P->mpi_info->counter() + rank, P->mpi_info->comm,
                   &P->col_coupler->mpi_requests[p+recv_neighbors]);
#endif
    }
#ifdef ESYS_MPI
    P->mpi_info->incCounter(size);
    MPI_Waitall(send_neighbors+recv_neighbors, P->col_coupler->mpi_requests,
                P->col_coupler->mpi_stati);
#endif

    delete[] send_degree;
    m = util::cumsum(recv_rows, recv_ptr);
    recv_ptr[recv_rows] = m;
    index_t* recv_idx = new index_t[m];
    double* recv_val = new double[m * block_size];

    // Next, send/receive the index array
    j = 0;
    for (p=0; p<recv_neighbors; p++) { // Receiving
        m = recv->offsetInShared[p];
        n = recv->offsetInShared[p+1];
        i = recv_ptr[n] - recv_ptr[m];
#ifdef ESYS_MPI
        if (i > 0) {
            MPI_Irecv(&recv_idx[j], i, MPI_INT, recv->neighbour[p],
                    P->mpi_info->counter() + recv->neighbour[p],
                    P->mpi_info->comm, &P->col_coupler->mpi_requests[p]);
        }
#endif
        j += i;
    }

    j = 0;
    for (p=0; p<send_neighbors; p++) { /* Sending */
        m = send->offsetInShared[p];
        n = send->offsetInShared[p+1];
        i = ptr[n] - ptr[m];
        if (i > 0) {
#ifdef ESYS_MPI
            MPI_Issend(&idx[j], i, MPI_INT, send->neighbour[p],
                       P->mpi_info->counter() + rank, P->mpi_info->comm,
                       &P->col_coupler->mpi_requests[p+recv_neighbors]);
#endif
            j += i;
        }
    }
#ifdef ESYS_MPI
    P->mpi_info->incCounter(size);
    MPI_Waitall(send_neighbors+recv_neighbors, P->col_coupler->mpi_requests,
                P->col_coupler->mpi_stati);
#endif

    // Last, send/receive the data array
    j = 0;
    for (p=0; p<recv_neighbors; p++) { /* Receiving */
        m = recv->offsetInShared[p];
        n = recv->offsetInShared[p+1];
        i = recv_ptr[n] - recv_ptr[m];
#ifdef ESYS_MPI
        if (i > 0)
            MPI_Irecv(&recv_val[j], i*block_size, MPI_DOUBLE, recv->neighbour[p],
                P->mpi_info->counter() + recv->neighbour[p],
                P->mpi_info->comm, &P->col_coupler->mpi_requests[p]);
#endif
        j += (i*block_size);
    }

    j = 0;
    for (p=0; p<send_neighbors; p++) { /* Sending */
        m = send->offsetInShared[p];
        n = send->offsetInShared[p+1];
        i = ptr[n] - ptr[m];
        if (i >0) {
#ifdef ESYS_MPI
            MPI_Issend(&val[j], i * block_size, MPI_DOUBLE, send->neighbour[p],
                       P->mpi_info->counter() + rank, P->mpi_info->comm,
                       &P->col_coupler->mpi_requests[p+recv_neighbors]);
#endif
            j += i * block_size;
        }
    }
#ifdef ESYS_MPI
    P->mpi_info->incCounter(size);
    MPI_Waitall(send_neighbors+recv_neighbors, P->col_coupler->mpi_requests,
                P->col_coupler->mpi_stati);
#endif

    // Clean up and return with received ptr, index and data arrays
    delete[] ptr;
    delete[] idx;
    delete[] val;
    *p_ptr = recv_ptr;
    *p_idx = recv_idx;
    *p_val = recv_val;
}

SystemMatrix_ptr Preconditioner_AMG_buildInterpolationOperator(
        SystemMatrix_ptr A, SystemMatrix_ptr P,
        SystemMatrix_ptr R)
{
   escript::JMPI& mpi_info=A->mpi_info;
   SystemMatrix_ptr out;
   SystemMatrixPattern_ptr pattern;
   escript::Distribution_ptr input_dist, output_dist;
   Connector_ptr col_connector, row_connector;
   const dim_t row_block_size=A->row_block_size;
   const dim_t col_block_size=A->col_block_size;
   const dim_t block_size = A->block_size;
   const dim_t P_block_size = P->block_size;
   const double ZERO = 0.0;
   double *RAP_main_val=NULL, *RAP_couple_val=NULL, *RAP_ext_val=NULL;
   double rtmp, *RAP_val, *RA_val, *R_val, *temp_val=NULL, *t1_val, *t2_val;
   index_t size=mpi_info->size, rank=mpi_info->rank;
   index_t *RAP_main_ptr=NULL, *RAP_couple_ptr=NULL, *RAP_ext_ptr=NULL;
   index_t *RAP_main_idx=NULL, *RAP_couple_idx=NULL, *RAP_ext_idx=NULL;
   index_t *row_couple_ptr=NULL, *row_couple_idx=NULL;
   index_t *Pcouple_to_Pext=NULL, *Pext_to_RAP=NULL, *Pcouple_to_RAP=NULL;
   index_t *temp=NULL, *global_id_P=NULL, *global_id_RAP=NULL;
   index_t *shared=NULL, *P_marker=NULL, *A_marker=NULL;
   index_t sum, i, j, k, iptr, irb, icb, ib;
   index_t num_Pmain_cols, num_Pcouple_cols, num_Pext_cols;
   index_t num_A_cols, row_marker, num_RAPext_cols, num_Acouple_cols, offset;
   index_t i_r, i_c, i1, i2, j1, j1_ub, j2, j2_ub, j3, j3_ub, num_RAPext_rows;
   index_t row_marker_ext, *where_p=NULL;
   index_t **send_ptr=NULL, **send_idx=NULL;
   dim_t p, num_neighbors;
   dim_t *recv_len=NULL, *send_len=NULL, *len=NULL;

/*   if (!(P->type & MATRIX_FORMAT_DIAGONAL_BLOCK))
     return Preconditioner_AMG_buildInterpolationOperatorBlock(A, P, R);*/

   /* two sparse matrices R_main and R_couple will be generate, as the
      transpose of P_main and P_col_couple, respectively. Note that,
      R_couple is actually the row_coupleBlock of R (R=P^T) */
   SparseMatrix_ptr R_main(P->mainBlock->getTranspose());
   SparseMatrix_ptr R_couple;
   if (size > 1)
     R_couple = P->col_coupleBlock->getTranspose();

   /* generate P_ext, i.e. portion of P that is stored on neighbour procs
      and needed locally for triple matrix product RAP
      to be specific, P_ext (on processor k) are group of rows in P, where
      the list of rows from processor q is given by
        A->col_coupler->connector->send->shared[rPtr...]
        rPtr=A->col_coupler->connector->send->OffsetInShared[k]
      on q.
      P_ext is represented by two sparse matrices P_ext_main and
      P_ext_couple */
   Preconditioner_AMG_extendB(A, P);

   /* count the number of cols in P_ext_couple, resize the pattern of
      sparse matrix P_ext_couple with new compressed order, and then
      build the col id mapping from P->col_coupleBlock to
      P_ext_couple */
   num_Pmain_cols = P->mainBlock->numCols;
   if (size > 1) {
     num_Pcouple_cols = P->col_coupleBlock->numCols;
     num_Acouple_cols = A->col_coupleBlock->numCols;
     sum = P->remote_coupleBlock->pattern->ptr[P->remote_coupleBlock->numRows];
   } else {
     num_Pcouple_cols = 0;
     num_Acouple_cols = 0;
     sum = 0;
   }
   num_A_cols = A->mainBlock->numCols + num_Acouple_cols;
   offset = P->col_distribution->first_component[rank];
   num_Pext_cols = 0;
   if (P->global_id) {
     /* first count the number of cols "num_Pext_cols" in both P_ext_couple
        and P->col_coupleBlock */
     iptr = 0;
     if (num_Pcouple_cols || sum > 0) {
        temp = new index_t[num_Pcouple_cols+sum];
        #pragma omp parallel for lastprivate(iptr) schedule(static)
        for (iptr=0; iptr<sum; iptr++){
          temp[iptr] = P->remote_coupleBlock->pattern->index[iptr];
        }
        for (j=0; j<num_Pcouple_cols; j++, iptr++){
          temp[iptr] = P->global_id[j];
        }
     }
     if (iptr) {
          qsort(temp, (size_t)iptr, sizeof(index_t), util::comparIndex);
        num_Pext_cols = 1;
        i = temp[0];
        for (j=1; j<iptr; j++) {
          if (temp[j] > i) {
            i = temp[j];
            temp[num_Pext_cols++] = i;
          }
        }
     }
     /* resize the pattern of P_ext_couple */
     if(num_Pext_cols){
        global_id_P = new index_t[num_Pext_cols];
        #pragma omp parallel for private(i) schedule(static)
        for (i=0; i<num_Pext_cols; i++)
          global_id_P[i] = temp[i];
     }
     if (num_Pcouple_cols || sum > 0)
        delete[] temp;
     #pragma omp parallel for private(i, where_p) schedule(static)
     for (i=0; i<sum; i++) {
        where_p = (index_t *)bsearch(
                        &(P->remote_coupleBlock->pattern->index[i]),
                        global_id_P, num_Pext_cols,
                        sizeof(index_t), util::comparIndex);
        P->remote_coupleBlock->pattern->index[i] =
                        (index_t)(where_p -global_id_P);
     }

     /* build the mapping */
     if (num_Pcouple_cols) {
        Pcouple_to_Pext = new index_t[num_Pcouple_cols];
        iptr = 0;
        for (i=0; i<num_Pext_cols; i++)
          if (global_id_P[i] == P->global_id[iptr]) {
            Pcouple_to_Pext[iptr++] = i;
            if (iptr == num_Pcouple_cols) break;
          }
     }
   }

   /* alloc and initialise the makers */
   sum = num_Pext_cols + num_Pmain_cols;
   P_marker = new index_t[sum];
   A_marker = new index_t[num_A_cols];
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   /* Now, count the size of RAP_ext. Start with rows in R_couple */
   sum = 0;
   for (i_r=0; i_r<num_Pcouple_cols; i_r++){
     row_marker = sum;
     /* then, loop over elements in row i_r of R_couple */
     j1_ub = R_couple->pattern->ptr[i_r+1];
     for (j1=R_couple->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_couple->pattern->index[j1];
        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];

                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
          }
        }
     }
   }

   /* Now we have the number of non-zero elements in RAP_ext, allocate
      PAP_ext_ptr, RAP_ext_idx and RAP_ext_val */
   RAP_ext_ptr = new index_t[num_Pcouple_cols+1];
   RAP_ext_idx = new index_t[sum];
   RAP_ext_val = new double[sum * block_size];
   RA_val = new double[block_size];
   RAP_val = new double[block_size];

   /* Fill in the RAP_ext_ptr, RAP_ext_idx, RAP_val */
   sum = num_Pext_cols + num_Pmain_cols;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;
   sum = 0;
   RAP_ext_ptr[0] = 0;
   for (i_r=0; i_r<num_Pcouple_cols; i_r++){
     row_marker = sum;
     /* then, loop over elements in row i_r of R_couple */
     j1_ub = R_couple->pattern->ptr[i_r+1];
     for (j1=R_couple->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_couple->pattern->index[j1];
        R_val = &(R_couple->val[j1*P_block_size]);

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];
          temp_val = &(A->col_coupleBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++)
            for (icb=0; icb<col_block_size; icb++)
                RA_val[irb+row_block_size*icb] = ZERO;
          for (irb=0; irb<P_block_size; irb++) {
            rtmp=R_val[irb];
            for (icb=0; icb<col_block_size; icb++) {
                RA_val[irb+row_block_size*icb] += rtmp * temp_val[irb+col_block_size*icb];
            }
          }

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb] = ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = i_c + offset;
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = global_id_P[i_c - num_Pmain_cols];
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          /* If entry RA[i_r, i2] is visited, no new RAP entry is created.
             Only the contributions are added. */
          } else {
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          temp_val = &(A->mainBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++)
            for (icb=0; icb<col_block_size; icb++)
                RA_val[irb+row_block_size*icb]=ZERO;
          for (irb=0; irb<P_block_size; irb++) {
            rtmp=R_val[irb];
            for (icb=0; icb<col_block_size; icb++) {
                RA_val[irb+row_block_size*icb] += rtmp * temp_val[irb+col_block_size*icb];
            }
          }

          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = i_c + offset;
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = global_id_P[i_c - num_Pmain_cols];
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
          } else {
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                    rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }
     }
     RAP_ext_ptr[i_r+1] = sum;
   }
   delete[] P_marker;
   delete[] Pcouple_to_Pext;

   /* Now we have part of RAP[r,c] where row "r" is the list of rows
      which is given by the column list of P->col_coupleBlock, and
      column "c" is the list of columns which possibly covers the
      whole column range of system matrix P. This part of data is to
      be passed to neighbouring processors, and added to corresponding
      RAP[r,c] entries in the neighbouring processors */
   Preconditioner_AMG_CopyRemoteData(P, &RAP_ext_ptr, &RAP_ext_idx,
                &RAP_ext_val, global_id_P, block_size);

   num_RAPext_rows = P->col_coupler->connector->send->numSharedComponents;
   sum = RAP_ext_ptr[num_RAPext_rows];
   num_RAPext_cols = 0;
   if (num_Pext_cols || sum > 0) {
     temp = new index_t[num_Pext_cols+sum];
     j1_ub = offset + num_Pmain_cols;
     for (i=0, iptr=0; i<sum; i++) {
        if (RAP_ext_idx[i] < offset || RAP_ext_idx[i] >= j1_ub)  /* XXX */
          temp[iptr++] = RAP_ext_idx[i];                  /* XXX */
     }
     for (j=0; j<num_Pext_cols; j++, iptr++){
        temp[iptr] = global_id_P[j];
     }

     if (iptr) {
          qsort(temp, (size_t)iptr, sizeof(index_t), util::comparIndex);
        num_RAPext_cols = 1;
        i = temp[0];
        for (j=1; j<iptr; j++) {
          if (temp[j] > i) {
            i = temp[j];
            temp[num_RAPext_cols++] = i;
          }
        }
     }
   }

   /* resize the pattern of P_ext_couple */
   if(num_RAPext_cols){
     global_id_RAP = new index_t[num_RAPext_cols];
     #pragma omp parallel for private(i) schedule(static)
     for (i=0; i<num_RAPext_cols; i++)
        global_id_RAP[i] = temp[i];
   }
   if (num_Pext_cols || sum > 0)
     delete[] temp;
   j1_ub = offset + num_Pmain_cols;
   #pragma omp parallel for private(i, where_p) schedule(static)
   for (i=0; i<sum; i++) {
     if (RAP_ext_idx[i] < offset || RAP_ext_idx[i] >= j1_ub){
        where_p = (index_t *) bsearch(&(RAP_ext_idx[i]), global_id_RAP,
/*XXX*/                 num_RAPext_cols, sizeof(index_t), util::comparIndex);
        RAP_ext_idx[i] = num_Pmain_cols + (index_t)(where_p - global_id_RAP);
     } else
        RAP_ext_idx[i] = RAP_ext_idx[i] - offset;
   }

   /* build the mapping */
   if (num_Pcouple_cols) {
     Pcouple_to_RAP = new index_t[num_Pcouple_cols];
     iptr = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (global_id_RAP[i] == P->global_id[iptr]) {
          Pcouple_to_RAP[iptr++] = i;
          if (iptr == num_Pcouple_cols) break;
        }
   }

   if (num_Pext_cols) {
     Pext_to_RAP = new index_t[num_Pext_cols];
     iptr = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (global_id_RAP[i] == global_id_P[iptr]) {
          Pext_to_RAP[iptr++] = i;
          if (iptr == num_Pext_cols) break;
        }
   }

   if (global_id_P){
     delete[] global_id_P;
     global_id_P = NULL;
   }

   /* alloc and initialise the makers */
   sum = num_RAPext_cols + num_Pmain_cols;
   P_marker = new index_t[sum];
#pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
#pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   // Now, count the size of RAP. Start with rows in R_main
   num_neighbors = P->col_coupler->connector->send->neighbour.size();
   std::vector<index_t> offsetInShared(P->col_coupler->connector->send->offsetInShared);
   shared = P->col_coupler->connector->send->shared;
   i = 0;
   j = 0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++){
     /* Mark the start of row for both main block and couple block */
     row_marker = i;
     row_marker_ext = j;

     /* Mark the diagonal element RAP[i_r, i_r], and other elements
        in RAP_ext */
     P_marker[i_r] = i;
     i++;
     for (j1 = 0; j1<num_neighbors; j1++) {
        for (j2 = offsetInShared[j1]; j2<offsetInShared[j1+1]; j2++) {
          if (shared[j2] == i_r) {
            for (k=RAP_ext_ptr[j2]; k<RAP_ext_ptr[j2+1]; k++) {
              i_c = RAP_ext_idx[k];
              if (i_c < num_Pmain_cols) {
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
              } else {
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
              }
            }
            break;
          }
        }
     }

     /* then, loop over elements in row i_r of R_main */
     j1_ub = R_main->pattern->ptr[i_r+1];
     for (j1=R_main->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_main->pattern->index[j1];

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];

                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
            }
          }
        }
     }
   }

   /* Now we have the number of non-zero elements in RAP_main and RAP_couple.
      Allocate RAP_main_ptr, RAP_main_idx and RAP_main_val for RAP_main,
      and allocate RAP_couple_ptr, RAP_couple_idx and RAP_couple_val for
      RAP_couple */
   RAP_main_ptr = new index_t[num_Pmain_cols+1];
   RAP_main_idx = new index_t[i];
   RAP_main_val = new double[i * block_size];
   RAP_couple_ptr = new index_t[num_Pmain_cols+1];
   RAP_couple_idx = new index_t[j];
   RAP_couple_val = new double[j * block_size];

   RAP_main_ptr[num_Pmain_cols] = i;
   RAP_couple_ptr[num_Pmain_cols] = j;

   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   /* Now, Fill in the data for RAP_main and RAP_couple. Start with rows
      in R_main */
   i = 0;
   j = 0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++){
     /* Mark the start of row for both main block and couple block */
     row_marker = i;
     row_marker_ext = j;
     RAP_main_ptr[i_r] = row_marker;
     RAP_couple_ptr[i_r] = row_marker_ext;

     /* Mark and setup the diagonal element RAP[i_r, i_r], and elements
        in row i_r of RAP_ext */
     P_marker[i_r] = i;
     for (ib=0; ib<block_size; ib++)
       RAP_main_val[i*block_size+ib] = ZERO;
     RAP_main_idx[i] = i_r;
     i++;

     for (j1=0; j1<num_neighbors; j1++) {
        for (j2=offsetInShared[j1]; j2<offsetInShared[j1+1]; j2++) {
          if (shared[j2] == i_r) {
            for (k=RAP_ext_ptr[j2]; k<RAP_ext_ptr[j2+1]; k++) {
              i_c = RAP_ext_idx[k];
              if (i_c < num_Pmain_cols) {
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), &(RAP_ext_val[k*block_size]), block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  t1_val = &(RAP_ext_val[k*block_size]);
                  t2_val = &(RAP_main_val[P_marker[i_c]*block_size]);
                  for (ib=0; ib<block_size; ib++)
                    t2_val[ib] += t1_val[ib];
                }
              } else {
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), &(RAP_ext_val[k*block_size]), block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  t1_val = &(RAP_ext_val[k*block_size]);
                  t2_val = &(RAP_couple_val[P_marker[i_c]*block_size]);
                  for (ib=0; ib<block_size; ib++)
                    t2_val[ib] += t1_val[ib];
                }
              }
            }
            break;
          }
        }
     }

     /* then, loop over elements in row i_r of R_main */
     j1_ub = R_main->pattern->ptr[i_r+1];
     for (j1=R_main->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_main->pattern->index[j1];
        R_val = &(R_main->val[j1*P_block_size]);

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];
          temp_val = &(A->col_coupleBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++)
            for (icb=0; icb<col_block_size; icb++)
                RA_val[irb+row_block_size*icb]=ZERO;
          for (irb=0; irb<P_block_size; irb++) {
            rtmp=R_val[irb];
            for (icb=0; icb<col_block_size; icb++) {
                RA_val[irb+col_block_size*icb] += rtmp * temp_val[irb+col_block_size*icb];
            }
          }


          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }


                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          /* If entry RA[i_r, i2] is visited, no new RAP entry is created.
             Only the contributions are added. */
          } else {
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          temp_val = &(A->mainBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++)
            for (icb=0; icb<col_block_size; icb++)
                RA_val[irb+row_block_size*icb]=ZERO;
          for (irb=0; irb<P_block_size; irb++) {
            rtmp=R_val[irb];
            for (icb=0; icb<col_block_size; icb++) {
                RA_val[irb+row_block_size*icb] += rtmp * temp_val[irb+col_block_size*icb];
            }
          }

          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          } else {
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp=temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*P_block_size]);
                for (irb=0; irb<row_block_size; irb++)
                  for (icb=0; icb<col_block_size; icb++)
                    RAP_val[irb+row_block_size*icb]=ZERO;
                for (icb=0; icb<P_block_size; icb++) {
                  rtmp = temp_val[icb];
                  for (irb=0; irb<row_block_size; irb++) {
                    RAP_val[irb+row_block_size*icb] += RA_val[irb+row_block_size*icb] * rtmp;
                  }
                }

                temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }
     }

     /* sort RAP_XXXX_idx and reorder RAP_XXXX_val accordingly */
     if (i > row_marker) {
        offset = i - row_marker;
        temp = new index_t[offset];
        #pragma omp parallel for schedule(static) private(iptr)
        for (iptr=0; iptr<offset; iptr++)
          temp[iptr] = RAP_main_idx[row_marker+iptr];
        if (offset > 0) {
            qsort(temp, (size_t)offset, sizeof(index_t), util::comparIndex);
        }
        temp_val = new double[offset * block_size];
        #pragma omp parallel for schedule(static) private(iptr,k)
        for (iptr=0; iptr<offset; iptr++){
          k = P_marker[temp[iptr]];
          memcpy(&(temp_val[iptr*block_size]), &(RAP_main_val[k*block_size]), block_size*sizeof(double));
          P_marker[temp[iptr]] = iptr + row_marker;
        }
        #pragma omp parallel for schedule(static) private(iptr)
        for (iptr=0; iptr<offset; iptr++){
          RAP_main_idx[row_marker+iptr] = temp[iptr];
          memcpy(&(RAP_main_val[(row_marker+iptr)*block_size]), &(temp_val[iptr*block_size]), block_size*sizeof(double));
        }
        delete[] temp;
        delete[] temp_val;
     }
     if (j > row_marker_ext) {
        offset = j - row_marker_ext;
        temp = new index_t[offset];
        #pragma omp parallel for schedule(static) private(iptr)
        for (iptr=0; iptr<offset; iptr++)
          temp[iptr] = RAP_couple_idx[row_marker_ext+iptr];
        if (offset > 0) {
            qsort(temp, (size_t)offset, sizeof(index_t), util::comparIndex);
        }
        temp_val = new double[offset * block_size];
        #pragma omp parallel for schedule(static) private(iptr, k)
        for (iptr=0; iptr<offset; iptr++){
          k = P_marker[temp[iptr] + num_Pmain_cols];
          memcpy(&(temp_val[iptr*block_size]), &(RAP_couple_val[k*block_size]), block_size*sizeof(double));
          P_marker[temp[iptr] + num_Pmain_cols] = iptr + row_marker_ext;
        }
        #pragma omp parallel for schedule(static) private(iptr)
        for (iptr=0; iptr<offset; iptr++){
          RAP_couple_idx[row_marker_ext+iptr] = temp[iptr];
          memcpy(&(RAP_couple_val[(row_marker_ext+iptr)*block_size]), &(temp_val[iptr*block_size]), block_size*sizeof(double));
        }
        delete[] temp;
        delete[] temp_val;
     }
   }

   delete[] RA_val;
   delete[] RAP_val;
   delete[] A_marker;
   delete[] Pext_to_RAP;
   delete[] Pcouple_to_RAP;
   delete[] RAP_ext_ptr;
   delete[] RAP_ext_idx;
   delete[] RAP_ext_val;
   R_main.reset();
   R_couple.reset();

   /* Check whether there are empty columns in RAP_couple */
   #pragma omp parallel for schedule(static) private(i)
   for (i=0; i<num_RAPext_cols; i++) P_marker[i] = 1;
   /* num of non-empty columns is stored in "k" */
   k = 0;
   j = RAP_couple_ptr[num_Pmain_cols];
   for (i=0; i<j; i++) {
     i1 = RAP_couple_idx[i];
     if (P_marker[i1]) {
        P_marker[i1] = 0;
        k++;
     }
   }

   /* empty columns is found */
   if (k < num_RAPext_cols) {
     temp = new index_t[k];
     k = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (!P_marker[i]) {
          P_marker[i] = k;
          temp[k] = global_id_RAP[i];
          k++;
        }
     #pragma omp parallel for schedule(static) private(i, i1)
     for (i=0; i<j; i++) {
        i1 = RAP_couple_idx[i];
        RAP_couple_idx[i] = P_marker[i1];
     }
     num_RAPext_cols = k;
     delete[] global_id_RAP;
     global_id_RAP = temp;
   }
   delete[] P_marker;

   /******************************************************/
   /* Start to create the coarse level System Matrix A_c */
   /******************************************************/
   /* first, prepare the sender/receiver for the col_connector */
   const std::vector<index_t> dist(P->pattern->input_distribution->first_component);
   recv_len = new dim_t[size];
   send_len = new dim_t[size];
   std::vector<int> neighbour;
   offsetInShared.clear();
   shared = new index_t[num_RAPext_cols];
   memset(recv_len, 0, sizeof(dim_t) * size);
   memset(send_len, 0, sizeof(dim_t) * size);
   offsetInShared.push_back(0);
   for (i=0, j=0, k=dist[j+1]; i<num_RAPext_cols; i++) {
     shared[i] = i + num_Pmain_cols;
     if (k <= global_id_RAP[i]) {
        if (recv_len[j] > 0) {
          neighbour.push_back(j);
          offsetInShared.push_back(i);
        }
        while (k <= global_id_RAP[i]) {
          j++;
          k = dist[j+1];
        }
     }
     recv_len[j] ++;
   }
   if (recv_len[j] > 0) {
     neighbour.push_back(j);
     offsetInShared.push_back(i);
   }

   SharedComponents_ptr recv(new SharedComponents(num_Pmain_cols, neighbour,
                             shared, offsetInShared));

#ifdef ESYS_MPI
    MPI_Alltoall(recv_len, 1, MPI_INT, send_len, 1, MPI_INT, mpi_info->comm);
    MPI_Request* mpi_requests = new MPI_Request[size*2];
    MPI_Status* mpi_stati = new MPI_Status[size*2];
#endif
   num_neighbors = 0;
   j = 0;
   neighbour.clear();
   offsetInShared.clear();
   offsetInShared.push_back(0);
   for (i=0; i<size; i++) {
     if (send_len[i] > 0) {
        neighbour.push_back(i);
        j += send_len[i];
        offsetInShared.push_back(j);
        num_neighbors++;
     }
   }
   delete[] shared;
   shared = new index_t[j];
   for (i=0, j=0; i<neighbour.size(); i++) {
     k = neighbour[i];
#ifdef ESYS_MPI
     MPI_Irecv(&shared[j], send_len[k] , MPI_INT, k, mpi_info->counter()+k,
                mpi_info->comm, &mpi_requests[i]);
#endif
     j += send_len[k];
   }
   for (i=0, j=0; i<recv->neighbour.size(); i++) {
     k = recv->neighbour[i];
#ifdef ESYS_MPI
     MPI_Issend(&(global_id_RAP[j]), recv_len[k], MPI_INT, k,
                mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[i+num_neighbors]);
#endif
     j += recv_len[k];
   }
#ifdef ESYS_MPI
   mpi_info->incCounter(size);
   MPI_Waitall(num_neighbors + recv->neighbour.size(), mpi_requests, mpi_stati);
#endif

   j = offsetInShared[num_neighbors];
   offset = dist[rank];
#pragma omp parallel for schedule(static) private(i)
   for (i=0; i<j; i++) shared[i] = shared[i] - offset;

   SharedComponents_ptr send(new SharedComponents(num_Pmain_cols, neighbour,
                                                  shared, offsetInShared));

   col_connector.reset(new Connector(send, recv));
   delete[] shared;

   /* now, create row distribution (output_distri) and col
      distribution (input_distribution) */
   input_dist.reset(new escript::Distribution(mpi_info, dist));
   output_dist.reset(new escript::Distribution(mpi_info, dist));

   /* then, prepare the sender/receiver for the row_connector, first, prepare
      the information for sender */
   sum = RAP_couple_ptr[num_Pmain_cols];
   len = new dim_t[size];
   send_ptr = new index_t*[size];
   send_idx = new index_t*[size];
   #pragma omp parallel for schedule(static) private(i)
   for (i=0; i<size; i++) {
     send_ptr[i] = new index_t[num_Pmain_cols];
     send_idx[i] = new index_t[sum];
     memset(send_ptr[i], 0, sizeof(index_t) * num_Pmain_cols);
   }
   memset(len, 0, sizeof(dim_t) * size);
   recv = col_connector->recv;
   sum=0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++) {
     i1 = RAP_couple_ptr[i_r];
     i2 = RAP_couple_ptr[i_r+1];
     if (i2 > i1) {
        /* then row i_r will be in the sender of row_connector, now check
           how many neighbours i_r needs to be send to */
        for (j=i1; j<i2; j++) {
          i_c = RAP_couple_idx[j];
          /* find out the corresponding neighbor "p" of column i_c */
          for (p=0; p<recv->neighbour.size(); p++) {
            if (i_c < recv->offsetInShared[p+1]) {
              k = recv->neighbour[p];
              if (send_ptr[k][i_r] == 0) sum++;
              send_ptr[k][i_r] ++;
              send_idx[k][len[k]] = global_id_RAP[i_c];
              len[k] ++;
              break;
            }
          }
        }
     }
   }
   if (global_id_RAP) {
     delete[] global_id_RAP;
     global_id_RAP = NULL;
   }

   /* now allocate the sender */
   shared = new index_t[sum];
   memset(send_len, 0, sizeof(dim_t) * size);
   neighbour.clear();
   offsetInShared.clear();
   offsetInShared.push_back(0);
   for (p=0, k=0; p<size; p++) {
     for (i=0; i<num_Pmain_cols; i++) {
        if (send_ptr[p][i] > 0) {
          shared[k] = i;
          k++;
          send_ptr[p][send_len[p]] = send_ptr[p][i];
          send_len[p]++;
        }
     }
     if (k > offsetInShared.back()) {
        neighbour.push_back(p);
        offsetInShared.push_back(k);
     }
   }
   send.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));

   /* send/recv number of rows will be sent from current proc
      recover info for the receiver of row_connector from the sender */
#ifdef ESYS_MPI
   MPI_Alltoall(send_len, 1, MPI_INT, recv_len, 1, MPI_INT, mpi_info->comm);
#endif
   neighbour.clear();
   offsetInShared.clear();
   offsetInShared.push_back(0);
   j = 0;
   for (i=0; i<size; i++) {
     if (i != rank && recv_len[i] > 0) {
        neighbour.push_back(i);
        j += recv_len[i];
        offsetInShared.push_back(j);
     }
   }
   num_neighbors = neighbour.size();
   delete[] shared;
   delete[] recv_len;
   shared = new index_t[j];
   k = offsetInShared.back();
#pragma omp parallel for schedule(static) private(i)
   for (i=0; i<k; i++) {
     shared[i] = i + num_Pmain_cols;
   }
   recv.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));
   row_connector.reset(new Connector(send, recv));
   delete[] shared;

   /* send/recv pattern->ptr for rowCoupleBlock */
   num_RAPext_rows = offsetInShared[num_neighbors];
   row_couple_ptr = new index_t[num_RAPext_rows+1];
   for (p = 0; p < num_neighbors; p++) {
     j = offsetInShared[p];
     i = offsetInShared[p+1];
#ifdef ESYS_MPI
     MPI_Irecv(&row_couple_ptr[j], i-j, MPI_INT, neighbour[p],
               mpi_info->counter()+neighbour[p],
               mpi_info->comm, &mpi_requests[p]);
#endif
   }
   send = row_connector->send;
   for (p=0; p<send->neighbour.size(); p++) {
#ifdef ESYS_MPI
     MPI_Issend(send_ptr[send->neighbour[p]], send_len[send->neighbour[p]],
                MPI_INT, send->neighbour[p], mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[p+num_neighbors]);
#endif
   }
#ifdef ESYS_MPI
   mpi_info->incCounter(size);
   MPI_Waitall(num_neighbors + send->neighbour.size(), mpi_requests, mpi_stati);
#endif
   delete[] send_len;

   sum = 0;
   for (i=0; i<num_RAPext_rows; i++) {
     k = row_couple_ptr[i];
     row_couple_ptr[i] = sum;
     sum += k;
   }
   row_couple_ptr[num_RAPext_rows] = sum;

   /* send/recv pattern->index for rowCoupleBlock */
   k = row_couple_ptr[num_RAPext_rows];
   row_couple_idx = new index_t[k];
   for (p=0; p<num_neighbors; p++) {
     j1 = row_couple_ptr[offsetInShared[p]];
     j2 = row_couple_ptr[offsetInShared[p+1]];
#ifdef ESYS_MPI
     MPI_Irecv(&row_couple_idx[j1], j2-j1, MPI_INT, neighbour[p],
                mpi_info->counter()+neighbour[p], mpi_info->comm,
                &mpi_requests[p]);
#endif
   }
   for (p = 0; p < send->neighbour.size(); p++) {
#ifdef ESYS_MPI
     MPI_Issend(send_idx[send->neighbour[p]], len[send->neighbour[p]],
                MPI_INT, send->neighbour[p],
                mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[p+num_neighbors]);
#endif
   }
#ifdef ESYS_MPI
    mpi_info->incCounter(size);
    MPI_Waitall(num_neighbors + send->neighbour.size(), mpi_requests, mpi_stati);
    delete[] mpi_requests;
    delete[] mpi_stati;
#endif

    offset = input_dist->first_component[rank];
    k = row_couple_ptr[num_RAPext_rows];
#pragma omp parallel for schedule(static) private(i)
    for (i=0; i<k; i++) {
        row_couple_idx[i] -= offset;
    }
#pragma omp parallel for schedule(static) private(i)
    for (i=0; i<size; i++) {
        delete[] send_ptr[i];
        delete[] send_idx[i];
    }
    delete[] send_ptr;
    delete[] send_idx;
    delete[] len;

    /* Now, we can create pattern for mainBlock and coupleBlock */
    Pattern_ptr main_pattern(new Pattern(MATRIX_FORMAT_DEFAULT,
               num_Pmain_cols, num_Pmain_cols, RAP_main_ptr, RAP_main_idx));
    Pattern_ptr col_couple_pattern(new Pattern(
               MATRIX_FORMAT_DEFAULT, num_Pmain_cols, num_RAPext_cols,
               RAP_couple_ptr, RAP_couple_idx));
    Pattern_ptr row_couple_pattern(new Pattern(
               MATRIX_FORMAT_DEFAULT, num_RAPext_rows, num_Pmain_cols,
               row_couple_ptr, row_couple_idx));

    /* next, create the system matrix */
    pattern.reset(new SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                output_dist, input_dist, main_pattern, col_couple_pattern,
                row_couple_pattern, col_connector, row_connector));
    out.reset(new SystemMatrix(A->type, pattern, row_block_size,
                               col_block_size, false, A->getRowFunctionSpace(),
                               A->getColumnFunctionSpace()));

    /* finally, fill in the data*/
    memcpy(out->mainBlock->val, RAP_main_val,
                out->mainBlock->len* sizeof(double));
    memcpy(out->col_coupleBlock->val, RAP_couple_val,
                out->col_coupleBlock->len * sizeof(double));

    delete[] RAP_main_val;
    delete[] RAP_couple_val;
    return out;
}


SystemMatrix_ptr Preconditioner_AMG_buildInterpolationOperatorBlock(
        SystemMatrix_ptr A, SystemMatrix_ptr P,
        SystemMatrix_ptr R)
{
   escript::JMPI mpi_info(A->mpi_info);
   SystemMatrix_ptr out;
   SystemMatrixPattern_ptr pattern;
   escript::Distribution_ptr input_dist, output_dist;
   SharedComponents_ptr send, recv;
   Connector_ptr col_connector, row_connector;
   const dim_t row_block_size=A->row_block_size;
   const dim_t col_block_size=A->col_block_size;
   const dim_t block_size = A->block_size;
   const double ZERO = 0.0;
   double *RAP_main_val=NULL, *RAP_couple_val=NULL, *RAP_ext_val=NULL;
   double rtmp, *RAP_val, *RA_val, *R_val, *temp_val=NULL;
   index_t size=mpi_info->size, rank=mpi_info->rank;
   index_t *RAP_main_ptr=NULL, *RAP_couple_ptr=NULL, *RAP_ext_ptr=NULL;
   index_t *RAP_main_idx=NULL, *RAP_couple_idx=NULL, *RAP_ext_idx=NULL;
   index_t *row_couple_ptr=NULL, *row_couple_idx=NULL;
   index_t *Pcouple_to_Pext=NULL, *Pext_to_RAP=NULL, *Pcouple_to_RAP=NULL;
   index_t *temp=NULL, *global_id_P=NULL, *global_id_RAP=NULL;
   index_t *shared=NULL, *P_marker=NULL, *A_marker=NULL;
   index_t sum, i, j, k, iptr, irb, icb, ib;
   index_t num_Pmain_cols, num_Pcouple_cols, num_Pext_cols;
   index_t num_A_cols, row_marker, num_RAPext_cols, num_Acouple_cols, offset;
   index_t i_r, i_c, i1, i2, j1, j1_ub, j2, j2_ub, j3, j3_ub, num_RAPext_rows;
   index_t row_marker_ext, *where_p=NULL;
   index_t **send_ptr=NULL, **send_idx=NULL;
   dim_t p, num_neighbors;
   dim_t *recv_len=NULL, *send_len=NULL, *len=NULL;
#ifdef ESYS_MPI
    MPI_Request* mpi_requests=NULL;
    MPI_Status* mpi_stati=NULL;
#else
    int *mpi_requests=NULL, *mpi_stati=NULL;
#endif

   /* two sparse matrices R_main and R_couple will be generated, as the
      transpose of P_main and P_col_couple, respectively. Note that,
      R_couple is actually the row_coupleBlock of R (R=P^T) */
   SparseMatrix_ptr R_main(P->mainBlock->getTranspose());
   SparseMatrix_ptr R_couple(P->col_coupleBlock->getTranspose());

   /* generate P_ext, i.e. portion of P that is stored on neighbor procs
      and needed locally for triple matrix product RAP
      to be specific, P_ext (on processor k) are group of rows in P, where
      the list of rows from processor q is given by
        A->col_coupler->connector->send->shared[rPtr...]
        rPtr=A->col_coupler->connector->send->offsetInShared[k]
      on q.
      P_ext is represented by two sparse matrices P_ext_main and
      P_ext_couple */
   Preconditioner_AMG_extendB(A, P);

   /* count the number of cols in P_ext_couple, resize the pattern of
      sparse matrix P_ext_couple with new compressed order, and then
      build the col id mapping from P->col_coupleBlock to
      P_ext_couple */
   num_Pmain_cols = P->mainBlock->numCols;
   num_Pcouple_cols = P->col_coupleBlock->numCols;
   num_Acouple_cols = A->col_coupleBlock->numCols;
   num_A_cols = A->mainBlock->numCols + num_Acouple_cols;
   sum = P->remote_coupleBlock->pattern->ptr[P->remote_coupleBlock->numRows];
   offset = P->col_distribution->first_component[rank];
   num_Pext_cols = 0;
   if (P->global_id) {
     /* first count the number of cols "num_Pext_cols" in both P_ext_couple
        and P->col_coupleBlock */
     iptr = 0;
     if (num_Pcouple_cols || sum > 0) {
        temp = new index_t[num_Pcouple_cols+sum];
        for (; iptr<sum; iptr++){
          temp[iptr] = P->remote_coupleBlock->pattern->index[iptr];
        }
        for (j=0; j<num_Pcouple_cols; j++, iptr++){
          temp[iptr] = P->global_id[j];
        }
     }
     if (iptr) {
          qsort(temp, (size_t)iptr, sizeof(index_t), util::comparIndex);
        num_Pext_cols = 1;
        i = temp[0];
        for (j=1; j<iptr; j++) {
          if (temp[j] > i) {
            i = temp[j];
            temp[num_Pext_cols++] = i;
          }
        }
     }

     /* resize the pattern of P_ext_couple */
     if(num_Pext_cols){
        global_id_P = new index_t[num_Pext_cols];
        for (i=0; i<num_Pext_cols; i++)
          global_id_P[i] = temp[i];
     }
     if (num_Pcouple_cols || sum > 0)
        delete[] temp;
     for (i=0; i<sum; i++) {
        where_p = (index_t *)bsearch(
                        &(P->remote_coupleBlock->pattern->index[i]),
                        global_id_P, num_Pext_cols,
                        sizeof(index_t), util::comparIndex);
        P->remote_coupleBlock->pattern->index[i] =
                        (index_t)(where_p -global_id_P);
     }

     /* build the mapping */
     if (num_Pcouple_cols) {
        Pcouple_to_Pext = new index_t[num_Pcouple_cols];
        iptr = 0;
        for (i=0; i<num_Pext_cols; i++)
          if (global_id_P[i] == P->global_id[iptr]) {
            Pcouple_to_Pext[iptr++] = i;
            if (iptr == num_Pcouple_cols) break;
          }
     }
   }

   /* alloc and initialise the makers */
   sum = num_Pext_cols + num_Pmain_cols;
   P_marker = new index_t[sum];
   A_marker = new index_t[num_A_cols];
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   /* Now, count the size of RAP_ext. Start with rows in R_couple */
   sum = 0;
   for (i_r=0; i_r<num_Pcouple_cols; i_r++){
     row_marker = sum;
     /* then, loop over elements in row i_r of R_couple */
     j1_ub = R_couple->pattern->ptr[i_r+1];
     for (j1=R_couple->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_couple->pattern->index[j1];
        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];

                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  sum++;
                }
            }
          }
        }
     }
   }

   /* Now we have the number of non-zero elements in RAP_ext, allocate
      PAP_ext_ptr, RAP_ext_idx and RAP_ext_val */
   RAP_ext_ptr = new index_t[num_Pcouple_cols+1];
   RAP_ext_idx = new index_t[sum];
   RAP_ext_val = new double[sum * block_size];
   RA_val = new double[block_size];
   RAP_val = new double[block_size];

   /* Fill in the RAP_ext_ptr, RAP_ext_idx, RAP_val */
   sum = num_Pext_cols + num_Pmain_cols;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;
   sum = 0;
   RAP_ext_ptr[0] = 0;
   for (i_r=0; i_r<num_Pcouple_cols; i_r++){
     row_marker = sum;
     /* then, loop over elements in row i_r of R_couple */
     j1_ub = R_couple->pattern->ptr[i_r+1];
     for (j1=R_couple->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_couple->pattern->index[j1];
        R_val = &(R_couple->val[j1*block_size]);

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];
          temp_val = &(A->col_coupleBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++) {
            for (icb=0; icb<col_block_size; icb++) {
                rtmp = ZERO;
                for (ib=0; ib<col_block_size; ib++) {
                  rtmp+= R_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                }
                RA_val[irb+row_block_size*icb]=rtmp;
            }
          }

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }


                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = i_c + offset;
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = global_id_P[i_c - num_Pmain_cols];
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          /* If entry RA[i_r, i2] is visited, no new RAP entry is created.
             Only the contributions are added. */
          } else {
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->remote_coupleBlock->pattern->index[j3] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          temp_val = &(A->mainBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++) {
            for (icb=0; icb<col_block_size; icb++) {
                rtmp = ZERO;
                for (ib=0; ib<col_block_size; ib++) {
                  rtmp+= R_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                }
                RA_val[irb+row_block_size*icb]=rtmp;
            }
          }

          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = i_c + offset;
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = sum;
                  memcpy(&(RAP_ext_val[sum*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_ext_idx[sum] = global_id_P[i_c - num_Pmain_cols];
                  sum++;
                } else {
                  temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
          } else {
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_Pext[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_ext_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }
     }
     RAP_ext_ptr[i_r+1] = sum;
   }
   delete[] P_marker;
   delete[] Pcouple_to_Pext;

   /* Now we have part of RAP[r,c] where row "r" is the list of rows
      which is given by the column list of P->col_coupleBlock, and
      column "c" is the list of columns which possibly covers the
      whole column range of system matrix P. This part of data is to
      be passed to neighbouring processors, and added to corresponding
      RAP[r,c] entries in the neighbouring processors */
   Preconditioner_AMG_CopyRemoteData(P, &RAP_ext_ptr, &RAP_ext_idx,
                &RAP_ext_val, global_id_P, block_size);

   num_RAPext_rows = P->col_coupler->connector->send->numSharedComponents;
   sum = RAP_ext_ptr[num_RAPext_rows];
   num_RAPext_cols = 0;
   if (num_Pext_cols || sum > 0) {
     temp = new index_t[num_Pext_cols+sum];
     j1_ub = offset + num_Pmain_cols;
     for (i=0, iptr=0; i<sum; i++) {
        if (RAP_ext_idx[i] < offset || RAP_ext_idx[i] >= j1_ub)  /* XXX */
          temp[iptr++] = RAP_ext_idx[i];                  /* XXX */
     }
     for (j=0; j<num_Pext_cols; j++, iptr++){
        temp[iptr] = global_id_P[j];
     }

     if (iptr) {
          qsort(temp, (size_t)iptr, sizeof(index_t), util::comparIndex);
        num_RAPext_cols = 1;
        i = temp[0];
        for (j=1; j<iptr; j++) {
          if (temp[j] > i) {
            i = temp[j];
            temp[num_RAPext_cols++] = i;
          }
        }
     }
   }

   /* resize the pattern of P_ext_couple */
   if(num_RAPext_cols){
     global_id_RAP = new index_t[num_RAPext_cols];
     for (i=0; i<num_RAPext_cols; i++)
        global_id_RAP[i] = temp[i];
   }
   if (num_Pext_cols || sum > 0)
     delete[] temp;
   j1_ub = offset + num_Pmain_cols;
   for (i=0; i<sum; i++) {
     if (RAP_ext_idx[i] < offset || RAP_ext_idx[i] >= j1_ub){
        where_p = (index_t *) bsearch(&(RAP_ext_idx[i]), global_id_RAP,
/*XXX*/                 num_RAPext_cols, sizeof(index_t), util::comparIndex);
        RAP_ext_idx[i] = num_Pmain_cols + (index_t)(where_p - global_id_RAP);
     } else
        RAP_ext_idx[i] = RAP_ext_idx[i] - offset;
   }

   /* build the mapping */
   if (num_Pcouple_cols) {
     Pcouple_to_RAP = new index_t[num_Pcouple_cols];
     iptr = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (global_id_RAP[i] == P->global_id[iptr]) {
          Pcouple_to_RAP[iptr++] = i;
          if (iptr == num_Pcouple_cols) break;
        }
   }

   if (num_Pext_cols) {
     Pext_to_RAP = new index_t[num_Pext_cols];
     iptr = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (global_id_RAP[i] == global_id_P[iptr]) {
          Pext_to_RAP[iptr++] = i;
          if (iptr == num_Pext_cols) break;
        }
   }

   if (global_id_P){
     delete[] global_id_P;
     global_id_P = NULL;
   }

   /* alloc and initialise the makers */
   sum = num_RAPext_cols + num_Pmain_cols;
   P_marker = new index_t[sum];
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   /* Now, count the size of RAP. Start with rows in R_main */
   num_neighbors = P->col_coupler->connector->send->neighbour.size();
   std::vector<index_t> offsetInShared = P->col_coupler->connector->send->offsetInShared;
   shared = P->col_coupler->connector->send->shared;
   i = 0;
   j = 0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++){
     /* Mark the start of row for both main block and couple block */
     row_marker = i;
     row_marker_ext = j;

     /* Mark the diagonal element RAP[i_r, i_r], and other elements
        in RAP_ext */
     P_marker[i_r] = i;
     i++;
     for (j1=0; j1<num_neighbors; j1++) {
        for (j2=offsetInShared[j1]; j2<offsetInShared[j1+1]; j2++) {
          if (shared[j2] == i_r) {
            for (k=RAP_ext_ptr[j2]; k<RAP_ext_ptr[j2+1]; k++) {
              i_c = RAP_ext_idx[k];
              if (i_c < num_Pmain_cols) {
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
              } else {
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
              }
            }
            break;
          }
        }
     }

     /* then, loop over elements in row i_r of R_main */
     j1_ub = R_main->pattern->ptr[i_r+1];
     for (j1=R_main->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_main->pattern->index[j1];

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];

          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];

                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  i++;
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  j++;
                }
            }
          }
        }
     }
   }

   /* Now we have the number of non-zero elements in RAP_main and RAP_couple.
      Allocate RAP_main_ptr, RAP_main_idx and RAP_main_val for RAP_main,
      and allocate RAP_couple_ptr, RAP_couple_idx and RAP_couple_val for
      RAP_couple */
   RAP_main_ptr = new index_t[num_Pmain_cols+1];
   RAP_main_idx = new index_t[i];
   RAP_main_val = new double[i * block_size];
   RAP_couple_ptr = new index_t[num_Pmain_cols+1];
   RAP_couple_idx = new index_t[j];
   RAP_couple_val = new double[j * block_size];

   RAP_main_ptr[num_Pmain_cols] = i;
   RAP_couple_ptr[num_Pmain_cols] = j;

   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<sum; i++) P_marker[i] = -1;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<num_A_cols; i++) A_marker[i] = -1;

   /* Now, Fill in the data for RAP_main and RAP_couple. Start with rows
      in R_main */
   i = 0;
   j = 0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++){
     /* Mark the start of row for both main block and couple block */
     row_marker = i;
     row_marker_ext = j;
     RAP_main_ptr[i_r] = row_marker;
     RAP_couple_ptr[i_r] = row_marker_ext;

     /* Mark and setup the diagonal element RAP[i_r, i_r], and elements
        in row i_r of RAP_ext */
     P_marker[i_r] = i;
     RAP_main_val[i] = ZERO;
     RAP_main_idx[i] = i_r;
     i++;

     for (j1=0; j1<num_neighbors; j1++) {
        for (j2=offsetInShared[j1]; j2<offsetInShared[j1+1]; j2++) {
          if (shared[j2] == i_r) {
            for (k=RAP_ext_ptr[j2]; k<RAP_ext_ptr[j2+1]; k++) {
              i_c = RAP_ext_idx[k];
              if (i_c < num_Pmain_cols) {
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), &(RAP_ext_val[k*block_size]), block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  temp_val = &(RAP_ext_val[k*block_size]);
                  RAP_val = &(RAP_main_val[P_marker[i_c]*block_size]);
                  for (ib=0; ib<block_size; ib++)
                    RAP_val[ib] += temp_val[ib];
                }
              } else {
                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), &(RAP_ext_val[k*block_size]), block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  temp_val = &(RAP_ext_val[k*block_size]);
                  RAP_val = &(RAP_couple_val[P_marker[i_c]*block_size]);
                  for (ib=0; ib<block_size; ib++)
                    RAP_val[ib] += temp_val[ib];
                }
              }
            }
            break;
          }
        }
     }

     /* then, loop over elements in row i_r of R_main */
     j1_ub = R_main->pattern->ptr[i_r+1];
     for (j1=R_main->pattern->ptr[i_r]; j1<j1_ub; j1++){
        i1 = R_main->pattern->index[j1];
        R_val = &(R_main->val[j1*block_size]);

        /* then, loop over elements in row i1 of A->col_coupleBlock */
        j2_ub = A->col_coupleBlock->pattern->ptr[i1+1];
        for (j2=A->col_coupleBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->col_coupleBlock->pattern->index[j2];
          temp_val = &(A->col_coupleBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++) {
            for (icb=0; icb<col_block_size; icb++) {
                rtmp = ZERO;
                for (ib=0; ib<col_block_size; ib++) {
                  rtmp+= R_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                }
                RA_val[irb+row_block_size*icb]=rtmp;
            }
          }


          /* check whether entry RA[i_r, i2] has been previously visited.
             RAP new entry is possible only if entry RA[i_r, i2] has not
             been visited yet */
          if (A_marker[i2] != i_r) {
            /* first, mark entry RA[i_r,i2] as visited */
            A_marker[i2] = i_r;

            /* then loop over elements in row i2 of P_ext_main */
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }


                /* check whether entry RAP[i_r,i_c] has been created.
                   If not yet created, create the entry and increase
                   the total number of elements in RAP_ext */
                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

            /* loop over elements in row i2 of P_ext_couple, do the same */
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          /* If entry RA[i_r, i2] is visited, no new RAP entry is created.
             Only the contributions are added. */
          } else {
            j3_ub = P->row_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->row_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->row_coupleBlock->pattern->index[j3];
                temp_val = &(P->row_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->remote_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->remote_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pext_to_RAP[P->remote_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->remote_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }

        /* now loop over elements in row i1 of A->mainBlock, repeat
           the process we have done to block A->col_coupleBlock */
        j2_ub = A->mainBlock->pattern->ptr[i1+1];
        for (j2=A->mainBlock->pattern->ptr[i1]; j2<j2_ub; j2++) {
          i2 = A->mainBlock->pattern->index[j2];
          temp_val = &(A->mainBlock->val[j2*block_size]);
          for (irb=0; irb<row_block_size; irb++) {
            for (icb=0; icb<col_block_size; icb++) {
                rtmp = ZERO;
                for (ib=0; ib<col_block_size; ib++) {
                  rtmp+= R_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                }
                RA_val[irb+row_block_size*icb]=rtmp;
            }
          }

          if (A_marker[i2 + num_Acouple_cols] != i_r) {
            A_marker[i2 + num_Acouple_cols] = i_r;
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker) {
                  P_marker[i_c] = i;
                  memcpy(&(RAP_main_val[i*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_main_idx[i] = i_c;
                  i++;
                } else {
                  temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                /* note that we need to map the column index in
                   P->col_coupleBlock back into the column index in
                   P_ext_couple */
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                if (P_marker[i_c] < row_marker_ext) {
                  P_marker[i_c] = j;
                  memcpy(&(RAP_couple_val[j*block_size]), RAP_val, block_size*sizeof(double));
                  RAP_couple_idx[j] = i_c - num_Pmain_cols;
                  j++;
                } else {
                  temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                  for (ib=0; ib<block_size; ib++) {
                    temp_val[ib] += RAP_val[ib];
                  }
                }
            }

          } else {
            j3_ub = P->mainBlock->pattern->ptr[i2+1];
            for (j3=P->mainBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = P->mainBlock->pattern->index[j3];
                temp_val = &(P->mainBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_main_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
            j3_ub = P->col_coupleBlock->pattern->ptr[i2+1];
            for (j3=P->col_coupleBlock->pattern->ptr[i2]; j3<j3_ub; j3++) {
                i_c = Pcouple_to_RAP[P->col_coupleBlock->pattern->index[j3]] + num_Pmain_cols;
                temp_val = &(P->col_coupleBlock->val[j3*block_size]);
                for (irb=0; irb<row_block_size; irb++) {
                  for (icb=0; icb<col_block_size; icb++) {
                    rtmp = ZERO;
                    for (ib=0; ib<col_block_size; ib++) {
                        rtmp+= RA_val[irb+row_block_size*ib]*temp_val[ib+col_block_size*icb];
                    }
                    RAP_val[irb+row_block_size*icb]=rtmp;
                  }
                }

                temp_val = &(RAP_couple_val[P_marker[i_c] * block_size]);
                for (ib=0; ib<block_size; ib++) {
                  temp_val[ib] += RAP_val[ib];
                }
            }
          }
        }
     }

     /* sort RAP_XXXX_idx and reorder RAP_XXXX_val accordingly */
     if (i > row_marker) {
        offset = i - row_marker;
        temp = new index_t[offset];
        for (iptr=0; iptr<offset; iptr++)
          temp[iptr] = RAP_main_idx[row_marker+iptr];
        if (offset > 0) {
            qsort(temp, (size_t)offset, sizeof(index_t), util::comparIndex);
        }
        temp_val = new double[offset * block_size];
        for (iptr=0; iptr<offset; iptr++){
          k = P_marker[temp[iptr]];
          memcpy(&(temp_val[iptr*block_size]), &(RAP_main_val[k*block_size]), block_size*sizeof(double));
          P_marker[temp[iptr]] = iptr + row_marker;
        }
        for (iptr=0; iptr<offset; iptr++){
          RAP_main_idx[row_marker+iptr] = temp[iptr];
          memcpy(&(RAP_main_val[(row_marker+iptr)*block_size]), &(temp_val[iptr*block_size]), block_size*sizeof(double));
        }
        delete[] temp;
        delete[] temp_val;
     }
     if (j > row_marker_ext) {
        offset = j - row_marker_ext;
        temp = new index_t[offset];
        for (iptr=0; iptr<offset; iptr++)
          temp[iptr] = RAP_couple_idx[row_marker_ext+iptr];
        if (offset > 0) {
            qsort(temp, (size_t)offset, sizeof(index_t), util::comparIndex);
        }
        temp_val = new double[offset * block_size];
        for (iptr=0; iptr<offset; iptr++){
          k = P_marker[temp[iptr] + num_Pmain_cols];
          memcpy(&(temp_val[iptr*block_size]), &(RAP_couple_val[k*block_size]), block_size*sizeof(double));
          P_marker[temp[iptr] + num_Pmain_cols] = iptr + row_marker_ext;
        }
        for (iptr=0; iptr<offset; iptr++){
          RAP_couple_idx[row_marker_ext+iptr] = temp[iptr];
          memcpy(&(RAP_couple_val[(row_marker_ext+iptr)*block_size]), &(temp_val[iptr*block_size]), block_size*sizeof(double));
        }
        delete[] temp;
        delete[] temp_val;
     }
   }

   delete[] RA_val;
   delete[] RAP_val;
   delete[] A_marker;
   delete[] Pext_to_RAP;
   delete[] Pcouple_to_RAP;
   delete[] RAP_ext_ptr;
   delete[] RAP_ext_idx;
   delete[] RAP_ext_val;
   R_main.reset();
   R_couple.reset();

   /* Check whether there are empty columns in RAP_couple */
   #pragma omp parallel for schedule(static) private(i)
   for (i=0; i<num_RAPext_cols; i++) P_marker[i] = 1;
   /* num of non-empty columns is stored in "k" */
   k = 0;
   j = RAP_couple_ptr[num_Pmain_cols];
   for (i=0; i<j; i++) {
     i1 = RAP_couple_idx[i];
     if (P_marker[i1]) {
        P_marker[i1] = 0;
        k++;
     }
   }

   /* empty columns is found */
   if (k < num_RAPext_cols) {
     temp = new index_t[k];
     k = 0;
     for (i=0; i<num_RAPext_cols; i++)
        if (!P_marker[i]) {
          P_marker[i] = k;
          temp[k] = global_id_RAP[i];
          k++;
        }
     for (i=0; i<j; i++) {
        i1 = RAP_couple_idx[i];
        RAP_couple_idx[i] = P_marker[i1];
     }
     num_RAPext_cols = k;
     delete[] global_id_RAP;
     global_id_RAP = temp;
   }
   delete[] P_marker;

   /******************************************************/
   /* Start to create the coarse level System Matrix A_c */
   /******************************************************/
   /* first, prepare the sender/receiver for the col_connector */
   const std::vector<index_t> dist(P->pattern->input_distribution->first_component);
   recv_len = new dim_t[size];
   send_len = new dim_t[size];
   std::vector<int> neighbour;
   offsetInShared.clear();
   shared = new index_t[num_RAPext_cols];
   memset(recv_len, 0, sizeof(dim_t) * size);
   memset(send_len, 0, sizeof(dim_t) * size);
   offsetInShared.push_back(0);
   for (i = 0, j = 0, k = dist[j+1]; i<num_RAPext_cols; i++) {
     shared[i] = i + num_Pmain_cols;
     if (k <= global_id_RAP[i]) {
        if (recv_len[j] > 0) {
          neighbour.push_back(j);
          offsetInShared.push_back(i);
          num_neighbors ++;
        }
        while (k <= global_id_RAP[i]) {
          j++;
          k = dist[j+1];
        }
     }
     recv_len[j] ++;
   }
   if (recv_len[j] > 0) {
     neighbour.push_back(j);
     offsetInShared.push_back(i);
   }
   recv.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));

#ifdef ESYS_MPI
    MPI_Alltoall(recv_len, 1, MPI_INT, send_len, 1, MPI_INT, mpi_info->comm);

    mpi_requests = new MPI_Request[size*2];
    mpi_stati = new MPI_Status[size*2];
#else
     mpi_requests=new int[size*2];
     mpi_stati=new int[size*2];
#endif
   num_neighbors = 0;
   j = 0;
   neighbour.clear();
   offsetInShared.clear();
   offsetInShared.push_back(0);
   for (i=0; i<size; i++) {
     if (send_len[i] > 0) {
        neighbour.push_back(i);
        j += send_len[i];
        offsetInShared.push_back(j);
        num_neighbors++;
     }
   }
   delete[] shared;
   shared = new index_t[j];
   for (i=0, j=0; i<num_neighbors; i++) {
     k = neighbour[i];
#ifdef ESYS_MPI
     MPI_Irecv(&shared[j], send_len[k] , MPI_INT, k, mpi_info->counter()+k,
               mpi_info->comm, &mpi_requests[i]);
#endif
     j += send_len[k];
   }
   for (i=0, j=0; i<recv->neighbour.size(); i++) {
     k = recv->neighbour[i];
#ifdef ESYS_MPI
     MPI_Issend(&global_id_RAP[j], recv_len[k], MPI_INT, k,
                mpi_info->counter()+rank, mpi_info->comm,
                &mpi_requests[i+num_neighbors]);
#endif
     j += recv_len[k];
   }
#ifdef ESYS_MPI
   mpi_info->incCounter(size);
   MPI_Waitall(num_neighbors + recv->neighbour.size(), mpi_requests, mpi_stati);
#endif

   j = offsetInShared[num_neighbors];
   offset = dist[rank];
   for (i=0; i<j; i++) shared[i] = shared[i] - offset;
   send.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));

   col_connector.reset(new Connector(send, recv));
   delete[] shared;

   /* now, create row distribution (output_distri) and col
      distribution (input_distribution) */
   input_dist.reset(new escript::Distribution(mpi_info, dist));
   output_dist.reset(new escript::Distribution(mpi_info, dist));

   /* then, prepare the sender/receiver for the row_connector, first, prepare
      the information for sender */
   sum = RAP_couple_ptr[num_Pmain_cols];
   len = new dim_t[size];
   send_ptr = new index_t*[size];
   send_idx = new index_t*[size];
   for (i=0; i<size; i++) {
     send_ptr[i] = new index_t[num_Pmain_cols];
     send_idx[i] = new index_t[sum];
     memset(send_ptr[i], 0, sizeof(index_t) * num_Pmain_cols);
   }
   memset(len, 0, sizeof(dim_t) * size);
   recv = col_connector->recv;
   sum=0;
   for (i_r=0; i_r<num_Pmain_cols; i_r++) {
     i1 = RAP_couple_ptr[i_r];
     i2 = RAP_couple_ptr[i_r+1];
     if (i2 > i1) {
        /* then row i_r will be in the sender of row_connector, now check
           how many neighbours i_r needs to be send to */
        for (j=i1; j<i2; j++) {
          i_c = RAP_couple_idx[j];
          /* find out the corresponding neighbour "p" of column i_c */
          for (p=0; p<recv->neighbour.size(); p++) {
            if (i_c < recv->offsetInShared[p+1]) {
              k = recv->neighbour[p];
              if (send_ptr[k][i_r] == 0) sum++;
              send_ptr[k][i_r] ++;
              send_idx[k][len[k]] = global_id_RAP[i_c];
              len[k] ++;
              break;
            }
          }
        }
     }
   }
   if (global_id_RAP) {
     delete[] global_id_RAP;
     global_id_RAP = NULL;
   }

   /* now allocate the sender */
   shared = new index_t[sum];
   memset(send_len, 0, sizeof(dim_t) * size);
   neighbour.clear();
   offsetInShared.clear();
   offsetInShared.push_back(0);
   for (p = 0, k = 0; p < size; p++) {
     for (i = 0; i < num_Pmain_cols; i++) {
        if (send_ptr[p][i] > 0) {
          shared[k] = i;
          k++;
          send_ptr[p][send_len[p]] = send_ptr[p][i];
          send_len[p]++;
        }
     }
     if (k > offsetInShared.back()) {
        neighbour.push_back(p);
        offsetInShared.push_back(k);
     }
   }
   send.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));

   /* send/recv number of rows will be sent from current proc
      recover info for the receiver of row_connector from the sender */
#ifdef ESYS_MPI
   MPI_Alltoall(send_len, 1, MPI_INT, recv_len, 1, MPI_INT, mpi_info->comm);
#endif
   neighbour.clear();
   offsetInShared.clear();
   num_neighbors = 0;
   offsetInShared.push_back(0);
   j = 0;
   for (i=0; i<size; i++) {
     if (i != rank && recv_len[i] > 0) {
        neighbour.push_back(i);
        j += recv_len[i];
        offsetInShared.push_back(j);
        num_neighbors ++;
     }
   }
   delete[] shared;
   delete[] recv_len;
   shared = new index_t[j];
   k = offsetInShared.back();
   for (i=0; i<k; i++) {
     shared[i] = i + num_Pmain_cols;
   }
   recv.reset(new SharedComponents(num_Pmain_cols, neighbour, shared,
                                   offsetInShared));
   row_connector.reset(new Connector(send, recv));
   delete[] shared;

   /* send/recv pattern->ptr for rowCoupleBlock */
   num_RAPext_rows = offsetInShared.back();
   row_couple_ptr = new index_t[num_RAPext_rows+1];
   for (p=0; p<num_neighbors; p++) {
     j = offsetInShared[p];
     i = offsetInShared[p+1];
#ifdef ESYS_MPI
     MPI_Irecv(&row_couple_ptr[j], i-j, MPI_INT, neighbour[p],
                mpi_info->counter()+neighbour[p],
                mpi_info->comm, &mpi_requests[p]);
#endif
   }
   send = row_connector->send;
   for (p=0; p<send->neighbour.size(); p++) {
#ifdef ESYS_MPI
     MPI_Issend(send_ptr[send->neighbour[p]], send_len[send->neighbour[p]],
                MPI_INT, send->neighbour[p],
                mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[p+num_neighbors]);
#endif
   }
#ifdef ESYS_MPI
   mpi_info->incCounter(size);
   MPI_Waitall(num_neighbors + send->neighbour.size(), mpi_requests, mpi_stati);
#endif
   delete[] send_len;

   sum = 0;
   for (i=0; i<num_RAPext_rows; i++) {
     k = row_couple_ptr[i];
     row_couple_ptr[i] = sum;
     sum += k;
   }
   row_couple_ptr[num_RAPext_rows] = sum;

   /* send/recv pattern->index for rowCoupleBlock */
   k = row_couple_ptr[num_RAPext_rows];
   row_couple_idx = new index_t[k];
   for (p=0; p<num_neighbors; p++) {
     j1 = row_couple_ptr[offsetInShared[p]];
     j2 = row_couple_ptr[offsetInShared[p+1]];
#ifdef ESYS_MPI
     MPI_Irecv(&row_couple_idx[j1], j2-j1, MPI_INT, neighbour[p],
                mpi_info->counter()+neighbour[p],
                mpi_info->comm, &mpi_requests[p]);
#endif
   }
   for (p=0; p<send->neighbour.size(); p++) {
#ifdef ESYS_MPI
     MPI_Issend(send_idx[send->neighbour[p]], len[send->neighbour[p]],
                MPI_INT, send->neighbour[p],
                mpi_info->counter()+rank,
                mpi_info->comm, &mpi_requests[p+num_neighbors]);
#endif
   }
#ifdef ESYS_MPI
   mpi_info->incCounter(size);
   MPI_Waitall(num_neighbors + send->neighbour.size(), mpi_requests, mpi_stati);
#endif

   offset = input_dist->first_component[rank];
   k = row_couple_ptr[num_RAPext_rows];
   for (i=0; i<k; i++) {
     row_couple_idx[i] -= offset;
   }

   for (i=0; i<size; i++) {
     delete[] send_ptr[i];
     delete[] send_idx[i];
   }
   delete[] send_ptr;
   delete[] send_idx;
   delete[] len;
   delete[] mpi_requests;
   delete[] mpi_stati;

   /* Now, we can create pattern for mainBlock and coupleBlock */
   Pattern_ptr main_pattern(new Pattern(MATRIX_FORMAT_DEFAULT,
               num_Pmain_cols, num_Pmain_cols, RAP_main_ptr, RAP_main_idx));
   Pattern_ptr col_couple_pattern(new Pattern(
               MATRIX_FORMAT_DEFAULT, num_Pmain_cols, num_RAPext_cols,
               RAP_couple_ptr, RAP_couple_idx));
   Pattern_ptr row_couple_pattern(new Pattern(
               MATRIX_FORMAT_DEFAULT, num_RAPext_rows, num_Pmain_cols,
               row_couple_ptr, row_couple_idx));

    /* next, create the system matrix */
    pattern.reset(new SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                  output_dist, input_dist, main_pattern, col_couple_pattern,
                  row_couple_pattern, col_connector, row_connector));
    out.reset(new SystemMatrix(A->type, pattern, row_block_size,
                               col_block_size, false, A->getRowFunctionSpace(),
                               A->getColumnFunctionSpace()));

    /* finally, fill in the data*/
    memcpy(out->mainBlock->val, RAP_main_val,
                out->mainBlock->len * sizeof(double));
    memcpy(out->col_coupleBlock->val, RAP_couple_val,
                out->col_coupleBlock->len * sizeof(double));

    delete[] RAP_main_val;
    delete[] RAP_couple_val;
    return out;
}

} // namespace paso

