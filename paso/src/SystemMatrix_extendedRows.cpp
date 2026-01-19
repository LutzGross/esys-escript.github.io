
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/
/* Paso: SystemMatrix                                       */
/*                                                          */
/*  Extend the ST sets of rows in row_coupleBlock           */
/*  Input: SystemMatrix A and ST info                       */
/*  Output:                                                 */
/*      degree_ST:                                          */
/*      offset_ST:                                          */
/*      ST:                                                 */
/****************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"
#include "PasoUtil.h"

#include <cstring> // memcpy

namespace paso {

template <>
void SystemMatrix<double>::extendedRowsForST(dim_t* degree_ST, index_t* offset_ST,
                                     index_t* ST)
{
    if (mpi_info->size == 1) return;

    // sending/receiving unknown's global ID
    index_t num_main_cols = mainBlock->numCols;
    double* cols = new double[num_main_cols];
    const index_t rank = mpi_info->rank;
    const index_t offset = col_distribution->first_component[rank];
    index_t i, j, k, p, z, z0, z1, size;

#pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_main_cols; ++i)
        cols[i] = offset + i;

    Coupler_ptr<real_t> coupler;
    if (global_id == NULL) {
        coupler.reset(new Coupler<real_t>(col_coupler->connector, 1, mpi_info));
        coupler->startCollect(cols);
    }

    const index_t my_n = mainBlock->numRows;
    double* rows = new double[my_n];
#pragma omp parallel for private(i) schedule(static)
    for (i=0; i<my_n; i++)
        rows[i] = degree_ST[i];

    index_t num_couple_cols = col_coupleBlock->numCols;
    size = num_main_cols + num_couple_cols;
    index_t overlapped_n = row_coupleBlock->numRows;
    index_t* recv_offset_ST = new index_t[overlapped_n+1];
    dim_t * recv_degree_ST = new dim_t[overlapped_n];
    index_t* send_ST = new index_t[offset_ST[my_n]];
    dim_t len = row_coupler->connector->send->numSharedComponents * size;
    index_t* send_buf = new index_t[len];

    // waiting for receiving unknown's global ID
    if (global_id == NULL) {
        coupler->finishCollect();
        global_id = new index_t[num_couple_cols];
#pragma omp parallel for private(i) schedule(static)
        for (i=0; i<num_couple_cols; ++i)
            global_id[i] = coupler->recv_buffer[i];
    }

    // sending/receiving the degree_ST
    coupler.reset(new Coupler<real_t>(row_coupler->connector, 1, mpi_info));
    coupler->startCollect(rows);

    // prepare ST with global ID
    index_t* B = new index_t[size*2];
    // find the point z in array of global_id that
    //    forall i < z, global_id[i] < offset; and
    //    forall i >= z, global_id[i] > offset + my_n
    for (i=0; i<num_couple_cols; i++)
        if (global_id[i] > offset) break;
    z = i;
    for (i=0; i<num_main_cols; i++) {
        p = offset_ST[i+1];
        z0 = 0;
        z1 = offset_ST[i];
        bool flag = (degree_ST[i] > 0 && (ST[p-1] < my_n || z == 0));
        for (j=offset_ST[i]; j<p; j++) {
            k = ST[j];
            if (!flag && k < my_n) {
                send_buf[z0] = k + offset;
                z0++;
            } else if (!flag && k - my_n >= z) {
                if (z0 >0)
                    memcpy(&(send_ST[z1]), &(send_buf[0]), z0*sizeof(index_t));
                z1 += z0;
                flag = true;
                send_ST[z1] = global_id[k - my_n];
                z1++;
            } else if (!flag && j == p-1) {
                send_ST[z1] = global_id[k - my_n];
                z1++;
                if (z0 >0)
                    memcpy(&(send_ST[z1]), &(send_buf[0]), z0*sizeof(index_t));
                flag = true;
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

    // waiting to receive the degree_ST
    coupler->finishCollect();
    delete[] rows;

    // preparing degree_ST and offset_ST for the to-be-received extended rows
#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < overlapped_n; i++)
        recv_degree_ST[i] = coupler->recv_buffer[i];
    recv_offset_ST[0] = 0;
    for (i = 0; i < overlapped_n; i++) {
        recv_offset_ST[i+1] = recv_offset_ST[i] + coupler->recv_buffer[i];
    }
    index_t* recv_ST = new index_t[recv_offset_ST[overlapped_n]];
    coupler.reset();

    // receiving ST for the extended rows
    z = 0;
    for (p=0; p<row_coupler->connector->recv->neighbour.size(); p++) {
        const index_t j_min = row_coupler->connector->recv->offsetInShared[p];
        const index_t j_max = row_coupler->connector->recv->offsetInShared[p+1];
        j = recv_offset_ST[j_max] - recv_offset_ST[j_min];
#ifdef ESYS_MPI
        MPI_Irecv(&recv_ST[z], j, MPI_INT,
                row_coupler->connector->recv->neighbour[p],
                mpi_info->counter()+row_coupler->connector->recv->neighbour[p],
                mpi_info->comm, &row_coupler->mpi_requests[p]);
#endif
        z += j;
    }

    /* sending ST for the extended rows */
    z0 = 0;
    for (p=0; p<row_coupler->connector->send->neighbour.size(); p++) {
        const index_t j_min = row_coupler->connector->send->offsetInShared[p];
        const index_t j_max = row_coupler->connector->send->offsetInShared[p+1];
        z = z0;
        for (j=j_min; j<j_max; j++) {
            const index_t row=row_coupler->connector->send->shared[j];
            if (degree_ST[row] > 0) {
                memcpy(&send_buf[z], &send_ST[offset_ST[row]], degree_ST[row] * sizeof(index_t));
                z += degree_ST[row];
            }
        }
#ifdef ESYS_MPI
        MPI_Issend(&send_buf[z0], z-z0, MPI_INT,
                 row_coupler->connector->send->neighbour[p],
                 mpi_info->counter()+mpi_info->rank, mpi_info->comm,
                 &row_coupler->mpi_requests[p+row_coupler->connector->recv->neighbour.size()]);
#endif
        z0 = z;
    }

    // first merge "cols" and "global_id" into array "B"
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

    // wait until everything is done
#ifdef ESYS_MPI
    mpi_info->incCounter(mpi_info->size);
    MPI_Waitall(row_coupler->connector->send->neighbour.size() +
                    row_coupler->connector->recv->neighbour.size(),
                    row_coupler->mpi_requests, row_coupler->mpi_stati);
#endif

    // filter the received ST (for extended rows) with cols in mainBlock as
    // well as cols in col_coupleBlock, their global ids are listed in "B"
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
        qsort(&ST[offset_ST[j]], (size_t) degree_ST[j], sizeof(index_t), util::comparIndex);
    }

    // release memory
    delete[] cols;
    delete[] B;
    delete[] send_buf;
    delete[] send_ST;
    delete[] recv_ST;
    delete[] recv_offset_ST;
    delete[] recv_degree_ST;
}

} // namespace paso

