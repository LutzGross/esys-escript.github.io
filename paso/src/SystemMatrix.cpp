
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
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

/* Paso: SystemMatrix */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"
#include "Options.h"
#include "PasoException.h"
#include "Preconditioner.h"
#include "Solver.h"

#include <escript/Data.h>

#include <cstring> // memcpy
#include <vector>

namespace paso {

template <>
void SystemMatrix<double>::setPreconditioner(Options* options)
{
    if (!solver_p) {
        SystemMatrix_ptr<double> mat(boost::dynamic_pointer_cast<SystemMatrix>(getPtr()));
        solver_p = Preconditioner_alloc(mat, options);
    }
}

template <>
void SystemMatrix<double>::solvePreconditioner(double* x, double* b)
{
    Preconditioner* prec=(Preconditioner*)solver_p;
    SystemMatrix_ptr<double> mat(boost::dynamic_pointer_cast<SystemMatrix>(getPtr()));
    Preconditioner_solve(prec, mat, x, b);
}

template <>
void SystemMatrix<double>::freePreconditioner()
{
    Preconditioner* prec = (Preconditioner*) solver_p;
    Preconditioner_free(prec);
    solver_p = NULL;
}

template <>
double SystemMatrix<double>::getGlobalSize() const
{
    double global_size=0;
    double my_size = mainBlock->getSize() + col_coupleBlock->getSize();
    if (mpi_info->size > 1) {
#ifdef ESYS_MPI
        MPI_Allreduce(&my_size, &global_size, 1, MPI_DOUBLE, MPI_SUM, mpi_info->comm);
#else
        global_size = my_size;
#endif
    } else {
        global_size = my_size;
    }
    return global_size;
}

template <>
index_t* SystemMatrix<double>::borrowMainDiagonalPointer() const
{
    int fail=0;
    index_t* out = mainBlock->borrowMainDiagonalPointer();
    if (out==NULL) fail=1;
#ifdef ESYS_MPI
    int fail_loc = fail;
    MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, mpi_info->comm);
#endif
    if (fail>0)
        throw PasoException("SystemMatrix::borrowMainDiagonalPointer: no main diagonal");
    return out;
}

template <>
void SystemMatrix<double>::makeZeroRowSums(double* left_over)
{
    const dim_t n = pattern->getNumOutput();
    const dim_t nblk = block_size;
    const dim_t blk = row_block_size;
    const index_t* main_ptr = borrowMainDiagonalPointer();

    rowSum(left_over);
    // left_over now holds the row sum

#pragma omp parallel for
    for (index_t ir=0; ir<n; ir++) {
        for (index_t ib=0; ib<blk; ib++) {
            const index_t irow = ib+blk*ir;
            const double rtmp2 = mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib];
            const double rtmp1 = rtmp2-left_over[irow];
            mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib] = rtmp1;
            left_over[irow]=rtmp2-rtmp1;
        }
    }
}

template <>
void SystemMatrix<double>::nullifyRows(double* mask_row, double main_diagonal_value)
{
    if (type & MATRIX_FORMAT_CSC) {
        throw PasoException("SystemMatrix::nullifyRows: Only CSR format is supported.");
    }

    if (col_block_size==1 && row_block_size==1) {
        startRowCollect(mask_row);
        mainBlock->nullifyRows_CSR_BLK1(mask_row, main_diagonal_value);
        col_coupleBlock->nullifyRows_CSR_BLK1(mask_row, 0.);
        double* remote_values = finishRowCollect();
        row_coupleBlock->nullifyRows_CSR_BLK1(remote_values, 0.);
    } else {
        startRowCollect(mask_row);
        mainBlock->nullifyRows_CSR(mask_row, main_diagonal_value);
        col_coupleBlock->nullifyRows_CSR(mask_row, 0.);
        double* remote_values = finishRowCollect();
        row_coupleBlock->nullifyRows_CSR(remote_values, 0.);
    }
}

template <>
void SystemMatrix<double>::copyColCoupleBlock()
{
    if (mpi_info->size == 1) {
        // nothing to do
        return;
    } else if (!row_coupleBlock) {
        throw PasoException("SystemMatrix::copyColCoupleBlock: "
                    "creation of row_coupleBlock pattern not supported yet.");
    } else if (row_coupler->in_use) {
        throw PasoException("SystemMatrix::copyColCoupleBlock: Coupler in use.");
    }

    const dim_t numNeighboursSend = row_coupler->connector->send->neighbour.size();
    const dim_t numNeighboursRecv = row_coupler->connector->recv->neighbour.size();
    // start receiving
    for (dim_t p = 0; p < numNeighboursRecv; p++) {
#ifdef ESYS_MPI
        const index_t irow1 = row_coupler->connector->recv->offsetInShared[p];
        const index_t irow2 = row_coupler->connector->recv->offsetInShared[p+1];
        const index_t a = row_coupleBlock->pattern->ptr[irow1];
        const index_t b = row_coupleBlock->pattern->ptr[irow2];

        MPI_Irecv(&row_coupleBlock->val[a*block_size], (b-a) * block_size,
                MPI_DOUBLE, row_coupler->connector->recv->neighbour[p],
                mpi_info->counter()+row_coupler->connector->recv->neighbour[p],
                mpi_info->comm, &row_coupler->mpi_requests[p]);

#endif
    }

    // start sending
    index_t z0 = 0;
    double* send_buffer = new double[col_coupleBlock->len];
    const size_t block_size_size = block_size*sizeof(double);

    for (dim_t p = 0; p < numNeighboursSend; p++) {
        // j_min, j_max defines the range of columns to be sent to processor p
        const index_t j_min = col_coupler->connector->recv->offsetInShared[p];
        const index_t j_max = col_coupler->connector->recv->offsetInShared[p+1];
        index_t z = z0;

        // run over the rows to be connected to processor p
        for (index_t rPtr=row_coupler->connector->send->offsetInShared[p];
                rPtr < row_coupler->connector->send->offsetInShared[p+1]; ++rPtr) {
            const index_t row = row_coupler->connector->send->shared[rPtr];

            // collect the entries in the col couple block referring to
            // columns on processor p
            for (index_t iPtr=col_coupleBlock->pattern->ptr[row];
                    iPtr < col_coupleBlock->pattern->ptr[row+1]; ++iPtr) {
                const index_t j = col_coupleBlock->pattern->index[iPtr];
                if (j_min <= j && j < j_max) {
                    memcpy(&send_buffer[z],
                           &col_coupleBlock->val[block_size*iPtr],
                           block_size_size);
                    z+=block_size;
                }
            }
        }
#ifdef ESYS_MPI
        MPI_Issend(&send_buffer[z0], z-z0, MPI_DOUBLE,
                   row_coupler->connector->send->neighbour[p],
                   mpi_info->counter()+mpi_info->rank,
                   mpi_info->comm,
                   &row_coupler->mpi_requests[p+numNeighboursRecv]);
#endif
        z0 = z;
    }

    // wait until everything is done
#ifdef ESYS_MPI
    mpi_info->incCounter(mpi_info->size);
    MPI_Waitall(numNeighboursSend+numNeighboursRecv, row_coupler->mpi_requests,
                row_coupler->mpi_stati);
#endif
    delete[] send_buffer;
}

template <>
void SystemMatrix<double>::applyBalanceInPlace(double* x, const bool RHS) const
{
    if (is_balanced) {
        if (RHS) {
            const dim_t nrow = getTotalNumRows();
#pragma omp parallel for
            for (index_t i=0; i<nrow; ++i) {
                x[i] *= balance_vector[i];
            }
        } else {
            const dim_t ncol = getTotalNumCols();
#pragma omp parallel for
            for (index_t i=0; i<ncol; ++i) {
                x[i] *= balance_vector[i];
            }
        }
    }
}

template <>
void SystemMatrix<double>::applyBalance(double* x_out, const double* x, bool RHS) const
{
    if (is_balanced) {
        if (RHS) {
            const dim_t nrow = getTotalNumRows();
#pragma omp parallel for
            for (index_t i=0; i<nrow; ++i) {
                x_out[i] = x[i] * balance_vector[i];
            }
        } else {
            const dim_t ncol = getTotalNumCols();
#pragma omp parallel for
            for (index_t i=0; i<ncol; ++i) {
                x_out[i] = x[i] * balance_vector[i];
            }
        }
    }
}

template <>
void SystemMatrix<double>::balance()
{
    const dim_t nrow = getTotalNumRows();

    if (!is_balanced) {
        if ((type & MATRIX_FORMAT_CSC) || (type & MATRIX_FORMAT_OFFSET1)) {
            throw PasoException("SystemMatrix_balance: No normalization "
                  "available for compressed sparse column or index offset 1.");
        }
        if (getGlobalNumRows() != getGlobalNumCols() ||
                row_block_size != col_block_size) {
            throw PasoException("SystemMatrix::balance: matrix needs to be a square matrix.");
        }
        // calculate absolute max value over each row
#pragma omp parallel for
        for (dim_t irow=0; irow<nrow; ++irow) {
            balance_vector[irow]=0;
        }
        mainBlock->maxAbsRow_CSR_OFFSET0(balance_vector);
        if (col_coupleBlock->pattern->ptr != NULL) {
            col_coupleBlock->maxAbsRow_CSR_OFFSET0(balance_vector);
        }

        // set balancing vector
        #pragma omp parallel for
        for (dim_t irow=0; irow<nrow; ++irow) {
            const double fac = balance_vector[irow];
            if (fac > 0) {
                balance_vector[irow]=sqrt(1./fac);
            } else {
                balance_vector[irow]=1.;
            }
        }
        ///// rescale matrix /////
        // start exchange
        startCollect(balance_vector);
        // process main block
        mainBlock->applyDiagonal_CSR_OFFSET0(balance_vector, balance_vector);
        // finish exchange
        double* remote_values = finishCollect();
        // process couple block
        if (col_coupleBlock->pattern->ptr != NULL) {
            col_coupleBlock->applyDiagonal_CSR_OFFSET0(balance_vector, remote_values);
        }
        if (row_coupleBlock->pattern->ptr != NULL) {
            row_coupleBlock->applyDiagonal_CSR_OFFSET0(remote_values, balance_vector);
        }
        is_balanced = true;
    }
}

template <>
SparseMatrix_ptr<double> SystemMatrix<double>::mergeSystemMatrix() const
{
    const index_t n = mainBlock->numRows;

    if (mpi_info->size == 1) {
        index_t* ptr = new index_t[n];
#pragma omp parallel for
        for (index_t i=0; i<n; i++)
            ptr[i] = i;
        SparseMatrix_ptr<double> out(mainBlock->getSubmatrix(n, n, ptr, ptr));
        delete[] ptr;
        return out;
    }

#ifdef ESYS_MPI
    const index_t size=mpi_info->size;
    const index_t rank=mpi_info->rank;

    // Merge main block and couple block to get the complete column entries
    // for each row allocated to current rank. Output (ptr, idx, val)
    // contains all info needed from current rank to merge a system matrix
    index_t *ptr, *idx;
    double  *val;
    mergeMainAndCouple(&ptr, &idx, &val);

    std::vector<MPI_Request> mpi_requests(size*2);
    std::vector<MPI_Status> mpi_stati(size*2);

    // Now, pass all info to rank 0 and merge them into one sparse matrix
    if (rank == 0) {
        // First, copy local ptr values into ptr_global
        const index_t global_n = getGlobalNumRows();
        index_t* ptr_global = new index_t[global_n+1];
        memcpy(ptr_global, ptr, (n+1)*sizeof(index_t));
        delete[] ptr;
        index_t iptr = n+1;
        index_t* temp_n = new index_t[size];
        index_t* temp_len = new index_t[size];
        temp_n[0] = iptr;

        // Second, receive ptr values from other ranks
        for (index_t i=1; i<size; i++) {
            const index_t remote_n = row_distribution->first_component[i+1] -
                                        row_distribution->first_component[i];
            MPI_Irecv(&ptr_global[iptr], remote_n, MPI_INT, i,
                        mpi_info->counter()+i, mpi_info->comm,
                        &mpi_requests[i]);
            temp_n[i] = remote_n;
            iptr += remote_n;
        }
        mpi_info->incCounter(size);
        MPI_Waitall(size-1, &mpi_requests[1], &mpi_stati[0]);

        // Then, prepare to receive idx and val from other ranks
        index_t len = 0;
        index_t offset = -1;
        for (index_t i=0; i<size; i++) {
            if (temp_n[i] > 0) {
                offset += temp_n[i];
                len += ptr_global[offset];
                temp_len[i] = ptr_global[offset];
            } else
                temp_len[i] = 0;
        }

        index_t* idx_global = new index_t[len];
        iptr = temp_len[0];
        offset = n+1;
        for (index_t i=1; i<size; i++) {
            len = temp_len[i];
            MPI_Irecv(&idx_global[iptr], len, MPI_INT, i,
                        mpi_info->counter()+i,
                        mpi_info->comm, &mpi_requests[i]);
            const index_t remote_n = temp_n[i];
            for (index_t j=0; j<remote_n; j++) {
                ptr_global[j+offset] = ptr_global[j+offset] + iptr;
            }
            offset += remote_n;
            iptr += len;
        }
        memcpy(idx_global, idx, temp_len[0]*sizeof(index_t));
        delete[] idx;
        MPI_Waitall(size-1, &mpi_requests[1], &mpi_stati[0]);
        mpi_info->incCounter(size);
        delete[] temp_n;

        // Then generate the sparse matrix
        const index_t rowBlockSize = mainBlock->row_block_size;
        const index_t colBlockSize = mainBlock->col_block_size;
        Pattern_ptr pat(new Pattern(mainBlock->pattern->type,
                        global_n, global_n, ptr_global, idx_global));
        SparseMatrix_ptr<double> out(new SparseMatrix<double>(mainBlock->type, pat,
                                   rowBlockSize, colBlockSize, false));

        // Finally, receive and copy the values
        iptr = temp_len[0] * block_size;
        for (index_t i=1; i<size; i++) {
            len = temp_len[i];
            MPI_Irecv(&out->val[iptr], len * block_size, MPI_DOUBLE, i,
                        mpi_info->counter()+i, mpi_info->comm,
                        &mpi_requests[i]);
            iptr += len*block_size;
        }
        memcpy(out->val, val, temp_len[0] * sizeof(double) * block_size);
        delete[] val;
        mpi_info->incCounter(size);
        MPI_Waitall(size-1, &mpi_requests[1], &mpi_stati[0]);
        delete[] temp_len;
        return out;

    } else { // it's not rank 0

        // First, send out the local ptr
        index_t tag = mpi_info->counter()+rank;
        MPI_Issend(&ptr[1], n, MPI_INT, 0, tag, mpi_info->comm,
                   &mpi_requests[0]);

        // Next, send out the local idx
        index_t len = ptr[n];
        tag += size;
        MPI_Issend(idx, len, MPI_INT, 0, tag, mpi_info->comm,
                   &mpi_requests[1]);

        // At last, send out the local val
        len *= block_size;
        tag += size;
        MPI_Issend(val, len, MPI_DOUBLE, 0, tag, mpi_info->comm,
                   &mpi_requests[2]);

        MPI_Waitall(3, &mpi_requests[0], &mpi_stati[0]);
        mpi_info->setCounter(tag + size - rank);
        delete[] ptr;
        delete[] idx;
        delete[] val;
    } // rank
#endif

    return SparseMatrix_ptr<double>();
}

template <>
SparseMatrix_ptr<cplx_t> SystemMatrix<cplx_t>::mergeSystemMatrix() const
{
    throw PasoException("SystemMatrix::mergeSystemMatrix(): complex not implemented.");
}

} // namespace paso

