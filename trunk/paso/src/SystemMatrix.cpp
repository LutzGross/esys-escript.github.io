
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

/* Paso: SystemMatrix */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"
#include "Preconditioner.h"

namespace paso {

/// Allocates a SystemMatrix of given type using the given matrix pattern.
/// Values are initialized with zero. 
/// If patternIsUnrolled and type & MATRIX_FORMAT_BLK1, it is assumed
/// that the pattern is already unrolled to match the requested block size
/// and offsets. Otherwise unrolling and offset adjustment will be performed.
SystemMatrix::SystemMatrix(SystemMatrixType ntype,
                           SystemMatrixPattern_ptr npattern, int rowBlockSize,
                           int colBlockSize, bool patternIsUnrolled) :
    type(ntype),
    logical_row_block_size(rowBlockSize),
    logical_col_block_size(colBlockSize),
    is_balanced(false),
    balance_vector(NULL), 
    global_id(NULL),
    solver_package(PASO_PASO),
    solver_p(NULL),
    trilinos_data(NULL)
{
    Esys_resetError();
    if (patternIsUnrolled) {
        if (!XNOR(ntype & MATRIX_FORMAT_OFFSET1, npattern->type & MATRIX_FORMAT_OFFSET1)) {
            Esys_setError(TYPE_ERROR, "SystemMatrix: requested offset and pattern offset do not match.");
        }
    }
    // do we need to apply unrolling?
    bool unroll  
          // we don't like non-square blocks
        = (rowBlockSize != colBlockSize)
#ifndef USE_LAPACK
          // or any block size bigger than 3
          || (colBlockSize > 3) 
# endif
          // or if block size one requested and the block size is not 1
          || ((ntype & MATRIX_FORMAT_BLK1) && colBlockSize > 1)
          // or the offsets don't match
          || ((ntype & MATRIX_FORMAT_OFFSET1) != (npattern->type & MATRIX_FORMAT_OFFSET1));

    SystemMatrixType pattern_format_out = (ntype & MATRIX_FORMAT_OFFSET1)
                             ? MATRIX_FORMAT_OFFSET1 : MATRIX_FORMAT_DEFAULT;

    mpi_info = Esys_MPIInfo_getReference(npattern->mpi_info);

    if (ntype & MATRIX_FORMAT_CSC) {
        if (unroll) {
            if (patternIsUnrolled) {
                pattern=npattern;
            } else {
                pattern = npattern->unrollBlocks(pattern_format_out,
                                                 colBlockSize, rowBlockSize);
            }
            row_block_size = 1;
            col_block_size = 1;
        } else {
            pattern = npattern->unrollBlocks(pattern_format_out, 1, 1);
            row_block_size = rowBlockSize;
            col_block_size = colBlockSize;
        }
        if (Esys_noError()) {
            row_distribution = pattern->input_distribution;
            col_distribution = pattern->output_distribution;
        }
    } else {
        if (unroll) {
            if (patternIsUnrolled) {
                pattern = npattern;
            } else {
                pattern = npattern->unrollBlocks(pattern_format_out,
                                                 rowBlockSize, colBlockSize);
            }
            row_block_size = 1;
            col_block_size = 1;
        } else {
            pattern = npattern->unrollBlocks(pattern_format_out, 1, 1);
            row_block_size = rowBlockSize;
            col_block_size = colBlockSize;
        }
        if (Esys_noError()) {
            row_distribution = pattern->output_distribution;
            col_distribution = pattern->input_distribution;
        }
    }
    if (Esys_noError()) {
        if (ntype & MATRIX_FORMAT_DIAGONAL_BLOCK) {
            block_size = MIN(row_block_size, col_block_size);
        } else {
            block_size = row_block_size*col_block_size;
        }
        col_coupler.reset(new paso::Coupler(pattern->col_connector, col_block_size));
        row_coupler.reset(new paso::Coupler(pattern->row_connector, row_block_size));
        if (ntype & MATRIX_FORMAT_TRILINOS_CRS) {
#ifdef TRILINOS
            trilinos_data = Paso_TRILINOS_alloc();
#endif
        } else {
            mainBlock.reset(new paso::SparseMatrix(type, pattern->mainPattern, row_block_size, col_block_size, true));
            col_coupleBlock.reset(new paso::SparseMatrix(type, pattern->col_couplePattern, row_block_size, col_block_size, true));
            row_coupleBlock.reset(new paso::SparseMatrix(type, pattern->row_couplePattern, row_block_size, col_block_size, true));
            const dim_t n_norm = MAX(mainBlock->numCols*col_block_size, mainBlock->numRows*row_block_size);
            balance_vector = new double[n_norm];
#pragma omp parallel for
            for (dim_t i=0; i<n_norm; ++i)
                balance_vector[i] = 1.;
        }
    }
}

// deallocates a SystemMatrix
SystemMatrix::~SystemMatrix()
{
    Paso_solve_free(this);
    Esys_MPIInfo_free(mpi_info);
    delete[] balance_vector;
    delete[] global_id;
#ifdef TRILINOS
    Paso_TRILINOS_free(trilinos_data);
#endif
#ifdef Paso_TRACE
    printf("SystemMatrix: system matrix has been deallocated.\n");
#endif
}

void SystemMatrix::setPreconditioner(Paso_Options* options)
{
    if (!solver_p) {
        solver_p = Paso_Preconditioner_alloc(shared_from_this(), options);
    }
}

void SystemMatrix::solvePreconditioner(double* x, double* b)
{
    Paso_Preconditioner* prec=(Paso_Preconditioner*)solver_p;
    Paso_Preconditioner_solve(prec, shared_from_this(), x, b);
}

void SystemMatrix::freePreconditioner()
{
    Paso_Preconditioner* prec = (Paso_Preconditioner*) solver_p;
    Paso_Preconditioner_free(prec);
    solver_p = NULL;
}

double SystemMatrix::getGlobalSize() const
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

index_t* SystemMatrix::borrowMainDiagonalPointer() const
{
    int fail=0;
    index_t* out = mainBlock->borrowMainDiagonalPointer();
    if (out==NULL) fail=1;
#ifdef ESYS_MPI
    int fail_loc = fail;
    MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, mpi_info->comm);
#endif
    if (fail>0)
        Esys_setError(VALUE_ERROR, "SystemMatrix::borrowMainDiagonalPointer: no main diagonal");
    return out;
}

void SystemMatrix::makeZeroRowSums(double* left_over) 
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
            const double rtmp1 = left_over[irow];
            const double rtmp2 = mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib];
            mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib] = -rtmp1;
            left_over[irow]=rtmp2+rtmp1;
        }
    }
}

void SystemMatrix::nullifyRows(double* mask_row, double main_diagonal_value)
{
    if ((type & MATRIX_FORMAT_CSC) || (type & MATRIX_FORMAT_TRILINOS_CRS)) {
        Esys_setError(SYSTEM_ERROR,
                "SystemMatrix::nullifyRows: Only CSR format is supported.");
        return;
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

void SystemMatrix::nullifyRowsAndCols(double* mask_row, double* mask_col,
                                      double main_diagonal_value)
{
    if (type & MATRIX_FORMAT_TRILINOS_CRS) {
        Esys_setError(SYSTEM_ERROR,
               "SystemMatrix::nullifyRowsAndCols: TRILINOS is not supported.");
        return;
    }

    if (mpi_info->size > 1) {
        if (type & MATRIX_FORMAT_CSC) {
            Esys_setError(SYSTEM_ERROR, "SystemMatrix::nullifyRowsAndCols: "
                                        "CSC is not supported with MPI.");
            return;
        }

        startColCollect(mask_col);
        startRowCollect(mask_row);
        if (col_block_size==1 && row_block_size==1) {
            mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, mask_col, main_diagonal_value);
            double* remote_values = finishColCollect();
            col_coupleBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, remote_values, 0.);
            remote_values = finishRowCollect();
            row_coupleBlock->nullifyRowsAndCols_CSR_BLK1(remote_values, mask_col, 0.);
        } else {
            mainBlock->nullifyRowsAndCols_CSR(mask_row, mask_col, main_diagonal_value);
            double* remote_values = finishColCollect();
            col_coupleBlock->nullifyRowsAndCols_CSR(mask_row, remote_values, 0.);
            remote_values = finishRowCollect();
            row_coupleBlock->nullifyRowsAndCols_CSR(remote_values, mask_col, 0.); 
        } 
    } else { 
        if (col_block_size==1 && row_block_size==1) {
            if (type & MATRIX_FORMAT_CSC) {
                mainBlock->nullifyRowsAndCols_CSC_BLK1(mask_row, mask_col, main_diagonal_value);
            } else {
                mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, mask_col, main_diagonal_value);
            }
        } else {
            if (type & MATRIX_FORMAT_CSC) {
                mainBlock->nullifyRowsAndCols_CSC(mask_row, mask_col, main_diagonal_value);
            } else {
                mainBlock->nullifyRowsAndCols_CSR(mask_row, mask_col, main_diagonal_value);
            }
        }
    }
}

void SystemMatrix::copyColCoupleBlock()
{
    if (mpi_info->size == 1) {
        // nothing to do
        return;
    } else if (!row_coupleBlock) {
        Esys_setError(VALUE_ERROR, "SystemMatrix::copyColCoupleBlock: "
                    "creation of row_coupleBlock pattern not supported yet.");
        return;
    } else if (row_coupler->in_use) {
        Esys_setError(SYSTEM_ERROR,
                "SystemMatrix::copyColCoupleBlock: Coupler in use.");
        return;
    }

    // start receiving
    for (dim_t p=0; p<row_coupler->connector->recv->numNeighbors; p++) {
#ifdef ESYS_MPI
        const index_t irow1 = row_coupler->connector->recv->offsetInShared[p];
        const index_t irow2 = row_coupler->connector->recv->offsetInShared[p+1];
        const index_t a = row_coupleBlock->pattern->ptr[irow1];
        const index_t b = row_coupleBlock->pattern->ptr[irow2];
         
        MPI_Irecv(&row_coupleBlock->val[a*block_size], (b-a) * block_size,
                MPI_DOUBLE, row_coupler->connector->recv->neighbor[p],
                mpi_info->msg_tag_counter+row_coupler->connector->recv->neighbor[p],
                mpi_info->comm, &row_coupler->mpi_requests[p]);

#endif
    }

    // start sending
    index_t z0 = 0;
    double* send_buffer = new double[col_coupleBlock->len];
    const size_t block_size_size = block_size*sizeof(double);
    
    for (dim_t p=0; p<row_coupler->connector->send->numNeighbors; p++) {
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
                   row_coupler->connector->send->neighbor[p], 
                   mpi_info->msg_tag_counter+mpi_info->rank,
                   mpi_info->comm,
                   &row_coupler->mpi_requests[p+row_coupler->connector->recv->numNeighbors]);
#endif
        z0 = z;
    }

    // wait until everything is done
#ifdef ESYS_MPI
    MPI_Waitall(row_coupler->connector->send->numNeighbors+row_coupler->connector->recv->numNeighbors,
                row_coupler->mpi_requests,
                row_coupler->mpi_stati);
#endif
    ESYS_MPI_INC_COUNTER(*mpi_info, mpi_info->size);
    delete[] send_buffer;
}

void SystemMatrix::applyBalanceInPlace(double* x, const bool RHS) const
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

void SystemMatrix::applyBalance(double* x_out, const double* x, bool RHS) const
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

void SystemMatrix::balance()
{
    const dim_t nrow = getTotalNumRows();

    if (!is_balanced) {
        if ((type & MATRIX_FORMAT_CSC) || (type & MATRIX_FORMAT_OFFSET1)) {
            Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_balance: No normalization available for compressed sparse column or index offset 1.");
        }
        if (getGlobalNumRows() != getGlobalNumCols() ||
                row_block_size != col_block_size) {
            Esys_setError(SYSTEM_ERROR,"SystemMatrix::balance: matrix needs to be a square matrix.");
        }
        if (Esys_noError()) {
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
}

index_t SystemMatrix::getSystemMatrixTypeId(index_t solver,
                                            index_t preconditioner,
                                            index_t package,
                                            bool symmetry,
                                            Esys_MPIInfo* mpi_info)
{
    index_t out = -1;
    index_t true_package = Paso_Options_getPackage(solver, package, symmetry, mpi_info);

    switch(true_package) {
        case PASO_PASO:
            out = MATRIX_FORMAT_DEFAULT;
        break;

        case PASO_MKL:
            out = MATRIX_FORMAT_BLK1 | MATRIX_FORMAT_OFFSET1;
        break;

        case PASO_UMFPACK:
            if (mpi_info->size > 1) {
                Esys_setError(VALUE_ERROR, "The selected solver UMFPACK "
                        "requires CSC format which is not supported with "
                        "more than one rank.");
            } else {
                out = MATRIX_FORMAT_CSC | MATRIX_FORMAT_BLK1;
            }
        break;

        case PASO_TRILINOS:
            // Distributed CRS
            out=MATRIX_FORMAT_TRILINOS_CRS | MATRIX_FORMAT_BLK1;
        break;

        default:
            Esys_setError(VALUE_ERROR, "unknown package code");
    }
    return out;
}

} // namespace paso

