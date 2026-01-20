
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


/****************************************************************************

 Paso: SystemMatrix: debugging tools

*****************************************************************************

 Author: Lutz Gross, l.gross@uq.edu.au

*****************************************************************************/

#include "SystemMatrix.h"

#include <cstring> // strcat
#include <iostream>

namespace paso {

// fills the matrix with values i+f1*j where i and j are the global row
// and column indices of the matrix entry
template <>
void SystemMatrix<double>::fillWithGlobalCoordinates(double f1)
{
    const dim_t n = getNumRows();
    const dim_t m = getNumCols();
    const index_t row_offset = row_distribution->getFirstComponent();
    const index_t col_offset = col_distribution->getFirstComponent();
    double* cols = new double[m];
    double* rows = new double[n];
    Coupler_ptr<real_t> col_couple(new Coupler<real_t>(col_coupler->connector, 1, mpi_info));
    Coupler_ptr<real_t> row_couple(new Coupler<real_t>(col_coupler->connector, 1, mpi_info));

#pragma omp parallel for
    for (dim_t i=0; i<n; ++i)
        rows[i]=row_offset+i;

    col_couple->startCollect(rows);

#pragma omp parallel for
    for (dim_t i=0; i<m; ++i)
        cols[i]=col_offset+i;

    row_couple->startCollect(cols);

    // main block
    for (dim_t q=0; q<n; ++q) {
        for (dim_t iPtr = mainBlock->pattern->ptr[q];
                iPtr < mainBlock->pattern->ptr[q+1]; ++iPtr) {
            const dim_t p = mainBlock->pattern->index[iPtr];
            for (dim_t ib=0; ib<block_size; ib++)
                mainBlock->val[iPtr*block_size+ib]=f1*rows[q]+cols[p];
        }
    }

    col_couple->finishCollect();
    if (col_coupleBlock != NULL) {
        for (dim_t q=0; q<col_coupleBlock->pattern->numOutput; ++q) {
            for (dim_t iPtr = col_coupleBlock->pattern->ptr[q];
                    iPtr < col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
                const dim_t p = col_coupleBlock->pattern->index[iPtr];
                for (dim_t ib=0; ib<block_size; ib++)
                    col_coupleBlock->val[iPtr*block_size+ib] =
                            f1*rows[q]+col_couple->recv_buffer[p];
            }
        }
    }

    row_couple->finishCollect();
    if (row_coupleBlock != NULL) {
        for (dim_t p=0; p<row_coupleBlock->pattern->numOutput; ++p) {
            for (dim_t iPtr = row_coupleBlock->pattern->ptr[p];
                    iPtr < row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                const dim_t q = row_coupleBlock->pattern->index[iPtr];
                for (dim_t ib=0; ib<block_size; ib++)
                    row_coupleBlock->val[iPtr*block_size+ib] =
                        row_couple->recv_buffer[p]*f1+cols[q];
            }
        }
    }

    delete[] cols;
    delete[] rows;
}

template <>
void SystemMatrix<double>::print() const
{
    const dim_t n = getNumRows();
    const int rank = mpi_info->rank;

    std::cerr << "rank " << rank << " Main Block:\n-----------\n";

    for (dim_t q=0; q<n; ++q) {
        std::cerr << "Row " << q << ": ";
        for (dim_t iPtr=mainBlock->pattern->ptr[q]; iPtr<mainBlock->pattern->ptr[q+1]; ++iPtr) {
            std::cerr << "(" << mainBlock->pattern->index[iPtr] << " ";
            for (dim_t ib=0; ib<block_size; ib++) {
                std::cerr << mainBlock->val[iPtr*block_size+ib] << " ";
            }
            std::cerr << "),";
        }
        std::cerr << std::endl;
    }
    if (col_coupleBlock != NULL && mpi_info->size>1) {
        std::cerr << "rank " << rank
                  << " Column Couple Block:\n------------------\n";
        for (dim_t q=0; q<col_coupleBlock->pattern->numOutput; ++q) {
            std::cerr << "Row " << q << ": ";
            for (dim_t iPtr=col_coupleBlock->pattern->ptr[q]; iPtr<col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
                if (global_id) {
                    std::cerr << "("
                        << global_id[col_coupleBlock->pattern->index[iPtr]]
                        << " " << col_coupleBlock->val[iPtr*block_size]
                        << "),";
                } else {
                    std::cerr << "("
                        << col_coupleBlock->pattern->index[iPtr]
                        << " " << col_coupleBlock->val[iPtr*block_size]
                        << "),";
                }
            }
            std::cerr << std::endl;
        }
    }
    if (row_coupleBlock != NULL && mpi_info->size>1) {
        std::cerr << "rank " << rank
                  << " Row Couple Block:\n--------------------\n";
        for (dim_t p=0; p<row_coupleBlock->pattern->numOutput; ++p) {
            std::cerr << "Row " << p << ": ";
            for (dim_t iPtr=row_coupleBlock->pattern->ptr[p]; iPtr<row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                std::cerr << "(" << row_coupleBlock->pattern->index[iPtr]
                    << " " << row_coupleBlock->val[iPtr*block_size] << "),";
            }
            std::cerr << std::endl;
        }
    }
    if (remote_coupleBlock != NULL && mpi_info->size>1) {
        std::cerr << "rank " << rank
                  << " Remote Couple Block:\n--------------------\n";
        for (dim_t p=0; p<remote_coupleBlock->pattern->numOutput; ++p) {
            std::cerr << "Row " << p << ": ";
            for (dim_t iPtr=remote_coupleBlock->pattern->ptr[p]; iPtr<remote_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                std::cerr << "(" << remote_coupleBlock->pattern->index[iPtr]
                    << " " << remote_coupleBlock->val[iPtr*block_size] << "),";
            }
            std::cerr << std::endl;
        }
    }
}

} // namespace paso

