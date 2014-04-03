
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

/* Paso: SparseMatrix

   Nullify rows and columns in the matrix

   The rows and columns are marked by positive values in
   mask_row and mask_col. Values on the main diagonal
   which are marked to set to zero by both mask_row and
   mask_col are set to main_diagonal_value
*/

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

namespace paso {

void SparseMatrix::nullifyRowsAndCols_CSC_BLK1(const double* mask_row,
                                               const double* mask_col,
                                               double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t icol=0; icol < pattern->numOutput; icol++) {
        #pragma ivdep
        for (index_t iptr=pattern->ptr[icol]-index_offset; iptr < pattern->ptr[icol+1]-index_offset; iptr++) {
            const index_t irow = pattern->index[iptr]-index_offset;
            if (mask_col[icol]>0. || mask_row[irow]>0.) {
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    }
}

void SparseMatrix::nullifyRowsAndCols_CSR_BLK1(const double* mask_row,
                                               const double* mask_col,
                                               double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t irow=0; irow < pattern->numOutput; irow++) {
        #pragma ivdep
        for (index_t iptr=pattern->ptr[irow]-index_offset; iptr < pattern->ptr[irow+1]-index_offset; iptr++) {
            const index_t icol = pattern->index[iptr]-index_offset;
            if (mask_col[icol]>0. || mask_row[irow]>0.) {
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    }
}

void SparseMatrix::nullifyRowsAndCols_CSC(const double* mask_row,
                                          const double* mask_col,
                                          double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t ic=0; ic < pattern->numOutput; ic++) {
        for (index_t iptr=pattern->ptr[ic]-index_offset; iptr < pattern->ptr[ic+1]-index_offset; iptr++) {
            for (index_t irb=0; irb < row_block_size; irb++) {
                const index_t irow=irb+row_block_size*(pattern->index[iptr]-index_offset);
                #pragma ivdep
                for (index_t icb=0; icb < col_block_size; icb++) {
                    const index_t icol=icb+col_block_size*ic;
                    if (mask_col[icol]>0. || mask_row[irow]>0.) {
                        const index_t l=iptr*block_size+irb+row_block_size*icb;
                        val[l] = (irow==icol ? main_diagonal_value : 0);
                    }
                }
            }
        }
    }
}

void SparseMatrix::nullifyRowsAndCols_CSR(const double* mask_row,
                                          const double* mask_col,
                                          double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t ir=0; ir < pattern->numOutput; ir++) {
        for (index_t iptr=pattern->ptr[ir]-index_offset; iptr < pattern->ptr[ir+1]-index_offset; iptr++) {
            for (index_t irb=0; irb < row_block_size; irb++) {
                const index_t irow=irb+row_block_size*ir;
                #pragma ivdep
                for (index_t icb=0; icb < col_block_size; icb++) {
                    const index_t icol=icb+col_block_size*(pattern->index[iptr]-index_offset);
                    if (mask_col[icol]>0. || mask_row[irow]>0.) {
                        const index_t l=iptr*block_size+irb+row_block_size*icb;
                        val[l] = (irow==icol ? main_diagonal_value : 0);
                    }
                }
            }
        }
    }
}

void SparseMatrix::nullifyRows_CSR_BLK1(const double* mask_row,
                                        double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t irow=0; irow < pattern->numOutput; irow++) {
        if (mask_row[irow]>0.) {
            #pragma ivdep
            for (index_t iptr=pattern->ptr[irow]-index_offset; iptr < pattern->ptr[irow+1]-index_offset; iptr++) {
                const index_t icol = pattern->index[iptr]-index_offset;
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    } 
}

void SparseMatrix::nullifyRows_CSR(const double* mask_row,
                                   double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
#pragma omp parallel for
    for (index_t ir=0; ir < pattern->numOutput; ir++) {
        for (index_t iptr=pattern->ptr[ir]-index_offset; iptr < pattern->ptr[ir+1]-index_offset; iptr++) {
            for (index_t irb=0; irb < row_block_size; irb++) {
                const index_t irow = irb+row_block_size*ir;
                if (mask_row[irow]>0.) {
                    #pragma ivdep
                    for (index_t icb=0; icb < col_block_size; icb++) {
                        const index_t icol=icb+col_block_size*(pattern->index[iptr]-index_offset);
                        const index_t l=iptr*block_size+irb+row_block_size*icb;
                        val[l] = (irow==icol ? main_diagonal_value : 0);
                    }
                }
            }
        }
    }
}

} // namespace paso

