
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
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

#include "SparseMatrix.h"

namespace paso {

template <>
void SparseMatrix<double>::nullifyRows_CSR_BLK1(const double* mask_row,
                                        double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t irow=0; irow < nOut; irow++) {
        if (mask_row[irow]>0.) {
            #pragma ivdep
            for (index_t iptr=pattern->ptr[irow]-index_offset; iptr < pattern->ptr[irow+1]-index_offset; iptr++) {
                const index_t icol = pattern->index[iptr]-index_offset;
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    }
}

template <>
void SparseMatrix<double>::nullifyRows_CSR(const double* mask_row,
                                   double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t ir=0; ir < nOut; ir++) {
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

