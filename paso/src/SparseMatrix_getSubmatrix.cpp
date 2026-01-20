
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

/* Paso: SparseMatrix */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SparseMatrix.h"
#include "BlockOps.h"

namespace paso {

/*****************************************************************************

    Returns the submatrix of A where rows are gathered by index row_list
    and columns are selected by non-negative values of new_col_index.
    If new_col_index[i]>-1 new_col_index[i] gives the column of i in
    the returned submatrix.
*/


template <>
SparseMatrix_ptr<double> SparseMatrix<double>::getSubmatrix(dim_t n_row_sub, dim_t n_col_sub,
                                            const index_t* row_list,
                                            const index_t* new_col_index) const
{
    SparseMatrix_ptr<double> out;
    if (type & MATRIX_FORMAT_CSC) {
        throw PasoException("SparseMatrix::getSubmatrix: gathering submatrices supports CSR matrix format only.");
    }

    const index_t index_offset = (type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    Pattern_ptr sub_pattern(pattern->getSubpattern(n_row_sub, n_col_sub,
                                                   row_list, new_col_index));
    // create the return object
    out.reset(new SparseMatrix<double>(type, sub_pattern, row_block_size,
                               col_block_size, true));
#pragma omp parallel for
    for (int i=0; i<n_row_sub; ++i) {
        const index_t subpattern_row = row_list[i];
        for (int k=pattern->ptr[subpattern_row]-index_offset;
                k < pattern->ptr[subpattern_row+1]-index_offset; ++k) {
            index_t tmp=new_col_index[pattern->index[k]-index_offset];
            if (tmp > -1) {
                #pragma ivdep
                for (index_t m=out->pattern->ptr[i]-index_offset;
                        m < out->pattern->ptr[i+1]-index_offset; ++m) {
                    if (out->pattern->index[m]==tmp+index_offset) {
                        BlockOps_Cpy_N(block_size, &out->val[m*block_size], &val[k*block_size]);
                        break;
                    }
                }
            }
        }
    }
    return out;
}

template <>
SparseMatrix_ptr<double> SparseMatrix<double>::getBlock(int blockid) const
{
    const dim_t blocksize = row_block_size;
    const dim_t n = numRows;
    SparseMatrix_ptr<double> out(new SparseMatrix<double>(type, pattern, 1, 1, 0));

    if (blocksize==1) {
        if (blockid==1) {
#pragma omp parallel for
            for (dim_t i=0; i<n; ++i) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[iptr];
                }
            }
        } else {
            throw PasoException("SparseMatrix::getBlock: Invalid block ID requested.");
        }
    } else if (blocksize==2) {
        if (blockid==1) {
#pragma omp parallel for
            for (dim_t i=0; i<n; i++) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[4*iptr];
                }
            }
        } else if (blockid==2) {
#pragma omp parallel for
            for (dim_t i=0; i<n; i++) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[4*iptr+3];
                }
            }
        } else {
            throw PasoException("SparseMatrix::getBlock: Invalid block ID requested.");
        }
    } else if (blocksize==3) {
        if (blockid==1) {
#pragma omp parallel for
            for (dim_t i=0; i<n; i++) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[9*iptr];
                }
            }
        } else if (blockid==2) {
#pragma omp parallel for
            for (dim_t i=0; i<n; i++) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[9*iptr+4];
                }
            }
        } else if (blockid==3) {
#pragma omp parallel for
            for (dim_t i=0; i<n; i++) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[9*iptr+8];
                }
            }
        } else {
            throw PasoException("SparseMatrix::getBlock: Invalid block ID requested.");
        }
    }
    return out;
}

} // namespace paso

