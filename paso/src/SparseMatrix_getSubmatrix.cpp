
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

/* Paso: SparseMatrix */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"

namespace paso {

/*****************************************************************************

    Returns the submatrix of A where rows are gathered by index row_list 
    and columns are selected by non-negative values of new_col_index.
    If new_col_index[i]>-1 new_col_index[i] gives the column of i in 
    the returned submatrix.
*/


SparseMatrix_ptr SparseMatrix::getSubmatrix(int n_row_sub, int n_col_sub,
                                            const index_t* row_list,
                                            const index_t* new_col_index) const
{
    SparseMatrix_ptr out;
    Esys_resetError();
    if (type & MATRIX_FORMAT_CSC) {
        Esys_setError(TYPE_ERROR, "SparseMatrix::getSubmatrix: gathering submatrices supports CSR matrix format only.");
        return out;
    }

    const index_t index_offset = (type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    Pattern_ptr sub_pattern(pattern->getSubpattern(n_row_sub, n_col_sub,
                                                   row_list, new_col_index));
    if (Esys_noError()) {
        // create the return object
        out.reset(new SparseMatrix(type, sub_pattern, row_block_size,
                                   col_block_size, true));
        if (Esys_noError()) {
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
                                Paso_copyShortDouble(block_size, &val[k*block_size], &out->val[m*block_size]);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return out;
}

SparseMatrix_ptr SparseMatrix::getBlock(int blockid) const
{
    const dim_t blocksize = row_block_size;
    const dim_t n = numRows;
    SparseMatrix_ptr out(new SparseMatrix(type, pattern, 1, 1, 0));

    if (blocksize==1) {
        if (blockid==1) {
#pragma omp parallel for
            for (dim_t i=0; i<n; ++i) {
                for (index_t iptr=pattern->ptr[i]; iptr<pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr] = val[iptr];
                }
            }
        } else {
            Esys_setError(VALUE_ERROR, "SparseMatrix::getBlock: Invalid block ID requested.");
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
            Esys_setError(VALUE_ERROR, "SparseMatrix::getBlock: Invalid block ID requested.");
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
            Esys_setError(VALUE_ERROR, "SparseMatrix::getBlock: Invalid block ID requested.");
        }
    }
    return out;
}

} // namespace paso

