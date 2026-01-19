
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


/****************************************************************************
 * Paso: SystemMatrix
 *
 *  Merge the MainBlock and CoupleBlock in the matrix
 *  Input: SystemMatrix A
 *  Output:
 *      p_ptr: the pointer to a vector of locations that start a row.
 *      p_idx: the pointer to the column indices for each of the rows,
 *             ordered by rows.
 *      p_val: the pointer to the data corresponding directly to the
 *             column entries in p_idx.
 ****************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"

namespace paso {

template <>
void SystemMatrix<double>::mergeMainAndCouple_CSR_OFFSET0_Block(index_t** p_ptr, index_t** p_idx, double** p_val) const
{
    const index_t main_num_rows = mainBlock->numRows;
    index_t* main_ptr = mainBlock->pattern->ptr;
    index_t* main_idx = mainBlock->pattern->index;
    double* main_val  = mainBlock->val;

    if (mpi_info->size == 1) {
        // allocate arrays "ptr", "index" and "val"
        const index_t num_vals = main_ptr[main_num_rows]-1;
        *p_ptr = new index_t[main_num_rows+1];
        *p_idx = new index_t[num_vals];
        *p_val = new double[num_vals * block_size];

#pragma omp parallel for
        for (index_t i=0; i<main_num_rows; i++) {
            const index_t j_lb = main_ptr[i];
            const index_t j_ub = main_ptr[i+1];
            (*p_ptr)[i] = j_lb;
            for (index_t ij_ptr=j_lb; ij_ptr<j_ub; ij_ptr++) {
                (*p_idx)[ij_ptr] = main_idx[ij_ptr];
                for (index_t ib=0; ib<block_size; ib++) {
                    (*p_val)[ij_ptr*block_size+ib] = main_val[ij_ptr*block_size+ib];
                }
            }
        }
        (*p_ptr)[main_num_rows] = main_ptr[main_num_rows];
        return;
    }

    const index_t couple_num_rows = col_coupleBlock->numRows;

    if (main_num_rows != couple_num_rows) {
        throw PasoException("SystemMatrix_mergeMainAndCouple_CSR_OFFSET0: number of rows do not match.");
    }

    double* rows = NULL;
    Coupler_ptr<real_t> coupler;
    if (global_id == NULL) {
        // prepare for global coordinates in colCoupleBlock, the results are
        // in coupler->recv_buffer
        rows = new double[main_num_rows];
        const index_t row_offset = row_distribution->getFirstComponent();
#pragma omp parallel for
        for (index_t i=0; i<main_num_rows; ++i)
            rows[i]=row_offset+i;
        coupler.reset(new Coupler<real_t>(col_coupler->connector, 1, mpi_info));
        coupler->startCollect(rows);
    }

    // initialisation, including allocate arrays "ptr", "index" and "val"
    index_t* couple_ptr = col_coupleBlock->pattern->ptr;
    index_t* couple_idx = col_coupleBlock->pattern->index;
    double*  couple_val = col_coupleBlock->val;
    const index_t col_offset = col_distribution->getFirstComponent();
    const index_t main_num_vals = main_ptr[main_num_rows]-main_ptr[0];
    const index_t couple_num_vals = couple_ptr[couple_num_rows]-couple_ptr[0];
    const index_t num_vals = main_num_vals + couple_num_vals;
    *p_ptr = new index_t[main_num_rows+1];
    *p_idx = new index_t[num_vals];
    *p_val = new double[num_vals*block_size];
    (*p_ptr)[0] = 0;

    if (global_id == NULL) {
        coupler->finishCollect();
        delete[] rows;
        const index_t num_cols = col_coupleBlock->numCols;
        global_id = new index_t[num_cols];
#pragma omp parallel for
        for (index_t i=0; i<num_cols; ++i)
            global_id[i] = coupler->recv_buffer[i];
        coupler.reset();
    }

    index_t i = 0;
    index_t j = 0;
    index_t idx = 0;
    index_t idx2 = 0;
    index_t ij_ptr = 0;

    // merge mainBlock and col_coupleBlock
    for (index_t row=1; row<=main_num_rows; row++) {
        const index_t i_ub = main_ptr[row];
        const index_t j_ub = couple_ptr[row];
        while (i < i_ub || j < j_ub) {
            ij_ptr = i + j;
            if (j < j_ub) {
                idx = global_id[couple_idx[j]];
            }
            if (i < i_ub) {
                idx2 = main_idx[i] + col_offset;
            }
            if (j == j_ub || (i < i_ub && idx2 < idx)){
                (*p_idx)[ij_ptr] = idx2;
                for (index_t ib=0; ib<block_size; ib++)
                    (*p_val)[ij_ptr*block_size+ib] = main_val[i*block_size+ib];
                i++;
            } else {
                (*p_idx)[ij_ptr] = idx;
                for (index_t ib=0; ib<block_size; ib++)
                    (*p_val)[ij_ptr*block_size+ib] = couple_val[j*block_size+ib];
                j++;
            }
        }
        (*p_ptr)[row] = ij_ptr+1;
    }
}

template <>
void SystemMatrix<double>::mergeMainAndCouple_CSR_OFFSET0(index_t** p_ptr, index_t** p_idx, double** p_val) const
{
    if (mainBlock->col_block_size != 1 || mainBlock->row_block_size != 1 ||
            col_coupleBlock->col_block_size != 1 ||
            col_coupleBlock->row_block_size != 1) {
        mergeMainAndCouple_CSR_OFFSET0_Block(p_ptr, p_idx, p_val);
        return;
    }

    const index_t main_num_rows = mainBlock->numRows;
    index_t* main_ptr = mainBlock->pattern->ptr;
    index_t* main_idx = mainBlock->pattern->index;
    double* main_val  = mainBlock->val;

    if (mpi_info->size == 1) {
        // allocate arrays "ptr", "index" and "val"
        const index_t num_vals = main_ptr[main_num_rows]-1;
        *p_ptr = new index_t[main_num_rows+1];
        *p_idx = new index_t[num_vals];
        *p_val = new double[num_vals];

#pragma omp parallel for
        for (index_t i=0; i<main_num_rows; i++) {
            const index_t j_lb = main_ptr[i];
            const index_t j_ub = main_ptr[i+1];
            (*p_ptr)[i] = j_lb;
            for (index_t ij_ptr=j_lb; ij_ptr<j_ub; ++ij_ptr) {
                (*p_idx)[ij_ptr] = main_idx[ij_ptr];
                (*p_val)[ij_ptr] = main_val[ij_ptr];
            }
        }
        (*p_ptr)[main_num_rows] = main_ptr[main_num_rows];
        return;
    }

    const index_t couple_num_rows = col_coupleBlock->numRows;

    if (main_num_rows != couple_num_rows) {
        throw PasoException("SystemMatrix::mergeMainAndCouple_CSR_OFFSET0: number of rows do not match.");
    }

    double* rows = NULL;
    Coupler_ptr<real_t> coupler;
    if (global_id == NULL) {
        // prepare for global coordinates in colCoupleBlock, the results are
        // in coupler->recv_buffer
        rows = new double[main_num_rows];
        const index_t row_offset = row_distribution->getFirstComponent();
#pragma omp parallel for
        for (index_t i=0; i<main_num_rows; ++i)
            rows[i] = row_offset+i;
        coupler.reset(new Coupler<real_t>(col_coupler->connector, 1, mpi_info));
        coupler->startCollect(rows);
    }

    // initialisation, including allocate arrays "ptr", "index" and "val"
    index_t* couple_ptr = col_coupleBlock->pattern->ptr;
    index_t* couple_idx = col_coupleBlock->pattern->index;
    double*  couple_val = col_coupleBlock->val;
    const index_t col_offset = col_distribution->getFirstComponent();
    const index_t main_num_vals = main_ptr[main_num_rows]-main_ptr[0];
    const index_t couple_num_vals = couple_ptr[couple_num_rows]-couple_ptr[0];
    const index_t num_vals = main_num_vals + couple_num_vals;
    *p_ptr = new index_t[main_num_rows+1];
    *p_idx = new index_t[num_vals];
    *p_val = new double[num_vals];
    (*p_ptr)[0] = 0;

    if (global_id == NULL) {
        coupler->finishCollect();
        delete[] rows;
        const index_t num_cols = col_coupleBlock->numCols;
        global_id = new index_t[num_cols];
#pragma omp parallel for
        for (index_t i=0; i<num_cols; ++i)
            global_id[i] = coupler->recv_buffer[i];
        coupler.reset();
    }

    index_t i = 0;
    index_t j = 0;
    index_t idx = 0;
    index_t idx2 = 0;
    index_t ij_ptr = 0;

    // merge mainBlock and col_coupleBlock
    for (index_t row=1; row<=main_num_rows; row++) {
        const index_t i_ub = main_ptr[row];
        const index_t j_ub = couple_ptr[row];
        while (i < i_ub || j < j_ub) {
            ij_ptr = i + j;
            if (j < j_ub) {
                idx = global_id[couple_idx[j]];
            }
            if (i < i_ub) {
                idx2 = main_idx[i] + col_offset;
            }
            if (j == j_ub || (i < i_ub && idx2 < idx)){
                (*p_idx)[ij_ptr] = idx2;
                (*p_val)[ij_ptr] = main_val[i];
                i++;
            } else {
                (*p_idx)[ij_ptr] = idx;
                (*p_val)[ij_ptr] = couple_val[j];
                j++;
            }
        }
        (*p_ptr)[row] = ij_ptr+1;
    }
}

template <>
void SystemMatrix<double>::mergeMainAndCouple_CSC_OFFSET1(index_t** p_ptr, index_t** p_idx, double** p_val) const
{
    throw PasoException("SystemMatrix_mergeMainAndCouple_CSC_OFFSET1: not implemented.");
}

template <>
void SystemMatrix<double>::mergeMainAndCouple(index_t** p_ptr, index_t** p_idx, double** p_val) const
{
    if (type & MATRIX_FORMAT_DEFAULT) {
        mergeMainAndCouple_CSR_OFFSET0(p_ptr, p_idx, p_val);
    } else if (type & MATRIX_FORMAT_CSC) {
        /* CSC part is for PASTIX */
        if (type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) {
            mergeMainAndCouple_CSC_OFFSET1(p_ptr, p_idx, p_val);
        } else {
            throw PasoException("SystemMatrix::mergeMainAndCouple: CSC with index 0 or block size larger than 1 is not supported.");
        }
    } else {
        throw PasoException("SystemMatrix::mergeMainAndCouple: CRS is not supported.");
    }
}

} // namespace paso

