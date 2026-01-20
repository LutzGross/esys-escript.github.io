
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

/* Paso: SparseMatrix */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SparseMatrix.h"
#include "BlockOps.h"
#include "MKL.h"
#include "Options.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "UMFPACK.h"
#include "MUMPS.h"
#include "mmio.h"

#include <boost/scoped_array.hpp>
#include <fstream>

/****************************************************************************/

namespace paso {

using escript::IndexList;

/* debug: print the entries */
/*
void print_entries(index_t *r, index_t *c, double *v, int nz)
{
    for(int i=0; i<nz; i++)
        printf("(%ld, %ld) == %e\n", (long)r[i], (long)c[i], v[i]);
}
*/

/* swap function */
void swap(index_t *r, index_t *c, double *v, int left, int right)
{
    double v_temp;
    index_t temp;

    temp = r[left];
    r[left] = r[right];
    r[right] = temp;

    temp = c[left];
    c[left] = c[right];
    c[right] = temp;

    v_temp = v[left];
    v[left] = v[right];
    v[right] = v_temp;
}

void q_sort(index_t *row, index_t *col, double *val, int begin, int end, int N)
{
    int l, r;
    index_t pivot, lval;

    if (end > begin) {
        pivot = N * row[begin] + col[begin];
        l = begin + 1;
        r = end;

        while (l < r) {
            lval = N * row[l] + col[l];
            if (lval < pivot)
                l++;
            else {
                r--;
                swap( row, col, val, l, r );
            }
        }
        l--;
        swap(row, col, val, begin, l);
        q_sort(row, col, val, begin, l, N);
        q_sort(row, col, val, r, end, N);
    }
}


template <>
SparseMatrix_ptr<double> SparseMatrix<double>::loadMM_toCSR(const char* filename)
{
    SparseMatrix_ptr<double> out;
    int i;
    MM_typecode matrixCode;

    // open the file
    std::ifstream f(filename);
    if (f.fail()) {
        throw PasoException("SparseMatrix::loadMM_toCSR: Cannot open file for reading.");
    }

    // process banner
    if (mm_read_banner(f, &matrixCode) != 0) {
        f.close();
        throw PasoException("SparseMatrix::loadMM_toCSR: Error processing MM banner.");
    }
    if (!(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode))) {
        f.close();
        throw PasoException("SparseMatrix::loadMM_toCSR: found Matrix Market type is not supported.");
    }

    // get matrix size
    int M, N, nz;

    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        f.close();
        throw PasoException("SparseMatrix::loadMM_toCSR: Could not parse matrix size.");
    }

    // prepare storage
    index_t* col_ind = new index_t[nz];
    index_t* row_ind = new index_t[nz];
    index_t* row_ptr = new index_t[M+1];
    double* val = new double[nz];

    // perform actual read of elements
    for (i=0; i<nz; i++) {
        f >> row_ind[i] >> col_ind[i] >> val[i];
        //scan_ret = fscanf(fileHandle_p, "%d %d %le\n", &row_ind[i], &col_ind[i], &val[i]);
        if (!f.good()) {
            delete[] val;
            delete[] row_ind;
            delete[] col_ind;
            delete[] row_ptr;
            f.close();
            return out;
        }
        row_ind[i]--;
        col_ind[i]--;
    }
    f.close();

    // sort the entries
    q_sort(row_ind, col_ind, val, 0, nz, N);

    // setup row_ptr
    int curr_row = 0;
    for(i=0; (i<nz && curr_row<M); curr_row++) {
        while(row_ind[i] != curr_row)
            i++;
        row_ptr[curr_row] = i;
    }
    row_ptr[M] = nz;

    Pattern_ptr mainPattern(new Pattern(MATRIX_FORMAT_DEFAULT, M, N,
                                        row_ptr, col_ind));
    out.reset(new SparseMatrix(MATRIX_FORMAT_DEFAULT, mainPattern, 1, 1, true));

    // copy values
    for (i=0; i<nz; i++)
        out->val[i] = val[i];

    delete[] val;
    delete[] row_ind;
    return out;
}

template <>
void SparseMatrix<double>::saveMM(const char* filename) const
{
    if (col_block_size != row_block_size) {
        throw PasoException("SparseMatrix::saveMM: currently only square blocks are supported.");
    }

    // open the file
    std::ofstream f(filename);
    if (f.fail()) {
        throw PasoException("SparseMatrix::saveMM: File could not be opened for writing");
    }
    if (type & MATRIX_FORMAT_CSC) {
        throw PasoException("SparseMatrix::saveMM does not support CSC.");
    } else {
        MM_typecode matcode;
        mm_initialize_typecode(&matcode);
        mm_set_matrix(&matcode);
        mm_set_coordinate(&matcode);
        mm_set_real(&matcode);

        const dim_t N = getNumRows();
        const dim_t M = getNumCols();
        mm_write_banner(f, matcode);
        mm_write_mtx_crd_size(f, N*row_block_size,
                              M*col_block_size, pattern->ptr[N]*block_size);

        const index_t offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);

        f.precision(15);

        if (type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
            for (dim_t i=0; i<N; i++) {
                for (dim_t iptr = pattern->ptr[i]-offset; iptr<pattern->ptr[i+1]-offset; ++iptr) {
                    const dim_t j=pattern->index[iptr]-offset;
                    for (dim_t ib=0; ib<block_size; ib++) {
                        const dim_t irow=ib+row_block_size*i;
                        const dim_t icol=ib+col_block_size*j;
                        f << irow+1 << " " << icol+1 << " "
                          << val[iptr*block_size+ib] << std::endl;
                    }
                }
            }
        } else {
            for (dim_t i=0; i<N; i++) {
                for (dim_t iptr = pattern->ptr[i]-offset; iptr<pattern->ptr[i+1]-offset; ++iptr) {
                    const dim_t j=pattern->index[iptr]-offset;
                    for (dim_t irb=0; irb<row_block_size; irb++) {
                        const dim_t irow=irb+row_block_size*i;
                        for (dim_t icb=0; icb<col_block_size; icb++) {
                            const dim_t icol=icb+col_block_size*j;
                            f << irow+1 << " " << icol+1 << " "
                                << val[iptr*block_size+irb+row_block_size*icb]
                                << std::endl;
                        }
                    }
                }
            }
        }
    }
    // close the file
    f.close();
}

template <>
void SparseMatrix<double>::addAbsRow_CSR_OFFSET0(double* array) const
{
    const dim_t nOut = pattern->numOutput;
#pragma omp parallel for
    for (dim_t ir=0; ir < nOut; ir++) {
        for (dim_t irb=0; irb < row_block_size; irb++) {
            const dim_t irow = irb+row_block_size*ir;
            double fac=0.;
            for (index_t iptr=pattern->ptr[ir]; iptr < pattern->ptr[ir+1]; iptr++) {
                for (dim_t icb=0; icb < col_block_size; icb++) {
                    const index_t idx = iptr*block_size+irb+row_block_size*icb;
                    fac += std::abs(val[idx]);
                }
            }
            array[irow]+=fac;
        }
    }
}

template <>
void SparseMatrix<double>::maxAbsRow_CSR_OFFSET0(double* array) const
{
    const dim_t nOut = pattern->numOutput;
#pragma omp parallel for
    for (dim_t ir=0; ir < nOut; ir++) {
        for (dim_t irb=0; irb < row_block_size; irb++) {
            const dim_t irow = irb+row_block_size*ir;
            double fac=0.;
            for (index_t iptr=pattern->ptr[ir]; iptr < pattern->ptr[ir+1]; iptr++) {
                for (dim_t icb=0; icb < col_block_size; icb++) {
                    const index_t idx = iptr*block_size+irb+row_block_size*icb;
                    fac=std::max(fac, std::abs(val[idx]));
                }
            }
            array[irow]=std::max(array[irow], fac);
        }
    }
}

template <>
void SparseMatrix<double>::addRow_CSR_OFFSET0(double* array) const
{
    const dim_t nOut = pattern->numOutput;
#pragma omp parallel for
    for (dim_t ir=0; ir < nOut; ir++) {
        for (dim_t irb=0; irb < row_block_size; irb++) {
            dim_t irow=irb+row_block_size*ir;
            double fac=0.;
            for (index_t iptr=pattern->ptr[ir]; iptr<pattern->ptr[ir+1]; iptr++) {
                for (dim_t icb=0; icb < col_block_size; icb++)
                    fac += val[iptr*block_size+irb+row_block_size*icb];

            }
            array[irow]+=fac;
        }
    }
}

template <>
void SparseMatrix<double>::copyBlockToMainDiagonal(const double* in)
{
    const dim_t n = pattern->numOutput;
    const dim_t nblk = block_size;
    const size_t nblk_size = sizeof(double)*nblk;
    const index_t* main_ptr = borrowMainDiagonalPointer();
#pragma omp parallel for
    for (index_t ir=0; ir < n;ir++) {
        memcpy((void*)&val[main_ptr[ir]*nblk], (void*)&in[nblk*ir], nblk_size);
    }
}

template <>
void SparseMatrix<double>::copyBlockFromMainDiagonal(double* out) const
{
    const dim_t n = pattern->numOutput;
    const dim_t nblk = block_size;
    const size_t nblk_size = sizeof(double)*nblk;
    const index_t* main_ptr = borrowMainDiagonalPointer();
#pragma omp parallel for
    for (index_t ir=0; ir < n; ir++) {
        memcpy((void*)&out[nblk*ir], (void*)&val[main_ptr[ir]*nblk], nblk_size);
    }
}

template <>
void SparseMatrix<double>::copyFromMainDiagonal(double* out) const
{
    const dim_t n = pattern->numOutput;
    const dim_t nblk = block_size;
    const dim_t blk = std::min(row_block_size, col_block_size);
    const index_t* main_ptr = borrowMainDiagonalPointer();
#pragma omp parallel for
    for (index_t ir=0; ir < n; ir++) {
        for (index_t ib=0; ib < blk; ib++) {
            out[ir*blk+ib] = val[main_ptr[ir]*nblk+ib+row_block_size*ib];
        }
    }
}

template <>
void SparseMatrix<double>::copyToMainDiagonal(const double* in)
{
    const dim_t n = pattern->numOutput;
    const dim_t nblk = block_size;
    const dim_t blk = std::min(row_block_size, col_block_size);
    const index_t* main_ptr = borrowMainDiagonalPointer();
#pragma omp parallel for
    for (index_t ir=0; ir < n; ir++) {
        for (index_t ib=0; ib < blk; ib++) {
            val[main_ptr[ir]*nblk+ib+row_block_size*ib] = in[ir*blk+ib];
        }
    }
}

template <>
void SparseMatrix<double>::applyDiagonal_CSR_OFFSET0(const double* left,
                                             const double* right)
{
    const dim_t row_block = row_block_size;
    const dim_t col_block = col_block_size;
    const dim_t n_block = row_block*col_block;
    const dim_t nOut = pattern->numOutput;

#pragma omp parallel for
    for (index_t ir=0; ir < nOut; ir++) {
        for (index_t irb=0; irb < row_block; irb++) {
            const index_t irow = irb+row_block*ir;
            const double rtmp = left[irow];
            for (index_t iptr=pattern->ptr[ir]; iptr < pattern->ptr[ir+1]; iptr++) {
                #pragma ivdep
                for (index_t icb=0; icb < col_block_size; icb++) {
                    const index_t icol = icb+col_block*pattern->index[iptr];
                    const index_t l = iptr*n_block + irb+row_block*icb;
                    val[l] *= rtmp*right[icol];
                }
            }
        }
    }
}

template <>
void SparseMatrix<double>::invMain(double* inv_diag, index_t* pivot) const
{
    int failed = 0;
    double A11;
    const dim_t n=numRows;
    const dim_t n_block=row_block_size;
    const dim_t m_block=col_block_size;
    dim_t i;
    index_t iPtr;
    index_t* main_ptr=pattern->borrowMainDiagonalPointer();
    // check matrix is square
    if (m_block != n_block) {
        throw PasoException("SparseMatrix::invMain: square block size expected.");
    }
    if (n_block == 1) {
#pragma omp parallel for private(i, iPtr, A11) schedule(static)
        for (i = 0; i < n; i++) {
            iPtr = main_ptr[i];
            A11 = val[iPtr];
            if (std::abs(A11) > 0.) {
                inv_diag[i]=1./A11;
            } else {
                failed=1;
            }
        }
    } else if (n_block==2) {
#pragma omp parallel for private(i, iPtr) schedule(static)
        for (i = 0; i < n; i++) {
            iPtr = main_ptr[i];
            BlockOps_invM_2(&inv_diag[i*4], &val[iPtr*4], &failed);
        }
    } else if (n_block==3) {
#pragma omp parallel for private(i, iPtr) schedule(static)
        for (i = 0; i < n; i++) {
            iPtr = main_ptr[i];
            BlockOps_invM_3(&inv_diag[i*9], &val[iPtr*9], &failed);
        }
    } else {
#pragma omp parallel for private(i, iPtr) schedule(static)
        for (i = 0; i < n; i++) {
            iPtr = main_ptr[i];
            BlockOps_Cpy_N(block_size, &inv_diag[i*block_size], &val[iPtr*block_size]);
            BlockOps_invM_N(n_block, &inv_diag[i*block_size], &pivot[i*n_block], &failed);
        }
    }
    if (failed > 0) {
        throw PasoException("SparseMatrix::invMain: non-regular main diagonal block.");
    }
}

template <>
void SparseMatrix<double>::applyBlockMatrix(double* block_diag, index_t* pivot,
                                    double* x, const double *b) const
{
    const dim_t n = numRows;
    const dim_t n_block = row_block_size;
    util::copy(n_block*n, x, b);
    BlockOps_solveAll(n_block, n, block_diag, pivot, x);
}

template <>
SparseMatrix_ptr<double> SparseMatrix<double>::getTranspose() const
{
    const dim_t m = numCols;
    const dim_t n = numRows;
    boost::scoped_array<IndexList> index_list(new IndexList[m]);

    for (dim_t i=0; i<n; ++i) {
        for (index_t iptr2=pattern->ptr[i]; iptr2<pattern->ptr[i+1]; ++iptr2) {
            const index_t j = pattern->index[iptr2];
            index_list[j].insertIndex(i);
        }
    }

    Pattern_ptr ATpattern(Pattern::fromIndexListArray(0,m,index_list.get(),0,n,0));
    SparseMatrix_ptr<double> AT(new SparseMatrix(type, ATpattern, col_block_size, row_block_size, false));

    if ( ((type & MATRIX_FORMAT_DIAGONAL_BLOCK) && (block_size == 1)) ||
         (row_block_size == 1 && col_block_size == 1)) {
#pragma omp parallel for
        for (dim_t i=0; i<m; ++i) {
            for (index_t iptr_AT=AT->pattern->ptr[i]; iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
                const index_t j = AT->pattern->index[iptr_AT];
                index_t jptr_A = pattern->ptr[j];
                const index_t* start_p = &pattern->index[jptr_A];
                const index_t* where_p=(index_t*)bsearch(&i, start_p,
                                          pattern->ptr[j+1]-jptr_A,
                                          sizeof(index_t), util::comparIndex);
                if (where_p != NULL) { // this should always be the case
                    jptr_A += (index_t)(where_p-start_p);
                    AT->val[iptr_AT] = val[jptr_A];
                }
            }
        }
    } else {
        if (type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
#pragma omp parallel for
            for (dim_t i=0; i<m; ++i) {
                for (index_t iptr_AT=AT->pattern->ptr[i]; iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
                    const index_t j = AT->pattern->index[iptr_AT];
                    index_t jptr_A = pattern->ptr[j];
                    const index_t* start_p = &pattern->index[jptr_A];
                    const index_t* where_p = (index_t*)bsearch(&i, start_p,
                                         pattern->ptr[j+1]-jptr_A,
                                         sizeof(index_t), util::comparIndex);
                    if (where_p != NULL) { // this should always be the case
                        jptr_A += (index_t)(where_p-start_p);
                        for (dim_t ib=0; ib < block_size; ++ib)
                            AT->val[iptr_AT*block_size+ib] = val[jptr_A*block_size+ib];
                    }
                }
            }
        } else {
#pragma omp parallel for
            for (dim_t i=0; i<m; ++i) {
                for (index_t iptr_AT=AT->pattern->ptr[i]; iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
                    const index_t j = AT->pattern->index[iptr_AT];
                    index_t jptr_A = pattern->ptr[j];
                    const index_t* start_p = &pattern->index[jptr_A];
                    const index_t* where_p=(index_t*)bsearch(&i, start_p,
                                       pattern->ptr[j + 1]-jptr_A,
                                       sizeof(index_t), util::comparIndex);
                    if (where_p != NULL) { // this should always be the case
                        jptr_A += (index_t)(where_p-start_p);
                        for (index_t irb=0; irb < row_block_size; ++irb) {
                            for (index_t icb=0 ; icb < col_block_size; ++icb) {
                                AT->val[iptr_AT*block_size+icb+col_block_size*irb] = val[jptr_A*block_size+irb+row_block_size*icb];
                            }
                        }
                    }
                }
            }
        }
    }
    return AT;
}

template <>
SparseMatrix_ptr<double> SparseMatrix<double>::unroll(SparseMatrixType newType) const
{
    const index_t out_type = (newType & MATRIX_FORMAT_BLK1) ? newType : newType + MATRIX_FORMAT_BLK1;
    SparseMatrix_ptr<double> out(new SparseMatrix(out_type, pattern, row_block_size, col_block_size, false));

    const dim_t n = numRows;
    const index_t A_offset = (type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    const index_t out_offset = (out_type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);

    if (out->type & MATRIX_FORMAT_CSC) {
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i) {
            for (index_t iptr=pattern->ptr[i]-A_offset; iptr<pattern->ptr[i+1]-A_offset; ++iptr) {
                const index_t j = pattern->index[iptr]-A_offset;
                for (dim_t icb=0; icb<col_block_size; ++icb) {
                    const index_t icol=j*col_block_size+icb;
                    const index_t* start_p=&out->pattern->index[out->pattern->ptr[icol]-out_offset];
                    const index_t l_col=out->pattern->ptr[icol+1]-out->pattern->ptr[icol];
                    for (dim_t irb=0; irb<row_block_size; ++irb) {
                        const index_t irow=row_block_size*i+irb+out_offset;
                        const index_t* where_p = (index_t*)bsearch(&irow,
                                    start_p, l_col, sizeof(index_t),
                                    util::comparIndex);
                        if (where_p != NULL)
                            out->val[out->pattern->ptr[icol]-out_offset+(index_t)(where_p-start_p)] =
                                val[block_size*iptr+irb+row_block_size*icb];
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i) {
            for (index_t iptr=pattern->ptr[i]-A_offset; iptr<pattern->ptr[i+1]-A_offset; ++iptr) {
                const index_t j = pattern->index[iptr]-A_offset;
                for (dim_t irb=0; irb<row_block_size; ++irb) {
                    const index_t irow=row_block_size*i+irb;
                    const index_t* start_p = &out->pattern->index[out->pattern->ptr[irow]-out_offset];
                    const index_t l_row=out->pattern->ptr[irow+1]-out->pattern->ptr[irow];
                    for (dim_t icb=0; icb<col_block_size; ++icb) {
                        const index_t icol=j*col_block_size+icb+out_offset;
                        const index_t* where_p = (index_t*)bsearch(&icol,
                                    start_p, l_row, sizeof(index_t),
                                    util::comparIndex);
                        if (where_p != NULL)
                            out->val[out->pattern->ptr[irow]-out_offset+(index_t)(where_p-start_p)] =
                                val[block_size*iptr+irb+row_block_size*icb];
                    }
                }
            }
        }
    }
    return out;
}

template <>
void SparseMatrix<cplx_t>::saveMM(const char* filename) const
{
    throw PasoException("SparseMatrix::saveMM(): complex not implemented.");
}

} // namespace paso

