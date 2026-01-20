
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
**
*****************************************************************************/


/****************************************************************************/

/*   Paso: SparseMatrix */

/****************************************************************************/

/*   Author: lgross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SPARSEMATRIX_H__
#define __PASO_SPARSEMATRIX_H__

#include "Pattern.h"

namespace paso {

template <typename T> struct SparseMatrix;
template <typename T> using SparseMatrix_ptr = boost::shared_ptr<SparseMatrix<T> >;
template <typename T> using const_SparseMatrix_ptr = boost::shared_ptr<const SparseMatrix<T> >;

typedef int SparseMatrixType;

// this struct holds a sparse matrix
template <typename T>
struct SparseMatrix : boost::enable_shared_from_this<SparseMatrix<T> >
{
    SparseMatrix(SparseMatrixType type, Pattern_ptr pattern,
                 dim_t rowBlockSize, dim_t colBlockSize,
                 bool patternIsUnrolled);

    ~SparseMatrix();

    void setValues(T value);

    void copyFromMainDiagonal(double* out) const;

    void copyToMainDiagonal(const double* in);

    void copyBlockFromMainDiagonal(double* out) const;

    void copyBlockToMainDiagonal(const double* in);

    void applyBlockMatrix(double* block_diag, index_t* pivot, double* x,
                          const double* b) const;

    void invMain(double* inv_diag, index_t* pivot) const;

    SparseMatrix_ptr<double> unroll(SparseMatrixType type) const;
    SparseMatrix_ptr<double> getSubmatrix(dim_t n_row_sub,
                                  dim_t n_col_sub,
                                  const index_t* row_list,
                                  const index_t* new_col_index) const;

    SparseMatrix_ptr<double> getBlock(int blockid) const;

    SparseMatrix_ptr<double> getTranspose() const;

    void saveHB_CSC(const char* filename) const;

    void saveMM(const char* filename) const;

    inline index_t* borrowMainDiagonalPointer() const
    {
        return pattern->borrowMainDiagonalPointer();
    }

    inline index_t* borrowColoringPointer() const
    {
       return pattern->borrowColoringPointer();
    }

    inline dim_t getNumColors() const
    {
       return pattern->getNumColors();
    }

    inline dim_t maxDeg() const
    {
       return pattern->maxDeg();
    }

    inline dim_t getTotalNumRows() const
    {
       return numRows * row_block_size;
    }

    inline dim_t getTotalNumCols() const
    {
       return numCols * col_block_size;
    }

    inline dim_t getNumRows() const
    {
       return numRows;
    }

    inline dim_t getNumCols() const
    {
       return numCols;
    }

    inline double getSize() const
    {
        return (double)len;
    }

    inline double getSparsity() const
    {
        return getSize() / ((double)getTotalNumRows()*getTotalNumCols());
    }

    static SparseMatrix_ptr<double> loadMM_toCSR(const char* filename);


    void nullifyRowsAndCols_CSC_BLK1(const double* mask_row,
                                     const double* mask_col,
                                     double main_diagonal_value);

    void nullifyRowsAndCols_CSR_BLK1(const double* mask_row,
                                     const double* mask_col,
                                     double main_diagonal_value);

    void nullifyRowsAndCols_CSC(const double* mask_row, const double* mask_col,
                                double main_diagonal_value);

    void nullifyRowsAndCols_CSR(const double* mask_row, const double* mask_col,
                                double main_diagonal_value);

    void nullifyRows_CSR_BLK1(const double* mask_row,
                              double main_diagonal_value);

    void nullifyRows_CSR(const double* mask_row, double main_diagonal_value);

    void maxAbsRow_CSR_OFFSET0(double* array) const;

    void addAbsRow_CSR_OFFSET0(double* array) const;

    void addRow_CSR_OFFSET0(double* array) const;

    void applyDiagonal_CSR_OFFSET0(const double* left, const double* right);

    SparseMatrixType type;
    dim_t row_block_size;
    dim_t col_block_size;
    dim_t block_size;
    dim_t numRows;
    dim_t numCols;
    Pattern_ptr pattern;
    dim_t len;

    /// this is used for classical CSR or CSC
    T* val;

    /// package controlling the solver pointer
    index_t solver_package;

    /// pointer to data needed by a solver
    void* solver_p;
};

//  interfaces:

void SparseMatrix_MatrixVector_CSC_OFFSET0(const double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSC_OFFSET1(const double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0(const double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           const double beta, double* out);

template <typename T>
void SparseMatrix_MatrixVector_CSR_OFFSET1(const double alpha,
                                           const_SparseMatrix_ptr<T> A,
                                           const T* in,
                                           const double beta, T* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(const double alpha,
                                                const_SparseMatrix_ptr<double> A,
                                                const double* in,
                                                const double beta, double* out);

SparseMatrix_ptr<double> SparseMatrix_MatrixMatrix(const_SparseMatrix_ptr<double> A,
                                           const_SparseMatrix_ptr<double> B);

SparseMatrix_ptr<double> SparseMatrix_MatrixMatrixTranspose(const_SparseMatrix_ptr<double> A,
                                                    const_SparseMatrix_ptr<double> B,
                                                    const_SparseMatrix_ptr<double> T);

} // namespace paso

#include "Options.h"
#include "MKL.h"
#include "MUMPS.h"
#include "UMFPACK.h"

namespace paso {

struct Preconditioner_LocalSmoother;
void PASO_DLL_API Preconditioner_LocalSmoother_free(Preconditioner_LocalSmoother * in);

template <>
void PASO_DLL_API SparseMatrix<double>::saveHB_CSC(const char* filename) const;
template <>
void PASO_DLL_API SparseMatrix<double>::saveMM(const char* filename) const;

template <>
void PASO_DLL_API SparseMatrix<cplx_t>::saveHB_CSC(const char* filename) const;
template <>
void PASO_DLL_API SparseMatrix<cplx_t>::saveMM(const char* filename) const;

/* Allocates a SparseMatrix of given type using the given matrix pattern.
   Values are initialized with zero.
   If patternIsUnrolled and type & MATRIX_FORMAT_BLK1, it is assumed that the
   pattern is already unrolled to match the requested block size
   and offsets. Otherwise unrolling and offset adjustment will be performed.
*/
template <typename T>
SparseMatrix<T>::SparseMatrix(SparseMatrixType ntype, Pattern_ptr npattern,
                           dim_t rowBlockSize, dim_t colBlockSize,
                           bool patternIsUnrolled) :
    type(ntype),
    val(NULL),
    solver_package(PASO_PASO),
    solver_p(NULL)
{
    if (patternIsUnrolled) {
        if ((ntype & MATRIX_FORMAT_OFFSET1) != (npattern->type & MATRIX_FORMAT_OFFSET1)) {
            throw PasoException("SparseMatrix: requested offset and pattern offset do not match.");
        }
    }
    // do we need to apply unrolling?
    bool unroll
          // we don't like non-square blocks
        = (rowBlockSize != colBlockSize)
#ifndef ESYS_HAVE_LAPACK
          // or any block size bigger than 3
          || (colBlockSize > 3)
#endif
          // or if block size one requested and the block size is not 1
          || ((ntype & MATRIX_FORMAT_BLK1) && (colBlockSize > 1))
          // or if offsets don't match
          || ((ntype & MATRIX_FORMAT_OFFSET1) != (npattern->type & MATRIX_FORMAT_OFFSET1));

    SparseMatrixType pattern_format_out = (ntype & MATRIX_FORMAT_OFFSET1)
                             ? MATRIX_FORMAT_OFFSET1 : MATRIX_FORMAT_DEFAULT;

    // === compressed sparse columns ===
    if (ntype & MATRIX_FORMAT_CSC) {
        if (unroll) {
            if (patternIsUnrolled) {
                pattern = npattern;
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
        numRows = pattern->numInput;
        numCols = pattern->numOutput;
    } else {
    // === compressed sparse row ===
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
        numRows = pattern->numOutput;
        numCols = pattern->numInput;
    }
    if (ntype & MATRIX_FORMAT_DIAGONAL_BLOCK) {
        block_size = std::min(row_block_size, col_block_size);
    } else {
        block_size = row_block_size*col_block_size;
    }
    len = (size_t)(pattern->len)*(size_t)(block_size);

    val = new T[len];
    setValues(0.);
}

template <typename T>
SparseMatrix<T>::~SparseMatrix()
{
    switch (solver_package) {
        case PASO_SMOOTHER:
            Preconditioner_LocalSmoother_free((Preconditioner_LocalSmoother*) solver_p);
            break;

        case PASO_MKL:
            MKL_free(this);
            break;

        case PASO_UMFPACK:
            UMFPACK_free(this);
            break;

        case PASO_MUMPS:
            MUMPS_free(this);
            break;
    }
    delete[] val;
}

template <typename T>
void SparseMatrix<T>::setValues(T value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    if (!pattern->isEmpty()) {
        const dim_t nOut = pattern->numOutput;
#pragma omp parallel for
        for (dim_t i=0; i < nOut; ++i) {
            for (index_t iptr=pattern->ptr[i]-index_offset; iptr < pattern->ptr[i+1]-index_offset; ++iptr) {
                for (dim_t j=0; j<block_size; ++j)
                    val[iptr*block_size+j] = value;
            }
        }
    }
}

/* CSR format with offset 1 */
template <typename T>
void SparseMatrix_MatrixVector_CSR_OFFSET1(double alpha,
                                           const_SparseMatrix_ptr<T> A,
                                           const T* in,
                                           double beta, T* out)
{
    const int totalRowSize = A->numRows * A->row_block_size;
    if (std::abs(beta) > 0) {
        if (beta != 1.) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < totalRowSize; irow++) {
                out[irow] *= beta;
            }
        }
    } else {
#pragma omp parallel for schedule(static)
        for (index_t irow=0; irow < totalRowSize; irow++) {
            out[irow] = 0.;
        }
    }

    if (std::abs(alpha) > 0) {
        const int nRows = A->pattern->numOutput;
        if (A->col_block_size==1 && A->row_block_size==1) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < nRows; irow++) {
                T reg = 0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[irow]-1;
                        iptr < A->pattern->ptr[irow+1]-1; ++iptr) {
                    reg += A->val[iptr] * in[A->pattern->index[iptr]-1];
                }
                out[irow] += alpha * reg;
            }
        } else if (A->col_block_size==2 && A->row_block_size==2) {
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                T reg1 = 0.;
                T reg2 = 0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ir]-1;
                        iptr < A->pattern->ptr[ir+1]-1; iptr++) {
                    const index_t ic=2*(A->pattern->index[iptr]-1);
                    reg1 += A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic];
                    reg2 += A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic];
                }
                out[  2*ir] += alpha * reg1;
                out[1+2*ir] += alpha * reg2;
            }
        } else if (A->col_block_size==3 && A->row_block_size==3) {
#pragma omp parallel for  schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                T reg1 = 0.;
                T reg2 = 0.;
                T reg3 = 0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ir]-1;
                        iptr < A->pattern->ptr[ir+1]-1; iptr++) {
                    const index_t ic=3*(A->pattern->index[iptr]-1);
                    reg1 += A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic];
                    reg2 += A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic];
                    reg3 += A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic];
                }
                out[  3*ir] += alpha * reg1;
                out[1+3*ir] += alpha * reg2;
                out[2+3*ir] += alpha * reg3;
            }
        } else {
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                for (index_t iptr=A->pattern->ptr[ir]-1;
                        iptr < A->pattern->ptr[ir+1]-1; iptr++) {
                    for (index_t irb=0; irb < A->row_block_size; irb++) {
                        T reg = 0.;
                        #pragma ivdep
                        for (index_t icb=0; icb < A->col_block_size; icb++) {
                            const index_t icol=icb+A->col_block_size*(A->pattern->index[iptr]-1);
                            reg += A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
                        }
                        const index_t irow=irb+A->row_block_size*ir;
                        out[irow] += alpha * reg;
                    }
                }
            }
        } // blocksizes
    } // alpha > 0
}

template <typename T>
void SparseMatrix<T>::nullifyRowsAndCols_CSC_BLK1(const double* mask_row,
                                               const double* mask_col,
                                               double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t icol=0; icol < nOut; icol++) {
        #pragma ivdep
        for (index_t iptr=pattern->ptr[icol]-index_offset; iptr < pattern->ptr[icol+1]-index_offset; iptr++) {
            const index_t irow = pattern->index[iptr]-index_offset;
            if (mask_col[icol]>0. || mask_row[irow]>0.) {
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    }
}

template <typename T>
void SparseMatrix<T>::nullifyRowsAndCols_CSC(const double* mask_row,
                                          const double* mask_col,
                                          double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t ic=0; ic < nOut; ic++) {
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

template <typename T>
void SparseMatrix<T>::nullifyRowsAndCols_CSR_BLK1(const double* mask_row,
                                               const double* mask_col,
                                               double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t irow=0; irow < nOut; irow++) {
        #pragma ivdep
        for (index_t iptr=pattern->ptr[irow]-index_offset; iptr < pattern->ptr[irow+1]-index_offset; iptr++) {
            const index_t icol = pattern->index[iptr]-index_offset;
            if (mask_col[icol]>0. || mask_row[irow]>0.) {
                val[iptr] = (irow==icol ? main_diagonal_value : 0);
            }
        }
    }
}

template <typename T>
void SparseMatrix<T>::nullifyRowsAndCols_CSR(const double* mask_row,
                                          const double* mask_col,
                                          double main_diagonal_value)
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int nOut = pattern->numOutput;
#pragma omp parallel for
    for (index_t ir=0; ir < nOut; ir++) {
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

} // namespace paso

#endif // __PASO_SPARSEMATRIX_H__

