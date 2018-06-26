
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
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

struct SparseMatrix;
typedef boost::shared_ptr<SparseMatrix> SparseMatrix_ptr;
typedef boost::shared_ptr<const SparseMatrix> const_SparseMatrix_ptr;

typedef int SparseMatrixType;

// this struct holds a sparse matrix
struct SparseMatrix : boost::enable_shared_from_this<SparseMatrix>
{
    SparseMatrix(SparseMatrixType type, Pattern_ptr pattern,
                 dim_t rowBlockSize, dim_t colBlockSize,
                 bool patternIsUnrolled);

    ~SparseMatrix();

    void setValues(double value);

    void copyFromMainDiagonal(double* out) const;

    void copyToMainDiagonal(const double* in);

    void copyBlockFromMainDiagonal(double* out) const;

    void copyBlockToMainDiagonal(const double* in);

    void applyBlockMatrix(double* block_diag, index_t* pivot, double* x,
                          const double* b) const;

    void invMain(double* inv_diag, index_t* pivot) const;

    SparseMatrix_ptr unroll(SparseMatrixType type) const;
    SparseMatrix_ptr getSubmatrix(dim_t n_row_sub,
                                  dim_t n_col_sub,
                                  const index_t* row_list,
                                  const index_t* new_col_index) const;

    SparseMatrix_ptr getBlock(int blockid) const;

    SparseMatrix_ptr getTranspose() const;

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

    static SparseMatrix_ptr loadMM_toCSR(const char* filename);


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
    double *val;

    /// package controlling the solver pointer
    index_t solver_package;

    /// pointer to data needed by a solver
    void* solver_p;
};

//  interfaces:

void SparseMatrix_MatrixVector_CSC_OFFSET0(const double alpha,
                                           const_SparseMatrix_ptr A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSC_OFFSET1(const double alpha,
                                           const_SparseMatrix_ptr A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0(const double alpha,
                                           const_SparseMatrix_ptr A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET1(const double alpha,
                                           const_SparseMatrix_ptr A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(const double alpha,
                                                const_SparseMatrix_ptr A,
                                                const double* in,
                                                const double beta, double* out);

SparseMatrix_ptr SparseMatrix_MatrixMatrix(const_SparseMatrix_ptr A,
                                           const_SparseMatrix_ptr B);

SparseMatrix_ptr SparseMatrix_MatrixMatrixTranspose(const_SparseMatrix_ptr A,
                                                    const_SparseMatrix_ptr B,
                                                    const_SparseMatrix_ptr T);

} // namespace paso

#endif // __PASO_SPARSEMATRIX_H__

