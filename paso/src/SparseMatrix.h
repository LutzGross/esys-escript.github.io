
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

/*   Paso: SparseMatrix and SystemVector */

/****************************************************************************/

/*   Author: lgross@uq.edu.au */

/****************************************************************************/

#ifndef INC_PASO_SPARSEMATRIX
#define INC_PASO_SPARSEMATRIX

#include "Common.h"
#include "Pattern.h"
#include "Options.h"
#include "Paso.h"

/****************************************************************************/

namespace paso {

typedef int SparseMatrixType;

// this struct holds a stiffness matrix
struct SparseMatrix {
    SparseMatrixType type;
    dim_t reference_counter;

    dim_t row_block_size;
    dim_t col_block_size;
    dim_t block_size;

    dim_t numRows;
    dim_t numCols;
    Pattern* pattern;
    dim_t len;

    /// this is used for classical CSR or CSC
    double *val;

    /// package controlling the solver pointer
    index_t solver_package;

    /// pointer to data needed by a solver
    void* solver_p;
};

/*  interfaces: */

SparseMatrix* SparseMatrix_alloc(SparseMatrixType type, Pattern* pattern,
                                 dim_t row_block_size, dim_t col_block_size,
                                 bool patternIsUnrolled);

SparseMatrix* SparseMatrix_getReference(SparseMatrix* mat);

dim_t SparseMatrix_getNumColors(const SparseMatrix* mat);

void SparseMatrix_applyDiagonal_CSR_OFFSET0(SparseMatrix* A,
                                            const double* left,
                                            const double* right);

index_t* SparseMatrix_borrowColoringPointer(const SparseMatrix* mat);

void SparseMatrix_free(SparseMatrix* mat);

void SparseMatrix_MatrixVector_CSC_OFFSET0(const double alpha,
                                           const SparseMatrix* A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSC_OFFSET1(const double alpha,
                                           const SparseMatrix* A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0(const double alpha,
                                           const SparseMatrix* A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET1(const double alpha,
                                           const SparseMatrix* A,
                                           const double* in,
                                           const double beta, double* out);

void SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(const double alpha,
                                                const SparseMatrix* A,
                                                const double* in,
                                                const double beta, double* out);

void SparseMatrix_maxAbsRow_CSR_OFFSET0(const SparseMatrix* A, double* array);

void SparseMatrix_addAbsRow_CSR_OFFSET0(const SparseMatrix* A, double* array);

void SparseMatrix_addRow_CSR_OFFSET0(const SparseMatrix* A, double* array);

void SparseMatrix_nullifyRowsAndCols_CSC_BLK1(SparseMatrix* A,
                                              const double* mask_row,
                                              const double* mask_col,
                                              double main_diagonal_value);

void SparseMatrix_nullifyRowsAndCols_CSR_BLK1(SparseMatrix* A,
                                              const double* mask_row,
                                              const double* mask_col,
                                              double main_diagonal_value);

void SparseMatrix_nullifyRowsAndCols_CSC(SparseMatrix* A,
                                         const double* mask_row,
                                         const double* mask_col,
                                         double main_diagonal_value);

void SparseMatrix_nullifyRowsAndCols_CSR(SparseMatrix* A,
                                         const double* mask_row,
                                         const double* mask_col,
                                         double main_diagonal_value);

void SparseMatrix_nullifyRows_CSR_BLK1(SparseMatrix* A, const double* mask_row,
                                       double main_diagonal_value);

void SparseMatrix_nullifyRows_CSR(SparseMatrix* A, const double* mask_row,
                                  double main_diagonal_value);

SparseMatrix* SparseMatrix_getSubmatrix(const SparseMatrix* A, dim_t n_row_sub,
                                        dim_t n_col_sub,
                                        const index_t* row_list,
                                        const index_t* new_col_index);

SparseMatrix* SparseMatrix_getBlock(const SparseMatrix* A, int blockid);

SparseMatrix* SparseMatrix_MatrixMatrix(const SparseMatrix* A,
                                        const SparseMatrix* B);

SparseMatrix* SparseMatrix_MatrixMatrixTranspose(const SparseMatrix* A,
                                                 const SparseMatrix* B,
                                                 const SparseMatrix* T);

SparseMatrix* SparseMatrix_unroll(SparseMatrixType type, const SparseMatrix* A);

SparseMatrix* SparseMatrix_getTranspose(const SparseMatrix* A);

void SparseMatrix_setValues(SparseMatrix* A, double value);

void SparseMatrix_saveHB_CSC(const SparseMatrix* A, FILE* handle);

void SparseMatrix_saveMM_CSC(const SparseMatrix* A, FILE* handle);

SparseMatrix* SparseMatrix_loadMM_toCSR(const char* fileName);

void SparseMatrix_saveMM(const SparseMatrix* A, const char* fileName);

index_t* SparseMatrix_borrowMainDiagonalPointer(const SparseMatrix* A);

void SparseMatrix_copyFromMainDiagonal(const SparseMatrix* A, double* out);

void SparseMatrix_copyToMainDiagonal(SparseMatrix* A, const double* in);

void SparseMatrix_copyBlockFromMainDiagonal(const SparseMatrix* A, double* out);

void SparseMatrix_copyBlockToMainDiagonal(SparseMatrix* A, const double* in);

void SparseMatrix_applyBlockMatrix(const SparseMatrix* A,
                                   double* block_diag, int* pivot,
                                   double* x, const double* b);

void SparseMatrix_invMain(const SparseMatrix* A, double* inv_diag,
                          int* pivot);

dim_t SparseMatrix_maxDeg(const SparseMatrix* A);
dim_t SparseMatrix_getTotalNumRows(const SparseMatrix* A);
dim_t SparseMatrix_getTotalNumCols(const SparseMatrix* A);
dim_t SparseMatrix_getNumRows(const SparseMatrix* A);
dim_t SparseMatrix_getNumCols(const SparseMatrix* A);
double SparseMatrix_getSize(const SparseMatrix* A);
double SparseMatrix_getSparsity(const SparseMatrix* A);

} // namespace paso

#endif // INC_PASO_SPARSEMATRIX

