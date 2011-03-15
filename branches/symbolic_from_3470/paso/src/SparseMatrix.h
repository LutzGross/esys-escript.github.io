
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: SparseMatrix and SystemVector */

/**************************************************************/

/*   Author: lgross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_SPARSEMATRIX
#define INC_PASO_SPARSEMATRIX

#include "Common.h"
#include "Pattern.h"
#include "Options.h"
#include "Paso.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */


typedef int Paso_SparseMatrixType;

typedef struct Paso_SparseMatrix {
  Paso_SparseMatrixType type;
  dim_t reference_counter;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  dim_t numRows;
  dim_t numCols;
  Paso_Pattern* pattern;
  dim_t len;

  double *val;         /* this is used for classical CSR or CSC */

  index_t solver_package;  /* package controling the solver pointer */
  void* solver_p;  /* pointer to data needed by a solver */

} Paso_SparseMatrix;

/*  interfaces: */

Paso_SparseMatrix* Paso_SparseMatrix_alloc(Paso_SparseMatrixType,Paso_Pattern*,dim_t,dim_t,const bool_t);
Paso_SparseMatrix* Paso_SparseMatrix_getReference(Paso_SparseMatrix*);
dim_t Paso_SparseMatrix_getNumColors(Paso_SparseMatrix*);
void Paso_SparseMatrix_applyDiagonal_CSR_OFFSET0(Paso_SparseMatrix* A, const double* left, const double* right);
index_t* Paso_SparseMatrix_borrowColoringPointer(Paso_SparseMatrix*);
void Paso_SparseMatrix_free(Paso_SparseMatrix*);
void Paso_SparseMatrix_MatrixVector_CSC_OFFSET0(const double alpha, const Paso_SparseMatrix* A, const double* in, const double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSC_OFFSET1(const double alpha, const Paso_SparseMatrix* A, const double* in, const double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(const double alpha, const Paso_SparseMatrix* A, const double* in, const double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSR_OFFSET1(const double alpha, const Paso_SparseMatrix* A, const double* in, const double beta, double* out);
void  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(const double alpha, const Paso_SparseMatrix* A, const double* in, const double beta, double* out);

void Paso_SparseMatrix_copy(Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_maxAbsRow_CSR_OFFSET0(const Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_addAbsRow_CSR_OFFSET0(const Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_addRow_CSR_OFFSET0(Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_nullifyRowsAndCols_CSC_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSC(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSR(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRows_CSR_BLK1(Paso_SparseMatrix* A, double* mask_row, double main_diagonal_value);
void Paso_SparseMatrix_saveHB_CSC(Paso_SparseMatrix *, FILE*);
Paso_SparseMatrix* Paso_SparseMatrix_getSubmatrix(Paso_SparseMatrix* A,dim_t,dim_t,index_t*,index_t*);
Paso_SparseMatrix* Paso_SparseMatrix_getBlock(Paso_SparseMatrix* A, int blockid);
Paso_SparseMatrix* Paso_SparseMatrix_MatrixMatrix(const Paso_SparseMatrix* A, const Paso_SparseMatrix* B);
void Paso_SparseMatrix_MatrixMatrix_DD(Paso_SparseMatrix *C, const Paso_SparseMatrix* A, const Paso_SparseMatrix* B);
void Paso_SparseMatrix_MatrixMatrix_DB(Paso_SparseMatrix *C, const Paso_SparseMatrix* A, const Paso_SparseMatrix* B);
void Paso_SparseMatrix_MatrixMatrix_BD(Paso_SparseMatrix *C, const Paso_SparseMatrix* A, const Paso_SparseMatrix* B);
void Paso_SparseMatrix_MatrixMatrix_BB(Paso_SparseMatrix *C, const Paso_SparseMatrix* A, const Paso_SparseMatrix* B);

Paso_SparseMatrix* Paso_SparseMatrix_unroll(const Paso_SparseMatrixType type, const Paso_SparseMatrix* A);
Paso_SparseMatrix* Paso_SparseMatrix_getTranspose(Paso_SparseMatrix* P);

void Paso_SparseMatrix_setValues(Paso_SparseMatrix*,double);
void Paso_SparseMatrix_saveMM_CSC(Paso_SparseMatrix *, FILE *);
void  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(const double alpha, const dim_t nRows, const dim_t row_block_size, const dim_t col_block_size, const index_t* ptr, const index_t* index, const double* val, const double* in, const double beta, double* out);
Paso_SparseMatrix* Paso_SparseMatrix_loadMM_toCSR( char *fileName_p );
void Paso_SparseMatrix_saveMM(Paso_SparseMatrix * A_p, char * fileName_p);
void Paso_SparseMatrix_nullifyRows_CSR(Paso_SparseMatrix*, double*, double);
index_t* Paso_SparseMatrix_borrowMainDiagonalPointer(Paso_SparseMatrix * A_p);
void Paso_SparseMatrix_copyFromMainDiagonal(Paso_SparseMatrix * A_p, double* out);
void Paso_SparseMatrix_copyToMainDiagonal(Paso_SparseMatrix * A_p, const double* in);
void Paso_SparseMatrix_copyBlockFromMainDiagonal(Paso_SparseMatrix * A_p, double* out);
void Paso_SparseMatrix_copyBlockToMainDiagonal(Paso_SparseMatrix * A_p, const double* in);
void Paso_SparseMatrix_applyBlockMatrix(Paso_SparseMatrix * A_p, double* block_diag, int* pivot, double*x, double *b);
void Paso_SparseMatrix_invMain(Paso_SparseMatrix * A_p, double* inv_diag, int* pivot);
dim_t Paso_SparseMatrix_maxDeg(Paso_SparseMatrix * A_p);
dim_t Paso_SparseMatrix_getTotalNumRows(const Paso_SparseMatrix* A);
dim_t Paso_SparseMatrix_getTotalNumCols(const Paso_SparseMatrix*A);
dim_t Paso_SparseMatrix_getNumRows(const Paso_SparseMatrix*A);
dim_t Paso_SparseMatrix_getNumCols(const Paso_SparseMatrix*A);
double Paso_SparseMatrix_getSize(const Paso_SparseMatrix*A);
double Paso_SparseMatrix_getSparsity(const Paso_SparseMatrix*A);


#endif /* #ifndef INC_PASO_SPARSEMATRIX */

