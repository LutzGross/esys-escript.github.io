/* $Id:$ */


/*
********************************************************************************
*               Copyright 2006, 2007 by ACcESS MNRF                            *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/*   Paso: SparseMatrix and SystemVector */

/**************************************************************/

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SPARSEMATRIX
#define INC_PASO_SPARSEMATRIX

#include "Common.h"
#include "Pattern.h"
#include "Options.h"
#include "Paso.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */

#define MATRIX_FORMAT_DEFAULT 0
#define MATRIX_FORMAT_CSC 1
#define MATRIX_FORMAT_SYM 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8
#define MATRIX_FORMAT_TRILINOS_CRS 16

typedef int Paso_SparseMatrixType;

typedef struct Paso_SparseMatrix {
  Paso_SparseMatrixType type;
  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;
  dim_t logical_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  dim_t numRows;
  dim_t numCols;
  Paso_Pattern* pattern;
  dim_t len;
  double *val;         /* this is used for classical CSR or CSC */
} Paso_SparseMatrix;

/*  interfaces: */

Paso_SparseMatrix* Paso_SparseMatrix_alloc(Paso_SparseMatrixType,Paso_Pattern*,dim_t,dim_t);
Paso_SparseMatrix* Paso_SparseMatrix_getReference(Paso_SparseMatrix*);
void Paso_SparseMatrix_free(Paso_SparseMatrix*);
void Paso_SparseMatrix_MatrixVector_CSC_OFFSET0(double alpha, Paso_SparseMatrix* A, double* in, double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSC_OFFSET1(double alpha, Paso_SparseMatrix* A, double* in, double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(double alpha, Paso_SparseMatrix* A, double* in, double beta, double* out);
void Paso_SparseMatrix_MatrixVector_CSR_OFFSET1(double alpha, Paso_SparseMatrix* A, double* in, double beta, double* out);
void Paso_SparseMatrix_copy(Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_addAbsRow_CSR_OFFSET0(Paso_SparseMatrix*,double*);
void Paso_SparseMatrix_nullifyRowsAndCols_CSC_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSC(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRowsAndCols_CSR(Paso_SparseMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SparseMatrix_saveHB(Paso_SparseMatrix *, FILE*);
void Paso_SparseMatrix_saveMM_CSC(Paso_SparseMatrix *, FILE*);
void Paso_SparseMatrix_saveMM_CSR(Paso_SparseMatrix *, FILE*);
Paso_SparseMatrix* Paso_SparseMatrix_getSubmatrix(Paso_SparseMatrix* A,dim_t,dim_t,index_t*,index_t*);
void Paso_SparseMatrix_setValues(Paso_SparseMatrix*,double);
/*
void Paso_SparseMatrix_add(Paso_SparseMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);
Paso_SparseMatrix* Paso_SparseMatrix_loadMM_toCSR(char *);
*/

#endif /* #ifndef INC_PASO_SPARSEMATRIX */

