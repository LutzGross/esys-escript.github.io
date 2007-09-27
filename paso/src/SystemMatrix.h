/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/*   Paso: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005,2006 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIX
#define INC_PASO_SYSTEMMATRIX

#include "Common.h"
#include "SystemMatrixPattern.h"
#include "Options.h"
#include "Paso_MPI.h"
#include "Paso.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */

#define MATRIX_FORMAT_DEFAULT 0
#define MATRIX_FORMAT_CSC 1
#define MATRIX_FORMAT_SYM 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8
#define MATRIX_FORMAT_TRILINOS_CRS 16

typedef int Paso_SystemMatrixType;

typedef struct Paso_SystemMatrix {
  Paso_SystemMatrixType type;
  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;
  dim_t logical_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  dim_t numRows;
  dim_t myNumRows;
  dim_t myFirstRow;
  dim_t maxNumRows;
  dim_t numCols;
  dim_t myNumCols;
  dim_t myFirstCol;
  dim_t maxNumCols;

  Paso_MPIInfo *mpi_info;
  Paso_SystemMatrixPattern* pattern;
  Paso_Distribution *row_distribution;
  Paso_Distribution *col_distribution;

  dim_t myLen;
  double *val;         /* this is used for classical CSR or CSC */
  void *trilinos_data; /* this is only used for a trilinos matrix */

  double *normalizer; /* vector with a inverse of the absolute row/col sum (set by Solver.c)*/
  bool_t normalizer_is_valid;
  index_t solver_package;  /* package controling the solver pointer */
  void* solver;  /* pointer to data needed by a solver */

} Paso_SystemMatrix;

/*  interfaces: */

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType,Paso_SystemMatrixPattern*,dim_t,dim_t);
Paso_SystemMatrix* Paso_SystemMatrix_reference(Paso_SystemMatrix*);
void Paso_SystemMatrix_dealloc(Paso_SystemMatrix*);

void Paso_SystemMatrix_setValues(Paso_SystemMatrix*,double);
void Paso_SystemMatrix_copy(Paso_SystemMatrix*,double*);
void Paso_SystemMatrix_add(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);
void Paso_SystemMatrix_MatrixVector(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);
void Paso_SystemMatrix_MatrixVector_CSC_OFFSET0(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);
void Paso_SystemMatrix_MatrixVector_CSC_OFFSET1(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out, double* buffer0, double* buffer1);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_S(double alpha, Paso_SystemMatrix* A, double* in, double* out);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0_P(double alpha, Paso_SystemMatrix* A, double* in, index_t min_index, index_t max_index, double* out);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET1(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix *, char *);
void Paso_SystemMatrix_saveHB(Paso_SystemMatrix *, char *);
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSR(char *);
void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRowsAndCols_CSC_BLK1(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRowsAndCols_CSR_BLK1(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRowsAndCols_CSC(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRowsAndCols_CSR(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRows_CSR(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value);
void Paso_SystemMatrix_nullifyRows_CSR_BLK1(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value);
void Paso_SystemMatrix_nullifyCols_CSR(Paso_SystemMatrix* A, double* mask_col, double main_diagonal_value, index_t min_index, index_t max_index);
void Paso_SystemMatrix_nullifyCols_CSR_BLK1(Paso_SystemMatrix* A, double* mask_col, double main_diagonal_value, index_t min_index, index_t max_index);
void Paso_SystemMatrix_setDefaults(Paso_Options*);
int Paso_SystemMatrix_getSystemMatrixTypeId(index_t solver, index_t package, bool_t symmetry);
Paso_SystemMatrix* Paso_SystemMatrix_getSubmatrix(Paso_SystemMatrix* A,dim_t,dim_t,index_t*,index_t*);
double* Paso_SystemMatrix_borrowNormalization(Paso_SystemMatrix* A);

void Paso_solve(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);
void Paso_solve_free(Paso_SystemMatrix* in);

#endif /* #ifndef INC_PASO_SYSTEMMATRIX */

