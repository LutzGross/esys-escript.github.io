
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005,2006 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIX
#define INC_PASO_SYSTEMMATRIX

#include "Common.h"
#include "SparseMatrix.h"
#include "SystemMatrixPattern.h"
#include "Options.h"
#include "Paso_MPI.h"
#include "Paso.h"
#include "Coupler.h"


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
  Paso_SystemMatrixPattern *pattern;

  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;
  dim_t logical_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  Paso_Distribution *row_distribution;
  Paso_Distribution *col_distribution;
  Paso_MPIInfo *mpi_info;

  Paso_Coupler* col_coupler;
  Paso_Coupler* row_coupler;
  /* this comes into play when PASO is used */
  Paso_SparseMatrix* mainBlock;                      /* main block */
  Paso_SparseMatrix* col_coupleBlock;                    /* coupling to naighbouring processors (row - col) */
  Paso_SparseMatrix* row_coupleBlock;                /* coupling to naighbouring processors (col - row) (uses CSC if coupleBlock uses CSR) */
  bool_t normalizer_is_valid;
  double *normalizer; /* vector with a inverse of the absolute row/col sum (set by Solver.c)*/
  index_t solver_package;  /* package controling the solver pointer */
  void* solver;  /* pointer to data needed by a solver */

  /* this is only used for a trilinos matrix */
  void *trilinos_data; 

} Paso_SystemMatrix;

/*  interfaces: */

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType,Paso_SystemMatrixPattern*,dim_t,dim_t);
Paso_SystemMatrix* Paso_SystemMatrix_reference(Paso_SystemMatrix*);
void Paso_SystemMatrix_free(Paso_SystemMatrix*);

void Paso_SystemMatrix_MatrixVector(const double alpha, Paso_SystemMatrix* A, const double* in, const double beta, double* out);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, Paso_SystemMatrix* A, const double* in, const double beta, double* out);
void Paso_solve(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);
void Paso_solve_free(Paso_SystemMatrix* in);
void  Paso_SystemMatrix_startCollect(Paso_SystemMatrix* A,const double* in);
double* Paso_SystemMatrix_finishCollect(Paso_SystemMatrix* A);
void  Paso_SystemMatrix_startColCollect(Paso_SystemMatrix* A,const double* in);
double* Paso_SystemMatrix_finishColCollect(Paso_SystemMatrix* A);
void  Paso_SystemMatrix_startRowCollect(Paso_SystemMatrix* A,const double* in);
double* Paso_SystemMatrix_finishRowCollect(Paso_SystemMatrix* A);
void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
double* Paso_SystemMatrix_borrowNormalization(Paso_SystemMatrix* A);
dim_t Paso_SystemMatrix_getTotalNumRows(const Paso_SystemMatrix* A);
dim_t Paso_SystemMatrix_getTotalNumCols(const Paso_SystemMatrix*);
dim_t Paso_SystemMatrix_getGlobalNumRows(Paso_SystemMatrix*);
dim_t Paso_SystemMatrix_getGlobalNumCols(Paso_SystemMatrix*);

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix *, char *);
void Paso_SystemMatrix_saveHB(Paso_SystemMatrix *, char *);
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSR(char *);
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSC(char *);
void Paso_SystemMatrix_setDefaults(Paso_Options*);
int Paso_SystemMatrix_getSystemMatrixTypeId(index_t solver, index_t package, bool_t symmetry);
dim_t Paso_SystemMatrix_getNumOutput(Paso_SystemMatrix* A);
void Paso_SystemMatrix_setValues(Paso_SystemMatrix*,double);
void Paso_SystemMatrix_add(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);
void Paso_SystemMatrix_rowSum(Paso_SystemMatrix* A, double* row_sum);
void Paso_SystemMatrix_nullifyRows(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value);
void Paso_SparseMatrix_nullifyRows_CSR(Paso_SparseMatrix*, double*, double);


#endif /* #ifndef INC_PASO_SYSTEMMATRIX */

