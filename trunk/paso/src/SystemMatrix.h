
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

/*   Paso: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005,2006 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIX
#define INC_PASO_SYSTEMMATRIX

#include "Common.h"
#include "SparseMatrix.h"
#include "SystemMatrixPattern.h"
#include "Options.h"
#include "esysUtils/Esys_MPI.h"
#include "Paso.h"
#include "Coupler.h"


/**************************************************************/

/*  this struct holds a stiffness matrix: */

typedef int Paso_SystemMatrixType;

typedef struct Paso_SystemMatrix {
  Paso_SystemMatrixType type;
  Paso_SystemMatrixPattern *pattern;

  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  Paso_Distribution *row_distribution;
  Paso_Distribution *col_distribution;
  Esys_MPIInfo *mpi_info;

  Paso_Coupler* col_coupler;
  Paso_Coupler* row_coupler;

  /* this comes into play when PASO is used */
  Paso_SparseMatrix* mainBlock;                      /* main block */
  Paso_SparseMatrix* col_coupleBlock;                    /* coupling to naighbouring processors (row - col) */
  Paso_SparseMatrix* row_coupleBlock;                /* coupling to naighbouring processors (col - row)  */
  bool_t normalizer_is_valid;
  double *normalizer; /* vector with a inverse of the absolute row/col sum (set by Solver.c)*/
  index_t solver_package;  /* package controling the solver pointer */
  void* solver_p;  /* pointer to data needed by a solver */

  /* this is only used for a trilinos matrix */
  void *trilinos_data; 

} Paso_SystemMatrix;

/*  interfaces: */

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType,Paso_SystemMatrixPattern*,dim_t,dim_t, const bool_t patternIsUnrolled);
Paso_SystemMatrix* Paso_SystemMatrix_getReference(Paso_SystemMatrix*);
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
void Paso_RHS_loadMM_toCSR( char *fileName_p, double *b, dim_t size);
void Paso_SystemMatrix_setDefaults(Paso_Options*);
int Paso_SystemMatrix_getSystemMatrixTypeId(const index_t solver,const index_t preconditioner, const  index_t package,const  bool_t symmetry, Esys_MPIInfo *mpi_info);
dim_t Paso_SystemMatrix_getNumOutput(Paso_SystemMatrix* A);
void Paso_SystemMatrix_setValues(Paso_SystemMatrix*,double);
void Paso_SystemMatrix_add(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);
void Paso_SystemMatrix_rowSum(Paso_SystemMatrix* A, double* row_sum);
void Paso_SystemMatrix_nullifyRows(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value);

void Paso_SystemMatrix_makeZeroRowSums(Paso_SystemMatrix * A_p, double* left_over); 
void Paso_SystemMatrix_copyBlockFromMainDiagonal(Paso_SystemMatrix * A_p, double* out);
void Paso_SystemMatrix_copyBlockToMainDiagonal(Paso_SystemMatrix * A_p, const double* in); 
void Paso_SystemMatrix_copyFromMainDiagonal(Paso_SystemMatrix * A_p, double* out);
void Paso_SystemMatrix_copyToMainDiagonal(Paso_SystemMatrix * A_p, const double* in); 

void Paso_SystemMatrix_solvePreconditioner(Paso_SystemMatrix* A,double* x,double* b);
void Paso_SystemMatrix_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_SystemMatrix_freePreconditioner(Paso_SystemMatrix* A);

  
#endif /* #ifndef INC_PASO_SYSTEMMATRIX */

