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

#ifndef INC_PASO_SYSTEM
#define INC_PASO_SYSTEM

#include "Common.h"
#include "SystemMatrixPattern.h"
#include "Options.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */

#define MATRIX_FORMAT_DEFAULT 0
#define MATRIX_FORMAT_CSC 1
#define MATRIX_FORMAT_SYM 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8

typedef int Paso_SystemMatrixType;

typedef struct Paso_SystemMatrix {
	/*
#ifdef PASO_MPI
	Paso_CommBuffer *CommBuffer;
	Paso_MPIInfo *MPIInfo;
	dim_t numLocal;
	dim_t numInternal;
	dim_t numBoundary;
	dim_t numExternal;
	dim_t *vtxdist;
#endif
	*/
  Paso_SystemMatrixType type;
  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;
  dim_t logical_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  dim_t num_rows;
  dim_t num_cols;

  Paso_SystemMatrixPattern* pattern;

  dim_t len;
  double *val;

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
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET1(double alpha, Paso_SystemMatrix* A, double* in, double beta, double* out);

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix *, char *);
void Paso_SystemMatrix_saveHB(Paso_SystemMatrix *, char *);
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSR(char *);
void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Paso_SystemMatrix_setDefaults(Paso_Options*);
int Paso_SystemMatrix_getSystemMatrixTypeId(index_t solver, index_t package, bool_t symmetry);
Paso_SystemMatrix* Paso_SystemMatrix_getSubmatrix(Paso_SystemMatrix* A,dim_t,index_t*,index_t*);
double* Paso_SystemMatrix_borrowNormalization(Paso_SystemMatrix* A);

#endif /* #ifndef INC_PASO_SYSTEM */

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
