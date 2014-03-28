
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

/*   Paso: SystemMatrix and SystemVector */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005,2006 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef INC_PASO_SYSTEMMATRIX
#define INC_PASO_SYSTEMMATRIX

#include "Common.h"
#include "SparseMatrix.h"
#include "SystemMatrixPattern.h"
#include "Options.h"
#include "esysUtils/Esys_MPI.h"
#include "Paso.h"
#include "Coupler.h"


typedef int Paso_SystemMatrixType;

//  this struct holds a stiffness matrix
struct Paso_SystemMatrix
{
  Paso_SystemMatrixType type;
  paso::SystemMatrixPattern *pattern;

  dim_t reference_counter;

  dim_t logical_row_block_size;
  dim_t logical_col_block_size;

  dim_t row_block_size;
  dim_t col_block_size;
  dim_t block_size;

  paso::Distribution_ptr row_distribution;
  paso::Distribution_ptr col_distribution;
  Esys_MPIInfo *mpi_info;

  paso::Coupler_ptr col_coupler;
  paso::Coupler_ptr row_coupler;

  /* this comes into play when PASO is used */
  paso::SparseMatrix* mainBlock;           /* main block */
  paso::SparseMatrix* col_coupleBlock;     /* coupling to neighbouring processors (row - col) */
  paso::SparseMatrix* row_coupleBlock;     /* coupling to neighbouring processors (col - row)  */
  paso::SparseMatrix* remote_coupleBlock;  /* coupling of rows-cols on neighbouring processors 
                                              don't assume that this is set */

  bool is_balanced;
  double *balance_vector; /* matrix may be balanced by a diagonal matrix D=diagonal(balance_vector)
                             if is_balanced is set, the matrix stored is D*A*D where A is the original matrix.
                             When the system of linear equations is solved we solve D*A*D*y=c.
                             So to solve A*x=b one needs to set c=D*b and x=D*y. */

  index_t *global_id; /* store the global ids for all cols in col_couplerBlock */

  
  
  index_t solver_package;  /* package controlling the solver pointer */
  void* solver_p;  /* pointer to data needed by a solver */

  /* this is only used for a trilinos matrix */
  void *trilinos_data; 
};

/*  interfaces: */

PASO_DLL_API
Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType,paso::SystemMatrixPattern*,dim_t,dim_t, const bool patternIsUnrolled);

PASO_DLL_API
Paso_SystemMatrix* Paso_SystemMatrix_getReference(Paso_SystemMatrix*);

PASO_DLL_API
void Paso_SystemMatrix_free(Paso_SystemMatrix*);


PASO_DLL_API
void Paso_SystemMatrix_MatrixVector(const double alpha, Paso_SystemMatrix* A, const double* in, const double beta, double* out);

PASO_DLL_API
void Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, Paso_SystemMatrix* A, const double* in, const double beta, double* out);

PASO_DLL_API
void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);

PASO_DLL_API
void Paso_SystemMatrix_applyBalanceInPlace(const Paso_SystemMatrix* A, double* x, const bool RHS);

PASO_DLL_API
void Paso_SystemMatrix_applyBalance(const Paso_SystemMatrix* A, double* x_out, const double* x, const bool RHS);

PASO_DLL_API
void Paso_SystemMatrix_balance(Paso_SystemMatrix* A);


PASO_DLL_API
void Paso_solve(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);

PASO_DLL_API
void Paso_solve_free(Paso_SystemMatrix* in);

PASO_DLL_API
void  Paso_SystemMatrix_startCollect(Paso_SystemMatrix* A,const double* in);

PASO_DLL_API
double* Paso_SystemMatrix_finishCollect(Paso_SystemMatrix* A);

PASO_DLL_API
void  Paso_SystemMatrix_startColCollect(Paso_SystemMatrix* A,const double* in);

PASO_DLL_API
double* Paso_SystemMatrix_finishColCollect(Paso_SystemMatrix* A);

PASO_DLL_API
void  Paso_SystemMatrix_startRowCollect(Paso_SystemMatrix* A,const double* in);

PASO_DLL_API
double* Paso_SystemMatrix_finishRowCollect(Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getTotalNumRows(const Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getTotalNumCols(const Paso_SystemMatrix*);

PASO_DLL_API
dim_t Paso_SystemMatrix_getGlobalNumRows(const Paso_SystemMatrix*);

PASO_DLL_API
dim_t Paso_SystemMatrix_getGlobalNumCols(const Paso_SystemMatrix*);

PASO_DLL_API
dim_t Paso_SystemMatrix_getGlobalTotalNumRows(const Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getGlobalTotalNumCols(const Paso_SystemMatrix* A);

PASO_DLL_API
double Paso_SystemMatrix_getGlobalSize(const Paso_SystemMatrix*A);

PASO_DLL_API
double Paso_SystemMatrix_getSparsity(const Paso_SystemMatrix*A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getNumRows(const Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getNumCols(const Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getRowOverlap(const Paso_SystemMatrix* A);

PASO_DLL_API
dim_t Paso_SystemMatrix_getColOverlap(const Paso_SystemMatrix* A);




PASO_DLL_API
void Paso_SystemMatrix_saveMM(Paso_SystemMatrix *, char *);

PASO_DLL_API
void Paso_SystemMatrix_saveHB(Paso_SystemMatrix *, char *);

PASO_DLL_API
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSR(char *);

PASO_DLL_API
Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSC(char *);

PASO_DLL_API
void Paso_RHS_loadMM_toCSR( char *fileName_p, double *b, dim_t size);


PASO_DLL_API
int Paso_SystemMatrix_getSystemMatrixTypeId(const index_t solver,const index_t preconditioner, const  index_t package,const  bool symmetry, Esys_MPIInfo *mpi_info);

PASO_DLL_API
dim_t Paso_SystemMatrix_getNumOutput(Paso_SystemMatrix* A);


PASO_DLL_API
void Paso_SystemMatrix_setValues(Paso_SystemMatrix*,double);

PASO_DLL_API
void Paso_SystemMatrix_add(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);

PASO_DLL_API
void Paso_SystemMatrix_rowSum(Paso_SystemMatrix* A, double* row_sum);

PASO_DLL_API
void Paso_SystemMatrix_nullifyRows(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value);


PASO_DLL_API
void Paso_SystemMatrix_makeZeroRowSums(Paso_SystemMatrix * A_p, double* left_over); 


PASO_DLL_API
void Paso_SystemMatrix_copyBlockFromMainDiagonal(Paso_SystemMatrix * A_p, double* out);

PASO_DLL_API
void Paso_SystemMatrix_copyBlockToMainDiagonal(Paso_SystemMatrix * A_p, const double* in); 

PASO_DLL_API
void Paso_SystemMatrix_copyFromMainDiagonal(Paso_SystemMatrix * A_p, double* out);

PASO_DLL_API
void Paso_SystemMatrix_copyToMainDiagonal(Paso_SystemMatrix * A_p, const double* in); 


PASO_DLL_API
void Paso_SystemMatrix_solvePreconditioner(Paso_SystemMatrix* A,double* x,double* b);

PASO_DLL_API
void Paso_SystemMatrix_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options);

PASO_DLL_API
void Paso_SystemMatrix_freePreconditioner(Paso_SystemMatrix* A);

PASO_DLL_API
void Paso_SystemMatrix_copyColCoupleBlock(Paso_SystemMatrix *A);

PASO_DLL_API
void Paso_SystemMatrix_copyRemoteCoupleBlock(Paso_SystemMatrix *A, const bool recreatePattern);

PASO_DLL_API
void Paso_SystemMatrix_fillWithGlobalCoordinates(Paso_SystemMatrix *A, const double f1);

PASO_DLL_API
void Paso_SystemMatrix_print(Paso_SystemMatrix *A);


PASO_DLL_API
void Paso_SystemMatrix_mergeMainAndCouple(Paso_SystemMatrix *A, index_t **p_ptr, index_t **p_idx, double **p_val);

PASO_DLL_API
void Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val);
void Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0_Block(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val);

PASO_DLL_API
void Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1(Paso_SystemMatrix *A, index_t **p_ptr, index_t **p_idx, double **p_val);

PASO_DLL_API
void Paso_SystemMatrix_copyMain_CSC_OFFSET1(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val);

void Paso_SystemMatrix_extendedRowsForST(Paso_SystemMatrix* A, dim_t* degree_ST, index_t* offset_ST, index_t* ST);

  
#endif /* #ifndef INC_PASO_SYSTEMMATRIX */

