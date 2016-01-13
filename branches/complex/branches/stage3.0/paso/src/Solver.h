
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_SOLVER
#define INC_SOLVER

#include "SystemMatrix.h"
#include "performance.h"
#include "Functions.h"

#define PASO_TRACE
/* error codes used in the solver */
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10
#define SOLVER_NEGATIVE_NORM_ERROR -11

#define TOLERANCE_FOR_SCALARS (double)(0.)
#define PASO_ONE (double)(1.0)
#define PASO_ZERO (double)(0.0)

/* static double ONE=1.; */
/* static double ZERO=0.;*/
/*static double TOLERANCE_FOR_SCALARS=0.;*/

/* jacobi  preconditioner */

typedef struct Paso_Solver_Jacobi {
  dim_t n_block;
  dim_t n;
  double* values;
  index_t* pivot;
} Paso_Solver_Jacobi;


/* ILU preconditioner */
struct Paso_Solver_ILU {
  dim_t n_block;
  dim_t n;
  index_t num_colors;
  index_t* colorOf;
  index_t* main_iptr;
  double* factors;
  Paso_Pattern* pattern;
};
typedef struct Paso_Solver_ILU Paso_Solver_ILU;

/* GS preconditioner */
struct Paso_Solver_GS {
  dim_t n_block;
  dim_t n;
  index_t num_colors;
  index_t* colorOf;
  index_t* main_iptr;
  double* diag;
  Paso_SparseMatrix * factors;
  Paso_Pattern* pattern;
  dim_t sweeps;
  double* x_old;
};
typedef struct Paso_Solver_GS Paso_Solver_GS;

/* RILU preconditioner */
struct Paso_Solver_RILU {
  dim_t n;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  double* inv_A_FF;
  index_t* A_FF_pivot;
  Paso_SparseMatrix * A_FC;
  Paso_SparseMatrix * A_CF;
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  struct Paso_Solver_RILU * RILU_of_Schur;
};
typedef struct Paso_Solver_RILU Paso_Solver_RILU;

/* AMG preconditioner */
struct Paso_Solver_AMG {
  dim_t n;
  dim_t level;
  bool_t coarsest_level;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  double* inv_A_FF;
  index_t* A_FF_pivot;
  Paso_SparseMatrix * A_FC;
  Paso_SparseMatrix * A_CF;
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  Paso_SparseMatrix * A;
  void* solver;
  Paso_Solver_Jacobi* GS;
  struct Paso_Solver_AMG * AMG_of_Schur;
};
typedef struct Paso_Solver_AMG Paso_Solver_AMG;


/* general preconditioner interface */

typedef struct Paso_Solver_Preconditioner {
  dim_t type;
  /* jacobi preconditioner */
  Paso_Solver_Jacobi* jacobi;
  /* ilu preconditioner */
  Paso_Solver_ILU* ilu;
  /* rilu preconditioner */
  Paso_Solver_RILU* rilu;
  /* Gauss-Seidel preconditioner */
  Paso_Solver_GS* gs;
  /* amg preconditioner */
  Paso_Solver_AMG* amg;

} Paso_Solver_Preconditioner;

void Paso_Solver(Paso_SystemMatrix*,double*,double*,Paso_Options*,Paso_Performance* pp);
void Paso_Solver_free(Paso_SystemMatrix*);
err_t Paso_Solver_BiCGStab( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_PCG( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_TFQMR( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_MINRES( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_GMRES(Paso_SystemMatrix * A, double * r, double * x, dim_t *num_iter, double * tolerance,dim_t length_of_recursion,dim_t restart, Paso_Performance* pp);
void Paso_Preconditioner_free(Paso_Solver_Preconditioner*);
void Paso_Solver_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_Solver_solvePreconditioner(Paso_SystemMatrix* A,double*,double*);
void Paso_Solver_applyBlockDiagonalMatrix(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x,double* b);

void Paso_Solver_ILU_free(Paso_Solver_ILU * in);
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveILU(Paso_Solver_ILU * ilu, double * x, double * b);

void Paso_Solver_GS_free(Paso_Solver_GS * in);
Paso_Solver_GS* Paso_Solver_getGS(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveGS(Paso_Solver_GS * gs, double * x, double * b);

void Paso_Solver_RILU_free(Paso_Solver_RILU * in);
Paso_Solver_RILU* Paso_Solver_getRILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveRILU(Paso_Solver_RILU * rilu, double * x, double * b);

void Paso_Solver_AMG_free(Paso_Solver_AMG * in);
Paso_Solver_AMG* Paso_Solver_getAMG(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Solver_solveAMG(Paso_Solver_AMG * amg, double * x, double * b);

void Paso_Solver_updateIncompleteSchurComplement(Paso_SparseMatrix* A_CC, Paso_SparseMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot, Paso_SparseMatrix *A_FC);
Paso_Solver_Jacobi* Paso_Solver_getJacobi(Paso_SparseMatrix * A_p);
void Paso_Solver_solveJacobi(Paso_Solver_Jacobi * prec, double * x, double * b);
void Paso_Solver_Jacobi_free(Paso_Solver_Jacobi * in);

err_t Paso_Solver_GMRES2(Paso_Function * F, const double* f0, const double* x0, double * x, dim_t *iter, double* tolerance, Paso_Performance* pp);
err_t Paso_Solver_NewtonGMRES(Paso_Function *F, double *x, Paso_Options* options, Paso_Performance* pp);

Paso_Function * Paso_Function_LinearSystem_alloc(Paso_SystemMatrix* A, double* b, Paso_Options* options);
err_t Paso_Function_LinearSystem_call(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp);
void Paso_Function_LinearSystem_free(Paso_Function * F);
err_t Paso_Function_LinearSystem_setInitialGuess(Paso_SystemMatrix* A, double* x, Paso_Performance *pp);

#endif /* #ifndef INC_SOLVER */
