/* $Id$ */

#ifndef INC_SOLVER
#define INC_SOLVER

#include "SystemMatrix.h"

#define FINLEY_SOLVER_TRACE
/* error codes used in the solver */
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10

static double ONE=1.;
static double ZERO=0.;
static double TOLERANCE_FOR_SCALARS=0.;

/* ILU preconditioner */
struct Paso_Solver_ILU {
  dim_t n;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  double* inv_A_FF;
  index_t* A_FF_pivot;
  Paso_SystemMatrix * A_FC;
  Paso_SystemMatrix * A_CF;
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  struct Paso_Solver_ILU * ILU_of_Schur;
};
typedef struct Paso_Solver_ILU Paso_Solver_ILU;


/* jacobi  preconditioner */

typedef struct Paso_Solver_Jacobi {
  dim_t n_block;
  dim_t n;
  double* values;
  index_t* pivot;
} Paso_Solver_Jacobi;

/* general preconditioner interface */

typedef struct Paso_Solver_Preconditioner {
  dim_t type;
  /* jacobi preconditioner */
  Paso_Solver_Jacobi* jacobi;
  /* ilu preconditioner */
  Paso_Solver_ILU* ilu;
} Paso_Solver_Preconditioner;

void Paso_Solver(Paso_SystemMatrix*,double*,double*,Paso_Options*);
void Paso_Solver_free(Paso_SystemMatrix*);
err_t Paso_Solver_BiCGStab( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance);
err_t Paso_Solver_PCG( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance);
err_t Paso_Solver_GMRES(Paso_SystemMatrix * A, double * r, double * x, dim_t *num_iter, double * tolerance,dim_t length_of_recursion,dim_t restart);
void Paso_Preconditioner_free(Paso_Solver_Preconditioner*);
void Paso_Solver_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_Solver_solvePreconditioner(Paso_SystemMatrix* A,double*,double*);
void Paso_Solver_applyBlockDiagonalMatrix(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x,double* b);
void Paso_Solver_ILU_free(Paso_Solver_ILU * in);
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SystemMatrix * A_p,bool_t verbose);
void Paso_Solver_solveILU(Paso_Solver_ILU * ilu, double * x, double * b);
void Paso_Solver_updateIncompleteSchurComplement(Paso_SystemMatrix* A_CC,Paso_SystemMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot,Paso_SystemMatrix *A_FC);
Paso_Solver_Jacobi* Paso_Solver_getJacobi(Paso_SystemMatrix * A_p);
void Paso_Solver_solveJacobi(Paso_Solver_Jacobi * prec, double * x, double * b);
void Paso_Solver_Jacobi_free(Paso_Solver_Jacobi * in);

#endif /* #ifndef INC_SOLVER */
