/* $Id$ */

#ifndef INC_SOLVER
#define INC_SOLVER

#include "System.h"

#define FINLEY_SOLVER_TRACE
/* error codes used in the solver */
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10

/* ILU preconditioner */
struct Finley_Solver_ILU {
  int n;
  int n_block;
  int n_F;
  int n_C;
  double* inv_A_FF;
  maybelong* A_FF_pivot;
  Finley_SystemMatrix * A_FC;
  Finley_SystemMatrix * A_CF;
  maybelong* rows_in_F;
  maybelong* rows_in_C;
  maybelong* mask_F;
  maybelong* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  struct Finley_Solver_ILU * ILU_of_Schur;
};
typedef struct Finley_Solver_ILU Finley_Solver_ILU;


/* jacobi  preconditioner */

typedef struct Finley_Solver_Jacobi {
  int n_block;
  int n;
  double* values;
  maybelong* pivot;
} Finley_Solver_Jacobi;

/* general preconditioner interface */

typedef struct Finley_Solver_Preconditioner {
  int type;
  /* jacobi preconditioner */
  Finley_Solver_Jacobi* jacobi;
  /* ilu preconditioner */
  Finley_Solver_ILU* ilu;
} Finley_Solver_Preconditioner;

void Finley_Solver(Finley_SystemMatrix*,double*,double*,Finley_SolverOptions*);
void Finley_Solver_free(Finley_SystemMatrix*);
int Finley_Solver_BiCGStab( Finley_SystemMatrix * A, double* B, double * X, int *iter, double * tolerance);
int Finley_Solver_PCG( Finley_SystemMatrix * A, double* B, double * X, int *iter, double * tolerance);
int Finley_Solver_GMRES(Finley_SystemMatrix * A, double * r, double * x, int *num_iter, double * tolerance,int length_of_recursion,int restart);
void Finley_Preconditioner_free(Finley_Solver_Preconditioner*);
void Finley_Solver_setPreconditioner(Finley_SystemMatrix* A,Finley_SolverOptions* options);
void Finley_Solver_solvePreconditioner(Finley_SystemMatrix* A,double*,double*);
void Finley_Solver_applyBlockDiagonalMatrix(int n_block,int n,double* D,maybelong* pivot,double* x,double* b);
void Finley_Solver_ILU_free(Finley_Solver_ILU * in);
Finley_Solver_ILU* Finley_Solver_getILU(Finley_SystemMatrix * A_p,int verbose);
void Finley_Solver_solveILU(Finley_Solver_ILU * ilu, double * x, double * b);
void Finley_Solver_updateIncompleteSchurComplement(Finley_SystemMatrix* A_CC,Finley_SystemMatrix *A_CF,double* invA_FF,maybelong* A_FF_pivot,Finley_SystemMatrix *A_FC);
Finley_Solver_Jacobi* Finley_Solver_getJacobi(Finley_SystemMatrix * A_p);
void Finley_Solver_solveJacobi(Finley_Solver_Jacobi * prec, double * x, double * b);
void Finley_Solver_Jacobi_free(Finley_Solver_Jacobi * in);



#endif /* #ifndef INC_SOLVER */
