/* $Id$ */

#ifndef INC_SOLVER
#define INC_SOLVER

#define FINLEY_SOLVER_TRACE
/* error codes used in the solver */
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10


typedef struct Finley_Solver_Preconditioner {
  int type;
  double* values;
  int numColors;
  maybelong *mainDiag;
  maybelong *color;
  int *pivot;
} Finley_Solver_Preconditioner;



void Finley_Solver(Finley_SystemMatrix*,double*,double*,Finley_SolverOptions*);
void Finley_Solver_free(Finley_SystemMatrix*);
int Finley_Solver_BiCGStab( Finley_SystemMatrix * A, double* B, double * X, int *iter, double * tolerance);
int Finley_Solver_PCG( Finley_SystemMatrix * A, double* B, double * X, int *iter, double * tolerance);
int Finley_Solver_GMRES(Finley_SystemMatrix * A, double * r, double * x, int *num_iter, double * tolerance,int length_of_recursion,int restart);
void Finley_Preconditioner_free(void* in);
void Finley_Solver_setPreconditioner(Finley_SystemMatrix* A,Finley_SolverOptions* options);
void Finley_Solver_solvePreconditioner(Finley_SystemMatrix* A,double*,double*);

void Finley_Solver_setJacobi(Finley_SystemMatrix*);
void Finley_Solver_solveJacobi(Finley_SystemMatrix*, double*, double*);
void Finley_Solver_setILU0(Finley_SystemMatrix*);
void Finley_Solver_solveILU0(Finley_SystemMatrix*, double*, double*);

void Finley_Solver_getMainDiagonal(Finley_SystemMatrix*,maybelong*);
void Finley_Solver_mis(Finley_SystemMatrixPattern* pattern_p, maybelong* mis_marker);

void Finley_Solver_coloring(Finley_SystemMatrixPattern*,maybelong*,maybelong*);

#endif /* #ifndef INC_SOLVER */
