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

static double ONE=1;
static double ZERO=0;
static double TOLERANCE_FOR_SCALARS=0;

/* ILU preconditioner */
struct Finley_Solver_ILU {
  dim_t n;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  double* inv_A_FF;
  index_t* A_FF_pivot;
  Finley_SystemMatrix * A_FC;
  Finley_SystemMatrix * A_CF;
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  struct Finley_Solver_ILU * ILU_of_Schur;
};
typedef struct Finley_Solver_ILU Finley_Solver_ILU;


/* jacobi  preconditioner */

typedef struct Finley_Solver_Jacobi {
  dim_t n_block;
  dim_t n;
  double* values;
  index_t* pivot;
} Finley_Solver_Jacobi;

/* general preconditioner interface */

typedef struct Finley_Solver_Preconditioner {
  dim_t type;
  /* jacobi preconditioner */
  Finley_Solver_Jacobi* jacobi;
  /* ilu preconditioner */
  Finley_Solver_ILU* ilu;
} Finley_Solver_Preconditioner;

void Finley_Solver(Finley_SystemMatrix*,double*,double*,Finley_SolverOptions*);
void Finley_Solver_free(Finley_SystemMatrix*);
err_t Finley_Solver_BiCGStab( Finley_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance);
err_t Finley_Solver_PCG( Finley_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance);
err_t Finley_Solver_GMRES(Finley_SystemMatrix * A, double * r, double * x, dim_t *num_iter, double * tolerance,dim_t length_of_recursion,dim_t restart);
void Finley_Preconditioner_free(Finley_Solver_Preconditioner*);
void Finley_Solver_setPreconditioner(Finley_SystemMatrix* A,Finley_SolverOptions* options);
void Finley_Solver_solvePreconditioner(Finley_SystemMatrix* A,double*,double*);
void Finley_Solver_applyBlockDiagonalMatrix(dim_t n_block,dim_t n,double* D,index_t* pivot,double* x,double* b);
void Finley_Solver_ILU_free(Finley_Solver_ILU * in);
Finley_Solver_ILU* Finley_Solver_getILU(Finley_SystemMatrix * A_p,bool_t verbose);
void Finley_Solver_solveILU(Finley_Solver_ILU * ilu, double * x, double * b);
void Finley_Solver_updateIncompleteSchurComplement(Finley_SystemMatrix* A_CC,Finley_SystemMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot,Finley_SystemMatrix *A_FC);
Finley_Solver_Jacobi* Finley_Solver_getJacobi(Finley_SystemMatrix * A_p);
void Finley_Solver_solveJacobi(Finley_Solver_Jacobi * prec, double * x, double * b);
void Finley_Solver_Jacobi_free(Finley_Solver_Jacobi * in);



#endif /* #ifndef INC_SOLVER */

/*
 * $Log$
 * Revision 1.7  2005/07/08 04:08:00  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.5  2005/06/29 02:34:59  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.4  2005/03/02 23:35:07  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 * Revision 1.1.1.1.2.3  2005/02/18 03:35:16  gross
 * another function added in prepartion for the reimplementation of ILU
 *
 * Revision 1.1.1.1.2.2  2004/12/07 10:12:06  gross
 * GMRES added
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:21  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:58  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
