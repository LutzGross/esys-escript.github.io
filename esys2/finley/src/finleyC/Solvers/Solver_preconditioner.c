/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix: sets-up the preconditioner           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Solver.h"

/***********************************************************************************/

/*  free space */

void Finley_Preconditioner_free(Finley_Solver_Preconditioner* in) {
    if (in!=NULL) {
      Finley_Solver_ILU_free(in->ilu);
      Finley_Solver_Jacobi_free(in->jacobi);
      MEMFREE(in);
    }
}
/*  call the iterative solver: */

void Finley_Solver_setPreconditioner(Finley_SystemMatrix* A,Finley_SolverOptions* options) {
    Finley_Solver_Preconditioner* prec=NULL;
    if (A->iterative==NULL) {
        /* allocate structure to hold preconditioner */
        prec=MEMALLOC(1,Finley_Solver_Preconditioner);
        if (Finley_checkPtr(prec)) return;
        prec->type=UNKNOWN;
        prec->ilu=NULL;
        prec->jacobi=NULL;
        A->iterative=prec;
        switch (options->preconditioner) {
           default:
           case ESCRIPT_JACOBI:
              if (options->verbose) printf("Jacobi preconditioner is used.\n");
              prec->jacobi=Finley_Solver_getJacobi(A);
              prec->type=ESCRIPT_JACOBI;
              break;
           case ESCRIPT_ILU0:
              if (options->verbose) printf("ILU preconditioner is used.\n");
              prec->ilu=Finley_Solver_getILU(A,options->verbose);
              prec->type=ESCRIPT_ILU0;
              break;
        }
        if (Finley_ErrorCode!=NO_ERROR) {
           Finley_Preconditioner_free(prec);
           A->iterative=NULL;
        }
    }
}

/* applies the preconditioner */
/* has to be called within a parallel reqion */
/* barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Finley_Solver_solvePreconditioner(Finley_SystemMatrix* A,double* x,double* b){
    Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) A->iterative;
    #pragma omp barrier
    switch (prec->type) {
        default:
        case ESCRIPT_JACOBI:
           Finley_Solver_solveJacobi(prec->jacobi,x,b);
           break;
        case ESCRIPT_ILU0:
           Finley_Solver_solveILU(prec->ilu,x,b);
           break;
    }
}
