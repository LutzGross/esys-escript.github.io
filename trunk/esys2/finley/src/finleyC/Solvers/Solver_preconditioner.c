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

void Finley_Preconditioner_free(void* in) {
#if ITERATIVE_SOLVER == NO_LIB
    if (in!=NULL) {
       Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) in;
       MEMFREE(prec->values);
       MEMFREE(prec->mainDiag);
       MEMFREE(prec->pivot);
       MEMFREE(prec->color);
       MEMFREE(prec);
    }
#endif
}
/*  call the iterative solver: */

void Finley_Solver_setPreconditioner(Finley_SystemMatrix* A,Finley_SolverOptions* options) {
#if ITERATIVE_SOLVER == NO_LIB
    Finley_Solver_Preconditioner* prec=NULL;
    if (A->iterative==NULL) {
        /* allocate structure to hold preconditioner */
        prec=MEMALLOC(1,Finley_Solver_Preconditioner);
        if (Finley_checkPtr(prec)) return;
        prec->type=UNKNOWN;
        prec->numColors=0;

        prec->values=NULL;
        prec->mainDiag=NULL;
        prec->pivot=NULL;
        prec->color=NULL;

        A->iterative=prec;
        switch (options->preconditioner) {
           default:
           case ESCRIPT_JACOBI:
              if (options->verbose) printf("Jacobi preconditioner is used.\n");
              Finley_Solver_setJacobi(A);
              prec->type=ESCRIPT_JACOBI;
              break;
           case ESCRIPT_ILU0:
              if (options->verbose) printf("ILU(0) preconditioner is used.\n");
              Finley_Solver_setILU0(A);
              prec->type=ESCRIPT_ILU0;
              break;
        }
        if (Finley_ErrorCode!=NO_ERROR) {
           Finley_Preconditioner_free(prec);
           A->iterative=NULL;
        }
    }
#endif
}

/* applies the preconditioner */
/* has to be called within a parallel reqion */
/* barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Finley_Solver_solvePreconditioner(Finley_SystemMatrix* A,double* x,double* b){
#if ITERATIVE_SOLVER == NO_LIB
    Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) A->iterative;
    #pragma omp barrier
    switch (prec->type) {
        default:
        case ESCRIPT_JACOBI:
           Finley_Solver_solveJacobi(A,x,b);
           break;
        case ESCRIPT_ILU0:
           Finley_Solver_solveILU0(A,x,b);
           break;
    }
#endif

}

/*
* $Log$
* Revision 1.4  2004/12/15 07:08:35  jgs
* *** empty log message ***
*
*
*
*/
