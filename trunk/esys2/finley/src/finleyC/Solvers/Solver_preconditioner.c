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
        prec=(Finley_Solver_Preconditioner*) MEMALLOC(sizeof(Finley_Solver_Preconditioner));
        if (Finley_checkPtr(prec)) return;
        prec->type=UNKNOWN;
        prec->values=NULL;
        prec->numColors=0;
        prec->mainDiag=NULL;
        prec->color=NULL;
        A->iterative=prec;
        switch (options->preconditioner) {
           default:
              printf("Information: Unsupported preconditioner selected. ILU0 is used.\n");
           case ILU0:
              Finley_Solver_setILU0(A);
              prec->type=ILU0;
              break;
           case JACOBI:
              Finley_Solver_setJacobi(A);
              prec->type=JACOBI;
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
        case ILU0:
           Finley_Solver_solveILU0(A,x,b);
           break;
        case JACOBI:
           Finley_Solver_solveJacobi(A,x,b);
           break;
        default:
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"Unknown preconditioner type.");
    }
#endif

}

/*
* $Log$
* Revision 1.1  2004/10/26 06:53:58  jgs
* Initial revision
*
* Revision 1.1  2004/07/02 04:21:14  gross
* Finley C code has been included
*
*
*/
