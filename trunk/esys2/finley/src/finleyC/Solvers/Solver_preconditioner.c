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
* Revision 1.2  2004/12/14 05:39:32  jgs
* *** empty log message ***
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:17  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
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
