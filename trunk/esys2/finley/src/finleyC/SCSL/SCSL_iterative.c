/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix: interface to SGI SCSL iterative solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include <stdlib.h>
#include "Finley.h"
#include "System.h"
#include "SCSL.h"
#if ITERATIVE_SOLVER == SGI_SCSL
#include <scsl_sparse.h>
#endif

/***********************************************************************************/

/*  free any extra stuff possibly used by the SCSL library */

void Finley_SCSL_iterative_free(Finley_SystemMatrix* A) {

}

/*  call the iterative solver: */

void Finley_SCSL_iterative(Finley_SystemMatrix* A,
                           double* out,double* in,Finley_SolverOptions* options) {
#if ITERATIVE_SOLVER == SGI_SCSL
    char text2[3];
    int iters,  method,precond,storage,maxiters;
    double drop_tolerance,drop_storage,finalres,convtol,time0;
    if (A->col_block_size!=1) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"SCSL linear solver can only be applied to block size 1.");
      return;
    }
    switch(A->type) {
    case CSR:
      storage=0;
      break;
    case CSC:
      storage=1;
      break;
    default:
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Matrix type %d is not supported by SCSL iterative solver.",A->type);
      return;
    } /* switch A->type */

    if (options->symmetric) {
      switch (options->method) {
      case ESCRIPT_PCG:
	method=0;
	break;
      case ESCRIPT_CR:
	method=1;
	break;
      default:
	method=0;
	break;
      }
    } else {
      switch (options->method) {
      case ESCRIPT_CGS:
	method=10;
	break;
      case ESCRIPT_BICGSTAB:
	method=11;
	break;
      default:
	method=11;
	break;
      }
    }
    switch (options->preconditioner) {
    case ESCRIPT_JACOBI:
      precond=0;
      break;
    case ESCRIPT_SSOR:
      precond=1;
      break;
    case ESCRIPT_ILU0:
      if (options->symmetric) precond = 2;
      else precond=0;
      break;
    case ESCRIPT_ILUT:
      if (options->symmetric) precond = 3;
      else precond=0;
      break;
    default:
      precond=0;
      break;
    }
    maxiters=options->iter_max;
    convtol=options->tolerance;

    drop_tolerance=options->drop_tolerance;
    DIterative_DropTol(drop_tolerance);
    drop_storage=options->drop_storage;
    DIterative_DropStorage(drop_storage);

    if (options->verbose) {
      setenv("ITERATIVE_VERBOSE","1",1);
    } else {
      unsetenv("ITERATIVE_VERBOSE");
    }
    if (options->reordering==ESCRIPT_NO_REORDERING) {
      sprintf(text2,"%d",0);
    } else {
      sprintf(text2,"%d",-1);
    }
    setenv("ITERATIVE_RCM",text2,1);
    setenv("ITERATIVE_COPY","1",1);

    time0=Finley_timer();
    DIterative(A->num_rows,A->pattern->ptr,A->pattern->index,A->val,storage,out,in,method,precond,maxiters,convtol,&iters,&finalres);
    options->iter=iters;
    options->final_residual=finalres;
    time0=Finley_timer()-time0;
    if (options->verbose) {
          printf("timing SCSL: solve: %.4e sec\n",time0);
          if (iters>0) printf("timing: per iteration: %.4e sec\n",time0/iters);
    }
#else
    Finley_ErrorCode=SYSTEM_ERROR;
    sprintf(Finley_ErrorMsg,"No SCSL iterative solver available.");
#endif
}
/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:31  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:20  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
