/* $Id: SCSL_iterative.c 150 2005-09-15 03:44:45Z jgs $ */

/**************************************************************/

/* Paso: SystemMatrix: interface to SGI SCSL iterative solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include <stdlib.h>
#include "Paso.h"
#include "SystemMatrix.h"
#include "SCSL.h"
#ifdef SCSL
#include <scsl_sparse.h>
#endif

/***********************************************************************************/

/*  free any extra stuff possibly used by the SCSL library */

void Paso_SCSL_iterative_free(Paso_SystemMatrix* A) {

}

/*  call the iterative solver: */

void Paso_SCSL_iterative(Paso_SystemMatrix* A,
                           double* out,double* in,Paso_Options* options) {
#ifdef SCSL
    char text2[3];
    int iters,  method,precond,storage,maxiters;
    double drop_tolerance,drop_storage,finalres,convtol,time0;
    if (A->col_block_size!=1) {
      Paso_setError(TYPE_ERROR,"__FILE__: block size 1 is required.");
    }
    switch(A->type) {
    case CSR:
      storage=0;
      break;
    case CSC:
      storage=1;
      break;
    default:
      Paso_setError(TYPE_ERROR,"__FILE__:Matrix type is not supported.");
      return;
    } /* switch A->type */

    method=Paso_Options_getSolver(options->method,PASO_PASO,options->symmetric);
    if (Paso_noError()) {
       switch (method) {
         case PASO_PCG:
	   method=0;
	   break;
         case PASO_CR:
	   method=1;
	   break;
         case PASO_CGS:
	   method=10;
	   break;
         case PASO_BICGSTAB:
	   method=11;
	   break;
         default:
	   method=11;
	   break;
       }
       switch (options->preconditioner) {
         case PASO_JACOBI:
           precond=0;
           break;
         case PASO_SSOR:
           precond=1;
           break;
         case PASO_ILU0:
           if (options->symmetric) precond = 2;
           else precond=0;
           break;
         case PASO_ILUT:
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
       if (options->reordering==PASO_NO_REORDERING) {
         sprintf(text2,"%d",0);
       } else {
         sprintf(text2,"%d",-1);
       }
       setenv("ITERATIVE_RCM",text2,1);
       setenv("ITERATIVE_COPY","1",1);

       time0=Paso_timer();
       DIterative(A->num_rows,A->pattern->ptr,A->pattern->index,A->val,storage,out,in,method,precond,maxiters,convtol,&iters,&finalres);
       options->iter=iters;
       options->final_residual=finalres;
       time0=Paso_timer()-time0;
       if (options->verbose) {
             printf("timing SCSL: solve: %.4e sec\n",time0);
             if (iters>0) printf("timing: per iteration: %.4e sec\n",time0/iters);
       }
    }
#else
    Paso_setError(SYSTEM_ERROR,"__FILE__: SCSL not available.");
#endif
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:40  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:49  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
