/* $Id$ */

/**************************************************************/

/* Finley: interface to the SCSL solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "System.h"
#include "SCSL/SCSL.h"

/**************************************************************/

void Finley_SCSL(Finley_SystemMatrix* A,
                               double* out,
                               double* in,
                               Finley_SolverOptions* options) {

  if (options->method==ESCRIPT_DIRECT || options->method==ESCRIPT_CHOLEVSKY) {
      #if DIRECT_SOLVER == SGI_SCSL
         Finley_SCSL_direct(A,out,in,options);
      #else
         Finley_ErrorCode=SYSTEM_ERROR;
         sprintf(Finley_ErrorMsg,"No direct solver available!");
      #endif
  }  else {
     #if ITERATIVE_SOLVER == SGI_SCSL
         Finley_SCSL_iterative(A,out,in,options);
     #else
         Finley_ErrorCode=SYSTEM_ERROR;
         sprintf(Finley_ErrorMsg,"No iterative solver available!");
     #endif
  }
}

/*  free memory possibly resereved for a recall */

void Finley_SCSL_free(Finley_SystemMatrix* in) {
    #if DIRECT_SOLVER == SGI_SCSL
       Finley_SCSL_direct_free(in);
    #endif
    #if ITERATIVE_SOLVER == SGI_SCSL
        Finley_SCSL_iterative_free(in);
    #endif
}
/*
 * $Log$
 * Revision 1.2  2004/12/15 07:08:34  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.1  2004/11/12 06:58:20  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 *
 */
