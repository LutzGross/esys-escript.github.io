/* $Id$ */

/**************************************************************/

/* Finley: interface to the direct solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "System.h"

#if DIRECT_SOLVER == SGI_SCSL
#include "SCSL/SCSL.h"
#else 
#include "Solvers/Solver.h"
#endif

#if ITERATIVE_SOLVER == SGI_SCSL
#include "SCSL/SCSL.h"
#else
#include "Solvers/Solver.h"
#endif

/**************************************************************/

void Finley_SystemMatrix_solve(Finley_SystemMatrix* A,
                               escriptDataC* out,
                               escriptDataC* in,
                               Finley_SolverOptions* options) {
  if (A->num_rows!=A->num_cols || A->col_block_size!=A->row_block_size) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"matrix has to be a square matrix.");
  } else if (A->num_rows*A->row_block_size!=getLength(in) ) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"number of rows and length of right hand side don't match.");
  } else if (A->num_cols*A->col_block_size!=getLength(out)  ) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"number of rows and length of solution vector don't match.");
  } else if (!isExpanded(in)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"input Data object has to be expanded");
  } else if (!isExpanded(out)) {
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"output Data object has to be expanded");
  }
  if (options->method==ESCRIPT_DIRECT || options->method==ESCRIPT_CHOLEVSKY) {
      #ifdef Finley_TRACE
         printf("Finley calls direct solver.\n");
      #endif
      #if DIRECT_SOLVER == SGI_SCSL
         Finley_SCSL(A,getSampleData(out,0),getSampleData(in,0),options);
      #else 
         Finley_Solver(A,getSampleData(out,0),getSampleData(in,0),options);
      #endif
  }  else {
     #ifdef Finley_TRACE
     printf("Iterative solver is called with tolerance %e.\n",options->tolerance);
     #endif
     #if ITERATIVE_SOLVER == SGI_SCSL
        Finley_SCSL(A,getSampleData(out,0),getSampleData(in,0),options);
     #else
        Finley_Solver(A,getSampleData(out,0),getSampleData(in,0),options);
     #endif
  }
}

/*  free memory possibly resereved for a recall */

void Finley_SystemMatrix_solve_free(Finley_SystemMatrix* in) { 
            #if DIRECT_SOLVER == SGI_SCSL 
                Finley_SCSL_free(in);
            #else
                Finley_Solver_free(in);
            #endif
            #if ITERATIVE_SOLVER == SGI_SCSL 
                Finley_SCSL_free(in);
            #else
                Finley_Solver_free(in);
            #endif
}
/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:31  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:19  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
