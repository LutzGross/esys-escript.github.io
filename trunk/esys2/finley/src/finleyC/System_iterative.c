/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix: interface to iterative solver        */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "escript/Data/DataC.h"
#include "System.h"

#if ITERATIVE_SOLVER == SGI_SCSL
#include "SCSL/SCSL.h"
#else
#include "Solvers/Solver.h"
#endif

/***********************************************************************************/

/*  call the iterative solver: */

void Finley_SystemMatrix_iterative(Finley_SystemMatrix* A,escriptDataC* out,escriptDataC* in,Finley_SolverOptions* options) {
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

  if (Finley_ErrorCode==NO_ERROR) {
     #ifdef Finley_TRACE
     printf("Iterative solver is called with tolerance %e.\n",options->tolerance);
     #endif
     #ifdef Finley_TRACE
     printf("direct solver is called.\n");
     #endif
     #if ITERATIVE_SOLVER == SGI_SCSL
     Finley_SCSL_iterative(A,getSampleData(out,0),getSampleData(in,0));
     #else
     Finley_Solver(A,getSampleData(out,0),getSampleData(in,0),options);
     #endif
  }
}

/*  free memory possibly resereved for a recall */

void Finley_SystemMatrix_iterative_free(Finley_SystemMatrix* in) {
  #if ITERATIVE_SOLVER == SGI_SCSL
  Finley_SCSL_iterative_free(in);
  #else
  Finley_Solver_free(in);
  #endif
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
