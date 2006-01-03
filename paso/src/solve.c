/* $Id$ */

/**************************************************************/

/* Paso: interface to the direct solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Solvers/Solver.h"

#ifdef SCSL
#include "SCSL.h"
#endif

#ifdef MKL
#include "MKL.h"
#endif

#ifdef UMFPACK
#include "UMFPACK.h"
#endif

/**************************************************************/

void Paso_solve(Paso_SystemMatrix* A,
                               double* out,
                               double* in,
                               Paso_Options* options) {

  Paso_resetError();
  if (A->num_rows!=A->num_cols || A->col_block_size!=A->row_block_size) {
       Paso_setError(VALUE_ERROR,"__FILE__: matrix has to be a square matrix.");
       return;
  }
  index_t package=Paso_Options_getPackage(options->method,options->package,options->symmetric);
  if (Paso_noError()) {
     switch(package) {

        case PASO_PASO:
          Paso_Solver(A,out,in,options);
          break;

        #ifdef SCSL
        case PASO_SCSL:
          Paso_SCSL(A,out,in,options);
          break;
        #endif

/*
        case PASO_MKL:
          Paso_MKL(A,out,in,options);
          break;
*/

/*
        case PASO_UMFPACK:
          Paso_UMFPACK(A,out,in,options);
          break;
*/

        default:
           Paso_setError(VALUE_ERROR,"__FILE__: unknown package code");
           return;
     }
  }
  return;
}

/*  free memory possibly resereved for a recall */

void Paso_solve_free(Paso_SystemMatrix* in) { 
          Paso_Solver_free(in);
          #ifdef SCSL
          Paso_SCSL_free(in);
          #endif
          /* Paso_MKL_free(A,out,in,options); */
          /* Paso_UMFPACK_free(A,out,in,options); */
          return;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
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
