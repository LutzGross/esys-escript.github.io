/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: interface to the direct solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "performance.h"
#include "Solver.h"

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
  Paso_Performance pp;
  Paso_resetError();
  if (A->num_rows!=A->num_cols || A->col_block_size!=A->row_block_size) {
       Paso_setError(VALUE_ERROR,"Paso_solve: matrix has to be a square matrix.");
       return;
  }
  Performance_open(&pp,options->verbose);
  index_t package=Paso_Options_getPackage(options->method,options->package,options->symmetric);
  if (Paso_noError()) {
     switch(package) {

        case PASO_PASO:
          Paso_Solver(A,out,in,options,&pp);
          A->solver_package=PASO_PASO;
          break;

        #ifdef SCSL
        case PASO_SCSL:
          Paso_SCSL(A,out,in,options,&pp);
          A->solver_package=PASO_SCSL;
          break;
        #endif

  
        #ifdef MKL
        case PASO_MKL:
          Paso_MKL(A,out,in,options,&pp);
          A->solver_package=PASO_MKL;
          break;
        #endif

/*
        case PASO_UMFPACK:
          Paso_UMFPACK(A,out,in,options);
          A->solver_package=PASO_UMFPACK;
          break;
*/

        default:
           Paso_setError(VALUE_ERROR,"Paso_solve: unknown package code");
           break;
     }
  }
  Performance_close(&pp,options->verbose);
  return;
}

/*  free memory possibly resereved for a recall */

void Paso_solve_free(Paso_SystemMatrix* in) { 

     switch(in->solver_package) {

        case PASO_PASO:
          Paso_Solver_free(in);
          break;

        #ifdef SCSL
        case PASO_SCSL:
          Paso_SCSL_free(in);
          break;
        #endif

  
        #ifdef MKL
        case PASO_MKL:
          Paso_MKL_free(in); 
          break;
        #endif

/*
        case PASO_UMFPACK:
          Paso_UMFPACK_free(in); 
          break;
*/
   }
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
