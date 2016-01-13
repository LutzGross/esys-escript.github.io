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

/* Paso: returns the package to be used                   */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"

/**************************************************************/

index_t Paso_Options_getPackage(index_t solver,index_t package, bool_t symmetry) {
  index_t out=PASO_PASO;
  if (package==PASO_DEFAULT) {
      if (solver==PASO_DIRECT) {
         #ifdef MKL
            out=PASO_MKL;
         #else
            #ifdef SCSL
              out=PASO_SCSL;
            #else
              #ifdef UMFPACK
                out=PASO_UMFPACK;
              #endif
            #endif
         #endif
      } else {
         out=PASO_PASO;
      }
  } else if (package==PASO_PASO) {
      out=PASO_PASO;
  } else if (package==PASO_SCSL) {
      out=PASO_SCSL;
  } else if (package==PASO_MKL) {
      out=PASO_MKL;
  } else if (package==PASO_UMFPACK) {
      out=PASO_UMFPACK;
  } else if (package==PASO_TRILINOS) {
      out=PASO_TRILINOS;
  } else {
      Paso_setError(VALUE_ERROR,"Unidentified package.");
  }
  return out;
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:46  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
