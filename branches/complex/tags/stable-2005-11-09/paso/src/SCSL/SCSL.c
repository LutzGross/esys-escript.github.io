/* $Id$ */

/**************************************************************/

/* Paso: interface to the SGI SCSL library */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: gross@@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "SCSL.h"

/**************************************************************/

/*  free any extra stuff possibly used by the SCSL library */

void Paso_SCSL_free(Paso_SystemMatrix* A) {
      Paso_SCSL_direct_free(A);
      Paso_SCSL_iterative_free(A);
}
/*  call the solver: */

void Paso_SCSL(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options) {

  index_t method=Paso_Options_getSolver(options->method,PASO_SCSL,options->symmetric);

  if (Paso_noError()) {
      if (method==PASO_CHOLEVSKY || method==PASO_DIRECT) {
          Paso_SCSL_direct(A,out,in,options);
      } else {
          Paso_SCSL_iterative(A,out,in,options);
      }
  }
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.3  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.2  2005/09/05 10:05:06  gross
 * naming error fixed
 *
 * Revision 1.1.2.1  2005/09/05 06:29:49  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
