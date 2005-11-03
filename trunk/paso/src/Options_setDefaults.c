/* $Id$ */

/**************************************************************/

/*   Paso: solver options */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"

/**************************************************************/

/* set the default values for solver options */

void Paso_Options_setDefaults(Paso_Options* options) {

  options->method=PASO_DEFAULT;
  options->package=PASO_DEFAULT;
  options->symmetric=FALSE;
  options->verbose=FALSE;
  options->reordering=PASO_NO_REORDERING;
  options->tolerance=1.E-8;
  options->final_residual=0;
  options->preconditioner=PASO_JACOBI;
  options->iter_max=1000;
  options->iter=0;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->restart=-1;
  options->truncation=20;
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:46  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */

