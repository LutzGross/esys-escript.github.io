/* $Id$ */

/**************************************************************/

/*   Finley: solver options */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */

/**************************************************************/

#include "Common.h"
#include "System.h"

/**************************************************************/

/* set the default values for solver options */

void Finley_SystemMatrix_setDefaults(Finley_SolverOptions* options) {
  options->method=ESCRIPT_DEFAULT_METHOD;
  options->symmetric=FALSE;
  options->verbose=TRUE;
  options->reordering=ESCRIPT_NO_REORDERING;
  options->tolerance=1.E-8;
  options->final_residual=0;
  options->preconditioner=ESCRIPT_JACOBI;
  options->iter_max=1000;
  options->iter=0;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->restart=-1;
  options->truncation=20;
}

/*
 * $Log$
 * Revision 1.4  2004/12/15 07:08:34  jgs
 * *** empty log message ***
 *
 *
 *
 */

