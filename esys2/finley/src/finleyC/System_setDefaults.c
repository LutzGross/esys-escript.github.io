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
  options->verbose=FALSE;
  options->reordering=NO_REORDERING;
  options->tolerance=1.E-8;
  options->final_residual=0;
  options->iterative_method=BICGSTAB;
  options->preconditioner=JACOBI;
  options->iter_max=1000;
  options->iter=0;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->iterative=FALSE;
}

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/08/28 12:58:08  gross
 * SimpleSolve is not running yet: problem with == of functionsspace
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */

