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
 * Revision 1.2  2004/12/14 05:39:31  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.2  2004/12/07 10:12:05  gross
 * GMRES added
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:19  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/08/28 12:58:08  gross
 * SimpleSolve is not running yet: problem with == of functionsspace
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */

