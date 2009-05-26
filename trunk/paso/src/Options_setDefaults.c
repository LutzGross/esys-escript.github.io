
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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
  options->absolute_tolerance=0.;
  options->inner_tolerance=0.9;
  options->adapt_inner_tolerance=TRUE;
  options->final_residual=0;
  options->preconditioner=PASO_JACOBI;
  options->iter_max=10000;
  options->inner_iter_max=10;
  options->iter=0;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->restart=-1;
  options->truncation=20;
  options->sweeps=2;
  options->couplingParam=0.05;
  options->AMGlevels=2;
}
