/* $Id$ */

/*
********************************************************************************
*               Copyright © 2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

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
  options->iter_max=10000;
  options->iter=0;
  options->drop_tolerance=0.01;
  options->drop_storage=2.;
  options->restart=-1;
  options->truncation=20;
}
