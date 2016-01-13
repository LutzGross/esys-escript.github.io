
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix: interface to PaStiX direct solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2010 */
/* Author: Lutz Gross, l.gao@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_PASTIX
#define INC_PASO_PASTIX

#include "SystemMatrix.h"
#include "performance.h"

#ifdef PASTIX
#include "pastix.h"
#include "cscd_utils.h"
#endif

void Paso_PASTIX_free(Paso_SystemMatrix* A);
void Paso_PASTIX(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options,Paso_Performance* pp);

#endif
