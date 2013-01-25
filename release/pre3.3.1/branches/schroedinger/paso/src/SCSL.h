
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

/* Paso: SystemMatrix: interface to SGI SCSL sparse solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SCSL
#define INC_PASO_SCSL

#include "SystemMatrix.h"
#include "performance.h"

void Paso_SCSL_free(Paso_SystemMatrix* A);
void Paso_SCSL(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options,Paso_Performance* pp);
void Paso_SCSL_iterative_free(Paso_SystemMatrix* A);
void Paso_SCSL_iterative(Paso_SystemMatrix* A, double* out,double* in,Paso_Options* options,Paso_Performance* pp);
void Paso_SCSL_direct_free(Paso_SystemMatrix* A);
void Paso_SCSL_direct(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options,Paso_Performance* pp);

#endif
