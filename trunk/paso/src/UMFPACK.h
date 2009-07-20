
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix: interface to intel UMFPACK sparse solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_UMFPACK
#define INC_PASO_UMFPACK

#include "SystemMatrix.h"
#include "performance.h"

#ifdef UMFPACK
#include <umfpack.h>
#endif

typedef struct {
    void *symbolic;
    void *numeric;
} Paso_UMFPACK_Handler;

void Paso_UMFPACK_free(Paso_SystemMatrix* A);
void Paso_UMFPACK(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options,Paso_Performance* pp);
void Paso_UMFPACK1(Paso_SparseMatrix* A, double* out, double* in, bool_t verbose);
#endif
