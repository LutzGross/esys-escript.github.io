
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/* Interface to Sandia TRILINOS sparse solver */

/* Author: k.steube@uq.edu.au */


#ifndef INC_PASO_TRILINOS
#define INC_PASO_TRILINOS

#include "Paso.h"
#include "performance.h"
#include "escript/system_dep.h"
#include "SystemMatrixPattern.h"
#include "SystemMatrix.h"
#include "Options.h"


void Paso_TRILINOS_alloc(void* trilinos_data, Paso_SystemMatrixPattern *pattern, dim_t row_block_size, dim_t col_block_size);

void Paso_TRILINOS(Paso_SystemMatrix* A,
                   double* out,
                   double* in,
                   Paso_Options* options,
                   Paso_Performance* pp);

void Paso_TRILINOS_free(void* in);

#endif /* ifndef INC_PASO_TRILINOS */
