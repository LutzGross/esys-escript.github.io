
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/* Paso: SystemMatrix: interface to intel UMFPACK sparse solver */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

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

void Paso_UMFPACK_free(paso::SparseMatrix* A);
void Paso_UMFPACK(paso::SparseMatrix* A, double* out, double* in, dim_t numRefinements, bool verbose);
#endif
