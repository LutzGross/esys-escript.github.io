
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to UMFPACK sparse solver */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_UMFPACK_H__
#define __PASO_UMFPACK_H__

#include "SparseMatrix.h"

#ifdef ESYS_HAVE_UMFPACK
#include <umfpack.h>
#endif

namespace paso {

struct UMFPACK_Handler {
    void *symbolic;
    void *numeric;
};

void UMFPACK_free(SparseMatrix* A);
void UMFPACK_solve(SparseMatrix_ptr A, double* out, double* in,
                   dim_t numRefinements, bool verbose);

} // namespace paso

#endif // __PASO_UMFPACK_H__

