
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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

void PASO_DLL_API UMFPACK_free(SparseMatrix<double>* A);
void PASO_DLL_API UMFPACK_free(SparseMatrix<cplx_t>* A);

void PASO_DLL_API UMFPACK_solve(SparseMatrix_ptr<double> A, double* out, double* in,
                   dim_t numRefinements, bool verbose);
void PASO_DLL_API UMFPACK_solve(SparseMatrix_ptr<cplx_t> A, cplx_t* out, cplx_t* in,
                   dim_t numRefinements, bool verbose);

} // namespace paso

#endif // __PASO_UMFPACK_H__

