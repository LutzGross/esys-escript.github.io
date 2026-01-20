
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to intel MKL sparse solver */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_MKL_H__
#define __PASO_MKL_H__

#include "SparseMatrix.h"

namespace paso {

#ifdef ESYS_INDEXTYPE_LONG
#define ES_PARDISO pardiso_64
#define ES_MKL_INT MKL_INT64
#else
#define ES_PARDISO pardiso
#define ES_MKL_INT MKL_INT
#endif

#ifdef ESYS_HAVE_MKL
#include <mkl_pardiso.h>
#endif


#define MKL_ERROR_NO 0
#define MKL_MTYPE_REAL_SYM -2
#define MKL_MTYPE_REAL_UNSYM 11

#define MKL_REORDERING_MINIMUM_DEGREE 0
#define MKL_REORDERING_NESTED_DISSECTION 2
#define MKL_REORDERING_NESTED_DISSECTION_OMP 3
#define MKL_PHASE_SYMBOLIC_FACTORIZATION 11
#define MKL_PHASE_FACTORIZATION 22
#define MKL_PHASE_SOLVE 33
#define MKL_PHASE_RELEASE_MEMORY -1


void PASO_DLL_API MKL_free(SparseMatrix<double>* A);
void PASO_DLL_API MKL_free(SparseMatrix<cplx_t>* A);

void PASO_DLL_API MKL_solve(SparseMatrix_ptr<double> A, double* out, double* in, index_t reordering,
               dim_t numRefinements, bool verbose);
void PASO_DLL_API MKL_solve(SparseMatrix_ptr<cplx_t> A, cplx_t* out, cplx_t* in, index_t reordering,
               dim_t numRefinements, bool verbose);

} // namespace paso

#endif // __PASO_MKL_H__

