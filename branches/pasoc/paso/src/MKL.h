
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


void MKL_free(SparseMatrix* A);
void MKL_solve(SparseMatrix_ptr A, double* out, double* in, index_t reordering,
               dim_t numRefinements, bool verbose);

} // namespace paso

#endif // __PASO_MKL_H__

