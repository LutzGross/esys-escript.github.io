
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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


/****************************************************************************

 * Paso: Merged solver for AMG

 ****************************************************************************/

/* Author: lgao@uq.edu.au, l.gross@uq.edu.au                 */

/****************************************************************************/

#ifndef __PASO_MERGEDSOLVER_H__
#define __PASO_MERGEDSOLVER_H__

#include "Paso.h"
#include "SystemMatrix.h"

namespace paso {

struct MergedSolver
{
    MergedSolver(const_SystemMatrix_ptr A, const Options* options);
    ~MergedSolver();

    void solve(double* local_x, const double* local_b);

    escript::JMPI mpi_info;
    SparseMatrix_ptr A;
    double* x;
    double* b;
    int* counts;
    int* offset;
    index_t reordering;
    index_t refinements;
    index_t verbose;
    index_t sweeps;
};

} // namespace paso

#endif // __PASO_MERGEDSOLVER_H__

