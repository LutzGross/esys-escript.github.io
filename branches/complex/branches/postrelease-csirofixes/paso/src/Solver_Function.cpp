
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

#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"

namespace paso {

/*
 * generate Linear System (mainly for test purposes)
 */
LinearSystem::LinearSystem(SystemMatrix_ptr A, double* _b, Options* options) :
    Function(A->mpi_info)
{
    A->setPreconditioner(options);
    n = A->getTotalNumRows();
    mat = A;
    b = _b;
    tmp = new double[n];
}

LinearSystem::~LinearSystem()
{
    delete[] tmp;
}

/*
 * evaluates value=P*(b-Ax)
 */
err_t LinearSystem::call(double* value, const double* arg, Performance* pp)
{
    // tmp = b
    util::copy(n, tmp, b);
    // tmp = (A*arg-tmp)
    SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, mat, arg, -PASO_ONE, tmp);
    // value = P*tmp
    mat->solvePreconditioner(value, tmp);
    return NO_ERROR;
}

} // namespace paso

