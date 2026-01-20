
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


#ifndef __PASO_FUNCTIONS_H__
#define __PASO_FUNCTIONS_H__

#include "Paso.h"
#include "performance.h"
#include "SystemMatrix.h"

namespace paso {

struct Function
{
    Function(const escript::JMPI& mpi_info);
    virtual ~Function();

    /// sets value=F(arg)
    virtual SolverResult call(double* value, const double* arg, Performance* pp) = 0;

    /// numerical calculation of the directional derivative J0w of F at x0 in
    /// the direction w. f0 is the value of F at x0. setoff is workspace
    SolverResult derivative(double* J0w, const double* w, const double* f0,
                            const double* x0, double* setoff, Performance* pp);

    /// returns the length of the vectors used by this function
    virtual dim_t getLen() = 0;

    const escript::JMPI mpi_info;
};

struct LinearSystem : public Function
{
    LinearSystem(SystemMatrix_ptr<double> A, double* b, Options* options);
    virtual ~LinearSystem();

    virtual SolverResult call(double* value, const double* arg, Performance* pp);

    virtual dim_t getLen() { return n; }

    SystemMatrix_ptr<double> mat;
    double* tmp;
    double* b;
    dim_t n;
};

} // namespace paso

#endif // __PASO_FUNCTIONS_H__

