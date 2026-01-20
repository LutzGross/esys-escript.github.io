
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


#ifndef __PASO_REACTIVESOLVER_H__
#define __PASO_REACTIVESOLVER_H__

#include "Transport.h"

namespace paso {

struct Performance;

struct ReactiveSolver
{
    ReactiveSolver(const_TransportProblem_ptr _tp) : tp(_tp) {}
    ~ReactiveSolver() {}

    inline void initialize(double _dt, Options*)
    {
        dt = _dt;
    }

    SolverResult solve(double* u, double* u_old, const double* source,
                Options* options, Performance* pp);

    static double getSafeTimeStepSize(const_TransportProblem_ptr tp);

    const_TransportProblem_ptr tp;
    double dt;
};


} // namespace paso

#endif // __PASO_REACTIVESOLVER_H__

