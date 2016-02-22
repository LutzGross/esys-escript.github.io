
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


#ifndef __PASO_REACTIVESOLVER_H__
#define __PASO_REACTIVESOLVER_H__

#include "performance.h"
#include "Transport.h"

namespace paso {

PASO_DLL_API
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

