
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


#ifndef __PASO_REACTIVESOLVER_H__
#define __PASO_REACTIVESOLVER_H__

#include "performance.h"
#include "Transport.h"

// exp(h)-1 ~ h + h**2/2 for abs(h) <  PASO_RT_EXP_LIM_MIN
#define PASO_RT_EXP_LIM_MIN sqrt(EPSILON)

// it is assumed that exp(h) with  h>PASO_RT_EXP_LIM_MAX is not reliable
#define PASO_RT_EXP_LIM_MAX log(1./sqrt(EPSILON))

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

    err_t solve(double* u, double* u_old, const double* source,
                Options* options, Performance* pp);

    static double getSafeTimeStepSize(const_TransportProblem_ptr tp);

    const_TransportProblem_ptr tp;
    double dt;
};



} // namespace paso

#endif // __PASO_REACTIVESOLVER_H__

