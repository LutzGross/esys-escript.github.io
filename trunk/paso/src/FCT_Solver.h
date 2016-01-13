
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

#ifndef __PASO_FCTSOLVER_H__
#define __PASO_FCTSOLVER_H__

#include "Transport.h"
#include "FluxLimiter.h"
#include "Solver.h"

namespace paso {

PASO_DLL_API
struct FCT_Solver
{
    FCT_Solver(const_TransportProblem_ptr tp, Options* options);

    ~FCT_Solver();

    err_t update(double* u, double* u_old, Options* options, Performance* pp);

    err_t updateNL(double* u, double* u_old, Options* options, Performance* pp);

    err_t updateLCN(double* u, double* u_old, Options* options, Performance* pp);

    void initialize(double dt, Options* options, Performance* pp);

    static double getSafeTimeStepSize(TransportProblem_ptr tp);

    static void setLowOrderOperator(TransportProblem_ptr tp);

    void setAntiDiffusionFlux_linearCN(SystemMatrix_ptr flux_matrix);

    void setAntiDiffusionFlux_BE(SystemMatrix_ptr flux_matrix);

    void setAntiDiffusionFlux_CN(SystemMatrix_ptr flux_matrix);

    void setMuPaLu(double* out, const_Coupler_ptr coupler, double a);

    inline double getTheta()
    {
        return method == PASO_BACKWARD_EULER ? 1. : 0.5;
    }

    const_TransportProblem_ptr transportproblem;
    esysUtils::JMPI mpi_info;
    FCT_FluxLimiter* flux_limiter;
    index_t method;
    double omega;
    double dt;
    double* b;
    double* z;
    double* du;
    Coupler_ptr u_coupler;
    Coupler_ptr u_old_coupler; /* last time step */
};


} // namespace paso

#endif // __PASO_FCTSOLVER_H__

