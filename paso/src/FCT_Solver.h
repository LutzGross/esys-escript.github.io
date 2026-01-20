
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

#ifndef __PASO_FCTSOLVER_H__
#define __PASO_FCTSOLVER_H__

#include "Transport.h"
#include "FluxLimiter.h"
#include "Solver.h"

namespace paso {

struct PASO_DLL_API FCT_Solver
{
    FCT_Solver(const_TransportProblem_ptr tp, Options* options);

    ~FCT_Solver();

    SolverResult update(double* u, double* u_old, Options* options, Performance* pp);

    SolverResult updateNL(double* u, double* u_old, Options* options, Performance* pp);

    SolverResult updateLCN(double* u, double* u_old, Options* options, Performance* pp);

    void initialize(double dt, Options* options, Performance* pp);

    static double getSafeTimeStepSize(const_TransportProblem_ptr tp);

    static void setLowOrderOperator(TransportProblem_ptr tp);

    void setAntiDiffusionFlux_linearCN(SystemMatrix_ptr<double> flux_matrix);

    void setAntiDiffusionFlux_BE(SystemMatrix_ptr<double> flux_matrix);

    void setAntiDiffusionFlux_CN(SystemMatrix_ptr<double> flux_matrix);

    void setMuPaLu(double* out, const_Coupler_ptr<real_t> coupler, double a);

    inline double getTheta()
    {
        return method == PASO_BACKWARD_EULER ? 1. : 0.5;
    }

    const_TransportProblem_ptr transportproblem;
    escript::JMPI mpi_info;
    FCT_FluxLimiter* flux_limiter;
    index_t method;
    double omega;
    double dt;
    double* b;
    double* z;
    double* du;
    Coupler_ptr<real_t> u_coupler;
    Coupler_ptr<real_t> u_old_coupler; /* last time step */
};


} // namespace paso

#endif // __PASO_FCTSOLVER_H__

