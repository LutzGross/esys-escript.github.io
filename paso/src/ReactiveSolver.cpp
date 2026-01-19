
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: reactive solver (D is a diagonal matrix)
 *
 *   - Mv_t=Dv+q   v(0)=u
 *
 *  to return v(dt)
 *
*/

/****************************************************************************/

/* Author: l.gross@uq.edu.au                                                */

/****************************************************************************/

#include "ReactiveSolver.h"
#include "PasoUtil.h"
#include "Solver.h"

namespace paso {

static const real_t EPSILON = escript::DataTypes::real_t_eps();

// exp(h)-1 ~ h + h**2/2 for abs(h) <  PASO_RT_EXP_LIM_MIN
static const real_t PASO_RT_EXP_LIM_MIN = sqrt(EPSILON);

// it is assumed that exp(h) with  h>PASO_RT_EXP_LIM_MAX is not reliable
static const real_t PASO_RT_EXP_LIM_MAX = log(1./sqrt(EPSILON));

SolverResult ReactiveSolver::solve(double* u, double* u_old,
                                   const double* source, Options* options,
                                   Performance* pp)
{
    const double EXP_LIM_MIN = PASO_RT_EXP_LIM_MIN;
    const double EXP_LIM_MAX = PASO_RT_EXP_LIM_MAX;
    const dim_t n = tp->transport_matrix->getTotalNumRows();
    int fail = 0;

#pragma omp parallel for
    for (dim_t i=0; i<n; ++i) {
        const double m_i = tp->lumped_mass_matrix[i];
        if (m_i > 0) {
            const double d_ii = tp->reactive_matrix[i];
            const double x_i = dt*d_ii/m_i;
            if (x_i >= EXP_LIM_MAX) {
                fail = 1;
            } else {
                const double F_i = source[i];
                const double e_i = exp(x_i);
                double u_i = e_i*u_old[i];
                if (std::abs(x_i) > EXP_LIM_MIN) {
                    u_i += F_i/d_ii*(e_i-1.);
                } else {
                    // second order approximation of (exp(x_i)-1)/x_i
                    u_i += F_i*dt/m_i * (1. + x_i/2);
                }
                u[i] = u_i;
            }
        } else {
            u[i] = u_old[i] + dt * source[i]; // constraints added
        }
    }
#ifdef ESYS_MPI
    index_t fail_loc = fail;
    MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, tp->mpi_info->comm);
#endif
    if (fail > 0) {
        return Divergence;
    } else {
        return NoError;
    }
}

double ReactiveSolver::getSafeTimeStepSize(const_TransportProblem_ptr tp)
{
    const real_t LARGE_POSITIVE_FLOAT = escript::DataTypes::real_t_max();
    const double EXP_LIM_MAX = PASO_RT_EXP_LIM_MAX;
    const dim_t n = tp->transport_matrix->getTotalNumRows();
    double dt_max = LARGE_POSITIVE_FLOAT;

    // calculate time step size
#pragma omp parallel
    {
        double dt_max_loc = LARGE_POSITIVE_FLOAT;
#pragma omp for
        for (dim_t i=0; i<n; ++i) {
            const double d_ii = tp->reactive_matrix[i];
            const double m_i = tp->lumped_mass_matrix[i];
            if (m_i > 0) { // no constraint
                if (d_ii > 0)
                    dt_max_loc = std::min(dt_max_loc, m_i/d_ii);
            }
        }
        #pragma omp critical
        {
            dt_max = std::min(dt_max, dt_max_loc);
        }
    }
#ifdef ESYS_MPI
    double dt_max_loc = dt_max;
    MPI_Allreduce(&dt_max_loc, &dt_max, 1, MPI_DOUBLE, MPI_MIN, tp->mpi_info->comm);
#endif

    if (dt_max < LARGE_POSITIVE_FLOAT ) {
        dt_max *= 0.5*EXP_LIM_MAX; // make sure there is no exp overflow
    } else {
        dt_max = LARGE_POSITIVE_FLOAT;
    }
    return dt_max;
}

} // namespace paso

