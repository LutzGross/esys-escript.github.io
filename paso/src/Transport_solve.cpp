
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
*
*****************************************************************************/


/****************************************************************************/

/* Paso: Transport solver
 *
 * solves Mu_t=Ku+q
 *
 * using an operator splitting approach (K=L+D where L is row sum zero and D
 * is a diagonal matrix):
 *
 *   - Mu_t=Du+q u(0)=u             (reactive part)
 *   - Mv_t=Lv   v(0)=u(dt/2)       (transport part uses flux correction (FCT))
 *   - Mu_t=Du+q u(dt/2)=v(dt/2)    (reactive part)
 *
 *  to return u(dt)
 *
 */
/****************************************************************************/

/* Author: l.gross@uq.edu.au */

/****************************************************************************/

#include "Transport.h"
#include "FCT_Solver.h"
#include "PasoUtil.h"
#include "ReactiveSolver.h"
#include "Solver.h"

#include <iostream>

namespace paso {

void TransportProblem::solve(double* u, double dt, double* u0, double* q,
                             Options* options)
{
    const real_t EPSILON = escript::DataTypes::real_t_eps();
    const real_t LARGE_POSITIVE_FLOAT = escript::DataTypes::real_t_max();
    const double reduction_after_divergence_factor = 0.5;
    const dim_t num_failures_max=50;

    Performance pp;
    ReactiveSolver* rsolver=NULL;
    FCT_Solver* fctsolver=NULL;

    dim_t i_substeps=0, n_substeps=1, num_failures=0;
    double *u_save=NULL, *u2=NULL;
    double  dt2,t=0, dt3;
    SolverResult errorCode=NoError;
    const dim_t n = transport_matrix->getTotalNumRows();
    options->time_step_backtracking_used = false;
    options->num_iter=0;

    if (dt <= 0.) {
        throw PasoException("TransportProblem::solve: dt must be positive.");
    } else if (getBlockSize() > 1) {
        throw PasoException("TransportProblem::solve: block size >1 "
                            "is not supported.");
    }
    if (options->verbose) {
        if (options->ode_solver == PASO_BACKWARD_EULER) {
            printf("TransportProblem::solve: Backward Euler is used (dt = %e)\n",dt);
        } else  if (options->ode_solver == PASO_LINEAR_CRANK_NICOLSON) {
            printf("TransportProblem::solve: linear Crank-Nicolson is used (dt = %e).\n",dt);
        } else  if (options->ode_solver == PASO_CRANK_NICOLSON) {
            printf("TransportProblem::solve: Crank-Nicolson is used (dt = %e).\n",dt);
        } else {
            throw PasoException("TransportProblem::solve: unknown ODE solver.");
        }
    }
    getSafeTimeStepSize();
    // allocate memory
    fctsolver = new FCT_Solver(shared_from_this(), options);
    rsolver = new ReactiveSolver(shared_from_this());
    u_save = new double[n];
    u2 = new double[n];

    // let the show begin!!!!
    const double dt_R = dt_max_R;
    const double dt_T = dt_max_T;
    dt2 = dt;
    if (dt_R < LARGE_POSITIVE_FLOAT)
        dt2 = std::min(dt_R*2, dt); // as we half the step size for the RT bit
    if (dt_T < LARGE_POSITIVE_FLOAT) {
        if (options->ode_solver == PASO_LINEAR_CRANK_NICOLSON || options->ode_solver == PASO_CRANK_NICOLSON) {
            dt2 = std::min(dt_T, dt);
        } // PASO_BACKWARD_EULER does not require a restriction
    }

    num_failures = 0;
    util::copy(n, u, u0); // copy initial value to return

    while((dt-t) > dt*sqrt(EPSILON)) {
        n_substeps = ceil((dt-t)/dt2);
        if (n_substeps <= 0) {
            throw PasoException("TransportProblem::solve: time stepping break down.");
        } else {
            dt3 = (dt-t)/n_substeps;
            if (options->verbose) {
                std::cout << "TransportProblem::solve: number of substeps = "
                    << n_substeps << " with dt = " << dt3 << "."
                    << std::endl;
            }
            // initialize the iteration matrix
            fctsolver->initialize(dt3, options, &pp);
            rsolver->initialize(dt3/2, options);
            errorCode = NoError;

            // start iteration
            for (i_substeps=0; i_substeps<n_substeps &&
                               errorCode==NoError; i_substeps++) {
                if (options->verbose) {
                    std::cout << "TransportProblem::solve: substep "
                        << i_substeps << " of " << n_substeps << " at t = "
                        << (t+dt3) << " (dt = " << dt3 << ")" << std::endl;
                }

                // create copy for restart in case of failure
                util::copy(n, u_save, u);
                // update u

                // Mu_t=Du+q u(0)=u
                errorCode = rsolver->solve(u2, u, q, options, &pp);

                // Mv_t=Lv   v(0)=u(dt/2)
                if (errorCode == NoError) {
                    errorCode = fctsolver->update(u, u2, options, &pp);

                }
                // Mu_t=Du+q u(dt/2)=v(dt/2)
                if (errorCode == NoError) {
                    errorCode = rsolver->solve(u2, u, q, options, &pp);
                }

                if (errorCode == NoError) {
                    num_failures = 0;
                    t += dt3;
                    util::copy(n, u, u2);
                }
            }
            if (errorCode == MaxIterReached || errorCode == Divergence) {
                // if num_failures_max failures in a row: give up
                if (num_failures >= num_failures_max) {
                    throw PasoException("TransportProblem::solve: "
                            "No convergence after time step reductions.");
                } else {
                    options->time_step_backtracking_used = true;
                    if (options->verbose) {
                        std::cout << "TransportProblem::solve: "
                            << "no convergence. Time step size is reduced."
                            << std::endl;
                    }
                    dt2 = dt3*reduction_after_divergence_factor;
                    num_failures++;
                    util::copy(n, u, u_save); // reset initial value
                }
            } else if (errorCode == InputError) {
                throw PasoException("TransportProblem::solve: input error for solver.");
            } else if (errorCode == MemoryError) {
                throw PasoException("TransportProblem::solve: memory allocation failed.");
            } else if (errorCode == Breakdown) {
                throw PasoException("TransportProblem::solve: solver break down.");
            } else if (errorCode == NegativeNormError) {
                throw PasoException("TransportProblem::solve: negative norm.");
            } else if (errorCode != NoError) {
                throw PasoException("TransportProblem::solve: general error.");
            }
        }
    } // end of time loop

    delete fctsolver;
    delete rsolver;
    delete[] u_save;
    delete[] u2;
}

double TransportProblem::getSafeTimeStepSize() const
{
    double dt_max=0.;
    const dim_t n = transport_matrix->getTotalNumRows();

    if (!valid_matrices) {
        // set row-sum of mass_matrix
        mass_matrix->rowSum(lumped_mass_matrix);
        // check for positive entries in lumped_mass_matrix and set
        // negative value for constraints
        int fail = 0;
#pragma omp parallel
        {
            int fail_loc = 0;
#pragma omp for
            for (index_t i=0; i<n; ++i) {
                const double m_i = lumped_mass_matrix[i];
                if (m_i > 0) {
                    if (constraint_mask[i] > 0)
                        lumped_mass_matrix[i]=-1.;
                } else {
                    fail_loc = 1;
                }
            }
            #pragma omp critical
            {
                fail = std::max(fail, fail_loc);
            }
        }
#ifdef ESYS_MPI
        int fail_loc = fail;
        MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, mpi_info->comm);
#endif
        if (fail > 0)
            throw PasoException("TransportProblem::getSafeTimeStepSize: "
                                "negative mass matrix entries detected.");
        // split off row-sum from transport_matrix
        transport_matrix->makeZeroRowSums(reactive_matrix);
        // get a copy of the main diagonal of the mass matrix
        mass_matrix->copyFromMainDiagonal(main_diagonal_mass_matrix);

        const double dt_R = ReactiveSolver::getSafeTimeStepSize(shared_from_this());
        const double dt_T = FCT_Solver::getSafeTimeStepSize(shared_from_this());
        dt_max_R = dt_R;
        dt_max_T = dt_T;
        valid_matrices = true;
        dt_max = std::min(2*dt_R, dt_T);
    } else {
        // factor 2 as we use operator splitting
        dt_max = std::min(2*dt_max_R, dt_max_T);
    }
    return dt_max;
}

} // namespace paso

