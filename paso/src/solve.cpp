
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
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to the direct solvers                      */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003                        */
/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/****************************************************************************/

#include "Paso.h"
#include "Options.h"
#include "performance.h"
#include "Preconditioner.h"
#include "Solver.h"
#include "MKL.h"
#include "UMFPACK.h"
#include "MUMPS.h"

namespace paso {

template <>
void SystemMatrix<double>::solve(double* out, double* in, Options* options) const
{
    Performance pp;
    index_t package;
    if (getGlobalNumCols() != getGlobalNumRows()
                    || col_block_size != row_block_size) {
        throw PasoException("solve: matrix has to be a square matrix.");
    }
    //options->show();
    Performance_open(&pp, options->verbose);
    package = Options::getPackage(options->method, options->package, options->symmetric, mpi_info);
    SolverResult res = NoError;

    switch (package) {
        case PASO_PASO:
            res = Solver(boost::const_pointer_cast<SystemMatrix>(
                    boost::dynamic_pointer_cast<const SystemMatrix>(getPtr())),
                    out, in, options, &pp);
            solver_package = PASO_PASO;
        break;

        case PASO_MKL:
            if (mpi_info->size > 1) {
                throw PasoException("solve: MKL package does not support MPI.");
            }
            options->converged = false;
            options->time = escript::gettime();
            Performance_startMonitor(&pp, PERFORMANCE_ALL);
            MKL_solve(mainBlock, out, in, options->reordering,
                      options->refinements, options->verbose);
            solver_package = PASO_MKL;
            Performance_stopMonitor(&pp, PERFORMANCE_ALL);
            options->time = escript::gettime()-options->time;
            options->set_up_time = 0;
            options->residual_norm = 0.;
            options->num_iter = 0;
            options->converged = true;
        break;

        case PASO_UMFPACK:
            if (mpi_info->size > 1) {
                throw PasoException("solve: UMFPACK package does not support MPI.");
            }
            options->converged = false;
            options->time = escript::gettime();
            Performance_startMonitor(&pp, PERFORMANCE_ALL);
            UMFPACK_solve(mainBlock, out, in, options->refinements, options->verbose);
            solver_package = PASO_UMFPACK;
            Performance_stopMonitor(&pp, PERFORMANCE_ALL);
            options->time = escript::gettime()-options->time;
            options->set_up_time = 0;
            options->residual_norm = 0.;
            options->num_iter = 0;
            options->converged = true;
        break;

        case PASO_MUMPS:
            if (mpi_info->size > 1) {
                throw PasoException("solve: MUMPS support for single MPI rank only.");
            }
            options->converged = false;
            options->time = escript::gettime();
            Performance_startMonitor(&pp, PERFORMANCE_ALL);
            MUMPS_solve(mainBlock, out, in, options->refinements, options->verbose, mpi_info);
            solver_package = PASO_MUMPS;
            Performance_stopMonitor(&pp, PERFORMANCE_ALL);
            options->time = escript::gettime()-options->time;
            options->set_up_time = 0;
            options->residual_norm = 0.;
            options->num_iter = 0;
            options->converged = true;
        break;

        default:
            throw PasoException("solve: unknown package code");
        break;
    }

    if (res == Divergence) {
        // cancel divergence errors
        if (options->accept_failed_convergence) {
            if (options->verbose)
                printf("paso: failed convergence error has been canceled as requested.\n");
        } else {
            throw PasoException("Solver: No improvement during iteration. Iterative solver gives up.");
        }
    } else if (res == MaxIterReached) {
        // cancel divergence errors
        if (options->accept_failed_convergence) {
            if (options->verbose)
                printf("paso: failed convergence error has been canceled as requested.\n");
        } else {
            throw PasoException("Solver: maximum number of iteration steps reached.\nReturned solution does not fulfil stopping criterion.");
        }
    } else if (res == InputError) {
        throw PasoException("Solver: illegal dimension in iterative solver.");
    } else if (res == NegativeNormError) {
        throw PasoException("Solver: negative energy norm (try other solver or preconditioner).");
    } else if (res == Breakdown) {
        throw PasoException("Solver: fatal break down in iterative solver.");
    } else if (res != NoError) {
        throw PasoException("Solver: Generic error in solver.");
    }
    Performance_close(&pp, options->verbose);
}

template <>
void SystemMatrix<cplx_t>::solve(cplx_t* out, cplx_t* in, Options* options) const
{
    Performance pp;
    index_t package;
    if (getGlobalNumCols() != getGlobalNumRows()
                    || col_block_size != row_block_size) {
        throw PasoException("solve: matrix has to be a square matrix.");
    }
    //options->show();
    Performance_open(&pp, options->verbose);
    package = Options::getPackage(options->method, options->package, options->symmetric, mpi_info);
    SolverResult res = NoError;

    switch (package) {
        case PASO_MUMPS:
            if (mpi_info->size > 1) {
                throw PasoException("solve: MUMPS support for single MPI rank only.");
            }
            options->converged = false;
            options->time = escript::gettime();
            Performance_startMonitor(&pp, PERFORMANCE_ALL);
            MUMPS_solve(mainBlock, out, in, options->refinements, options->verbose, mpi_info);
            solver_package = PASO_MUMPS;
            Performance_stopMonitor(&pp, PERFORMANCE_ALL);
            options->time = escript::gettime()-options->time;
            options->set_up_time = 0;
            options->residual_norm = 0.;
            options->num_iter = 0;
            options->converged = true;
        break;

        default:
            throw PasoException("solve: MUMPS required for complex matrices.");
        break;
    }

    if (res == Divergence) {
        // cancel divergence errors
        if (options->accept_failed_convergence) {
            if (options->verbose)
                printf("paso: failed convergence error has been canceled as requested.\n");
        } else {
            throw PasoException("Solver: No improvement during iteration. Iterative solver gives up.");
        }
    } else if (res == MaxIterReached) {
        // cancel divergence errors
        if (options->accept_failed_convergence) {
            if (options->verbose)
                printf("paso: failed convergence error has been canceled as requested.\n");
        } else {
            throw PasoException("Solver: maximum number of iteration steps reached.\nReturned solution does not fulfil stopping criterion.");
        }
    } else if (res == InputError) {
        throw PasoException("Solver: illegal dimension in iterative solver.");
    } else if (res == NegativeNormError) {
        throw PasoException("Solver: negative energy norm (try other solver or preconditioner).");
    } else if (res == Breakdown) {
        throw PasoException("Solver: fatal break down in iterative solver.");
    } else if (res != NoError) {
        throw PasoException("Solver: Generic error in solver.");
    }
    Performance_close(&pp, options->verbose);
}

} // namespace paso

