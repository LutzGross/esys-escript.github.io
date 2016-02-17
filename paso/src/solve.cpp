
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

namespace paso {

void SystemMatrix::solve(double* out, double* in, Options* options) const
{
    Performance pp;
    index_t package;
    Esys_resetError();
    if (getGlobalNumCols() != getGlobalNumRows()
                    || col_block_size != row_block_size) {
        Esys_setError(VALUE_ERROR,"solve: matrix has to be a square matrix.");
        return;
    }
    //options->show();
    Performance_open(&pp, options->verbose);
    package = Options::getPackage(options->method, options->package, options->symmetric, mpi_info);
    if (Esys_noError()) {
        switch(package) {
            case PASO_PASO:
                Solver(boost::const_pointer_cast<SystemMatrix>(shared_from_this()),
                       out, in, options, &pp);
                solver_package = PASO_PASO;
            break;

            case PASO_MKL:
                if (mpi_info->size > 1) {
                    Esys_setError(VALUE_ERROR,"solve: MKL package does not support MPI.");
                    return;
                }
                options->converged = false;
                options->time = Esys_timer();
                Performance_startMonitor(&pp, PERFORMANCE_ALL);
                MKL_solve(mainBlock, out, in, options->reordering,
                          options->refinements, options->verbose);
                solver_package = PASO_MKL;
                Performance_stopMonitor(&pp, PERFORMANCE_ALL);
                options->time = Esys_timer()-options->time;
                options->set_up_time = 0;
                options->residual_norm = 0.;
                options->num_iter = 0;
                if (Esys_MPIInfo_noError(mpi_info))
                    options->converged = true;
            break;

            case PASO_UMFPACK:
                if (mpi_info->size > 1) {
                    Esys_setError(VALUE_ERROR,"solve: UMFPACK package does not support MPI.");
                    return;
                }
                options->converged = false;
                options->time = Esys_timer();
                Performance_startMonitor(&pp, PERFORMANCE_ALL);
                UMFPACK_solve(mainBlock, out, in, options->refinements, options->verbose);
                solver_package = PASO_UMFPACK;
                Performance_stopMonitor(&pp, PERFORMANCE_ALL);
                options->time = Esys_timer()-options->time;
                options->set_up_time = 0;
                options->residual_norm = 0.;
                options->num_iter = 0;
                if (Esys_MPIInfo_noError(mpi_info))
                    options->converged = true;
            break;

            default:
                Esys_setError(VALUE_ERROR, "solve: unknown package code");
            break;
        }
    }
    /*
        cancel divergence errors
    */
    if (options->accept_failed_convergence) {
        if (Esys_getErrorType() == DIVERGED) {
            Esys_resetError();
            if (options->verbose)
                printf("paso: failed convergence error has been canceled as requested.\n");
        }
    }
    Performance_close(&pp, options->verbose);
}

void solve_free(SystemMatrix* in)
{
    if (!in) return;

    switch(in->solver_package) {
        case PASO_PASO:
            Solver_free(in);
            break;

        case PASO_SMOOTHER:
            Preconditioner_Smoother_free((Preconditioner_Smoother*) in->solver_p);
            break;

        case PASO_MKL:
            MKL_free(in->mainBlock.get());
            break;

        case PASO_UMFPACK:
            UMFPACK_free(in->mainBlock.get());
            break;
   }
}

} // namespace paso

