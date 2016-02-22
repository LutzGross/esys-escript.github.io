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

/* Paso: SystemMatrix: controls iterative linear system solvers  */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003/04                     */
/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/****************************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "Options.h"
#include "SystemMatrix.h"

#include <boost/math/special_functions/fpclassify.hpp>  // for isnan

#include <iostream>

namespace bm = boost::math;

namespace paso {

void Solver_free(SystemMatrix* A)
{
    A->freePreconditioner();
}

///  calls the iterative solver
void Solver(SystemMatrix_ptr A, double* x, double* b, Options* options,
            Performance* pp)
{
    const real_t EPSILON = escript::DataTypes::real_t_eps();
    double norm2_of_b,tol,tolerance,time_iter,net_time_start;
    double *r=NULL,norm2_of_residual,last_norm2_of_residual,norm_max_of_b;
    double norm2_of_b_local,norm_max_of_b_local,norm2_of_residual_local;
    double norm_max_of_residual_local,norm_max_of_residual;
    double last_norm_max_of_residual;
#ifdef ESYS_MPI
    double loc_norm;
#endif
    dim_t i,totIter=0,cntIter,method;
    bool finalizeIteration;
    SolverResult errorCode = NoError;
    const dim_t numSol = A->getTotalNumCols();
    const dim_t numEqua = A->getTotalNumRows();
    double *x0=NULL;

    Esys_resetError();
    tolerance=options->tolerance;
    if (tolerance < 100.* EPSILON) {
        Esys_setError(VALUE_ERROR,"Solver: Tolerance is too small.");
    }
    if (tolerance >1.) {
        Esys_setError(VALUE_ERROR,"Solver: Tolerance must be less than one.");
    }
    method=Options::getSolver(options->method, PASO_PASO, options->symmetric, A->mpi_info);
    /* check matrix type */
    if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1) ) {
        Esys_setError(TYPE_ERROR,"Solver: Iterative solver requires CSR format with unsymmetric storage scheme and index offset 0.");
    }
    if (A->col_block_size != A->row_block_size) {
        Esys_setError(TYPE_ERROR,"Solver: Iterative solver requires row and column block sizes to be equal.");
    }
    if (A->getGlobalNumCols() != A->getGlobalNumRows()) {
        Esys_setError(TYPE_ERROR,"Solver: Iterative solver requires a square matrix.");
        return;
    }
    time_iter=Esys_timer();
    /* this for testing only */
    if (method==PASO_NONLINEAR_GMRES) {
        LinearSystem* F = new LinearSystem(A, b, options);
        A->solvePreconditioner(x, b);
        errorCode = Solver_NewtonGMRES(F, x, options, pp);
        if (errorCode != NoError) {
            Esys_setError(SYSTEM_ERROR,"Solver_NewtonGMRES: an error has occurred.");
        }
        delete F;
        return;
    }

    r = new double[numEqua];
    x0 = new double[numEqua];
    A->balance();
    options->num_level=0;
    options->num_inner_iter=0;

    /* ========================= */
    if (Esys_noError()) {
        Performance_startMonitor(pp, PERFORMANCE_ALL);
        A->applyBalance(r, b, true);
        /* get the norm of the right hand side */
        norm2_of_b=0.;
        norm_max_of_b=0.;
        #pragma omp parallel private(norm2_of_b_local,norm_max_of_b_local)
        {
            norm2_of_b_local=0.;
            norm_max_of_b_local=0.;
            #pragma omp for private(i) schedule(static)
            for (i = 0; i < numEqua ; ++i) {
                norm2_of_b_local += r[i] * r[i];
                norm_max_of_b_local = MAX(ABS(r[i]),norm_max_of_b_local);
            }
            #pragma omp critical
            {
                norm2_of_b += norm2_of_b_local;
                norm_max_of_b = MAX(norm_max_of_b_local,norm_max_of_b);
            }
        }
#ifdef ESYS_MPI
        /* TODO: use one call */
        loc_norm = norm2_of_b;
        MPI_Allreduce(&loc_norm,&norm2_of_b, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
        loc_norm = norm_max_of_b;
        MPI_Allreduce(&loc_norm,&norm_max_of_b, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
#endif
        norm2_of_b=sqrt(norm2_of_b);
        /* if norm2_of_b==0 we are ready: x=0 */
        if (bm::isnan(norm2_of_b) || bm::isnan(norm_max_of_b)) {
            Esys_setError(VALUE_ERROR, "Solver: Matrix or right hand side contains undefined values.");
        } else if (norm2_of_b <= 0.) {
#pragma omp parallel for private(i) schedule(static)
            for (i = 0; i < numSol; i++) x[i]=0.;
            if (options->verbose)
                std::cout << "right hand side is identical to zero." << std::endl;
        } else {
            if (options->verbose) {
                std::cout << "Solver: l2/lmax-norm of right hand side is "
                    << norm2_of_b << "/" << norm_max_of_b << "." << std::endl
                    << "Solver: l2/lmax-stopping criterion is "
                    << norm2_of_b*tolerance << "/" << norm_max_of_b*tolerance
                    << "." << std::endl;
                switch (method) {
                    case PASO_BICGSTAB:
                        std::cout << "Solver: Iterative method is BiCGStab.\n";
                    break;
                    case PASO_PCG:
                        std::cout << "Solver: Iterative method is PCG.\n";
                    break;
                    case PASO_TFQMR:
                        std::cout << "Solver: Iterative method is TFQMR.\n";
                    break;
                    case PASO_MINRES:
                        std::cout << "Solver: Iterative method is MINRES.\n";
                    break;
                    case PASO_PRES20:
                        std::cout << "Solver: Iterative method is PRES20.\n";
                    break;
                    case PASO_GMRES:
                        if (options->restart > 0) {
                            std::cout << "Solver: Iterative method is GMRES("
                                << options->truncation << ","
                                << options->restart << ")." << std::endl;
                        } else {
                            std::cout << "Solver: Iterative method is GMRES("
                                << options->truncation << ")." << std::endl;
                        }
                    break;
                }
            }

            // construct the preconditioner
            Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER_INIT);
            A->setPreconditioner(options);
            Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER_INIT);
            options->set_up_time=Esys_timer()-time_iter;
            if (Esys_noError()) {
                // get an initial guess by evaluating the preconditioner
                A->solvePreconditioner(x, r);

                totIter = 1;
                finalizeIteration = false;
                last_norm2_of_residual=norm2_of_b;
                last_norm_max_of_residual=norm_max_of_b;
                net_time_start=Esys_timer();

                // main loop
                while (!finalizeIteration) {
                    cntIter = options->iter_max - totIter;
                    finalizeIteration = true;

                    // Set initial residual
                    if (totIter > 1) {
                        // in the first iteration r = balance * b already
                        A->applyBalance(r, b, true);
                    }

                    A->MatrixVector_CSR_OFFSET0(-1., x, 1., r);
                    norm2_of_residual = 0;
                    norm_max_of_residual = 0;
                    #pragma omp parallel private(norm2_of_residual_local,norm_max_of_residual_local)
                    {
                        norm2_of_residual_local = 0;
                        norm_max_of_residual_local = 0;
                        #pragma omp for private(i) schedule(static)
                        for (i = 0; i < numEqua; i++) {
                            norm2_of_residual_local+= r[i] * r[i];
                            norm_max_of_residual_local=MAX(ABS(r[i]),norm_max_of_residual_local);
                        }
                        #pragma omp critical
                        {
                            norm2_of_residual += norm2_of_residual_local;
                            norm_max_of_residual = MAX(norm_max_of_residual_local,norm_max_of_residual);
                        }
                    }
#ifdef ESYS_MPI
                    // TODO: use one call
                    loc_norm = norm2_of_residual;
                    MPI_Allreduce(&loc_norm,&norm2_of_residual, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                    loc_norm = norm_max_of_residual;
                    MPI_Allreduce(&loc_norm,&norm_max_of_residual, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
#endif
                    norm2_of_residual =sqrt(norm2_of_residual);
                    options->residual_norm=norm2_of_residual;

                    if (options->verbose)
                        std::cout << "Solver: Step " << totIter
                            << ": l2/lmax-norm of residual is "
                            << norm2_of_residual << "/" << norm_max_of_residual;

                    if (totIter > 1 &&
                            norm2_of_residual >= last_norm2_of_residual &&
                            norm_max_of_residual >= last_norm_max_of_residual) {

                        if (options->verbose) std::cout << " divergence!\n";
                        Esys_setError(DIVERGED, "Solver: No improvement during iteration. Iterative solver gives up.");

                    } else {
                        if (norm2_of_residual>tolerance*norm2_of_b ||
                                norm_max_of_residual>tolerance*norm_max_of_b ) {

                            tol=tolerance*MIN(norm2_of_b,0.1*norm2_of_residual/norm_max_of_residual*norm_max_of_b);
                            if (options->verbose)
                                std::cout << " (new tolerance = " << tol << ").\n";

                            last_norm2_of_residual=norm2_of_residual;
                            last_norm_max_of_residual=norm_max_of_residual;

                            // call the solver
                            switch (method) {
                                case PASO_BICGSTAB:
                                    errorCode = Solver_BiCGStab(A, r, x, &cntIter, &tol, pp);
                                break;
                                case PASO_PCG:
                                    errorCode = Solver_PCG(A, r, x, &cntIter, &tol, pp);
                                break;
                                case PASO_TFQMR:
                                    tol=tolerance*norm2_of_residual/norm2_of_b;
                                    errorCode = Solver_TFQMR(A, r, x0, &cntIter, &tol, pp);
                                    #pragma omp for private(i) schedule(static)
                                    for (i = 0; i < numEqua; i++) {
                                        x[i]+= x0[i];
                                    }
                                break;
                                case PASO_MINRES:
                                    //tol=tolerance*norm2_of_residual/norm2_of_b;
                                    errorCode = Solver_MINRES(A, r, x, &cntIter, &tol, pp);
                                break;
                                case PASO_PRES20:
                                    errorCode = Solver_GMRES(A, r, x, &cntIter, &tol, 5, 20, pp);
                                break;
                                case PASO_GMRES:
                                    errorCode = Solver_GMRES(A, r, x, &cntIter, &tol, options->truncation, options->restart, pp);
                                break;
                            }

                            totIter += cntIter;

                            // error handling
                            if (errorCode == NoError) {
                                finalizeIteration = false;
                            } else if (errorCode == MaxIterReached) {
                                Esys_setError(DIVERGED, "Solver: maximum number of iteration steps reached.\nReturned solution does not fulfil stopping criterion.");
                                if (options->verbose)
                                    std::cout << "Solver: Maximum number of "
                                        "iterations reached." << std::endl;
                            } else if (errorCode == InputError) {
                                Esys_setError(SYSTEM_ERROR, "Solver: illegal dimension in iterative solver.");
                                if (options->verbose)
                                    std::cout << "Solver: Internal error!\n";
                            } else if (errorCode == NegativeNormError) {
                                Esys_setError(VALUE_ERROR, "Solver: negative energy norm (try other solver or preconditioner).");
                                if (options->verbose)
                                    std::cout << "Solver: negative energy norm"
                                       " (try other solver or preconditioner)!\n";
                            } else if (errorCode == Breakdown) {
                                if (cntIter <= 1) {
                                    Esys_setError(ZERO_DIVISION_ERROR, "Solver: fatal break down in iterative solver.");
                                    if (options->verbose)
                                        std::cout << "Solver: Uncurable break "
                                            "down!" << std::endl;
                                } else {
                                    if (options->verbose)
                                        std::cout << "Solver: Breakdown at iter "
                                            << totIter << " (residual = "
                                            << tol << "). Restarting ...\n";
                                    finalizeIteration = false;
                                }
                            } else {
                                Esys_setError(SYSTEM_ERROR, "Solver: Generic error in solver.");
                                if (options->verbose)
                                    std::cout << "Solver: Generic error in solver!\n";
                            }
                        } else {
                            if (options->verbose)
                                std::cout << " convergence!" << std::endl;
                            options->converged = true;
                        }
                    }
                } // while
                options->net_time = Esys_timer()-net_time_start;
            }
            options->num_iter = totIter;
            A->applyBalanceInPlace(x, false);
        }
    }
    delete[] r;
    delete[] x0;
    options->time = Esys_timer()-time_iter;
    Performance_stopMonitor(pp, PERFORMANCE_ALL);
}

} // namespace paso

