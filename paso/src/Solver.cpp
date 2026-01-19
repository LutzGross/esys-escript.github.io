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

/* Paso: SystemMatrix: controls iterative linear system solvers  */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003/04                     */
/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/****************************************************************************/

#include "Solver.h"
#include "Options.h"
#include "SystemMatrix.h"

#include <boost/math/special_functions/fpclassify.hpp>  // for isnan

#include <iostream>

namespace bm = boost::math;

namespace paso {

void Solver_free(SystemMatrix<double>* A)
{
    A->freePreconditioner();
}

///  calls the iterative solver
SolverResult Solver(SystemMatrix_ptr<double> A, double* x, double* b, Options* options,
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

    tolerance=options->tolerance;
    if (tolerance < 100.* EPSILON) {
        throw PasoException("Solver: Tolerance is too small.");
    }
    if (tolerance >1.) {
        throw PasoException("Solver: Tolerance must be less than one.");
    }
    method=Options::getSolver(options->method, PASO_PASO, options->symmetric, A->mpi_info);
    /* check matrix type */
    if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1) ) {
        throw PasoException("Solver: Iterative solver requires CSR format with unsymmetric storage scheme and index offset 0.");
    }
    if (A->col_block_size != A->row_block_size) {
        throw PasoException("Solver: Iterative solver requires row and column block sizes to be equal.");
    }
    if (A->getGlobalNumCols() != A->getGlobalNumRows()) {
        throw PasoException("Solver: Iterative solver requires a square matrix.");
    }
    time_iter=escript::gettime();
    /* this for testing only */
    if (method==PASO_NONLINEAR_GMRES) {
        LinearSystem* F = new LinearSystem(A, b, options);
        A->solvePreconditioner(x, b);
        errorCode = Solver_NewtonGMRES(F, x, options, pp);
        if (errorCode != NoError) {
            throw PasoException("Solver_NewtonGMRES: an error has occurred.");
        }
        delete F;
        return errorCode;
    }

    r = new double[numEqua];
    x0 = new double[numEqua];
    A->balance();
    options->num_level=0;
    options->num_inner_iter=0;

    /* ========================= */
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
            norm_max_of_b_local = std::max(std::abs(r[i]),norm_max_of_b_local);
        }
        #pragma omp critical
        {
            norm2_of_b += norm2_of_b_local;
            norm_max_of_b = std::max(norm_max_of_b_local,norm_max_of_b);
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
        throw PasoException("Solver: Matrix or right hand side contains undefined values.");
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
        options->set_up_time=escript::gettime()-time_iter;
        // get an initial guess by evaluating the preconditioner
        A->solvePreconditioner(x, r);

        totIter = 1;
        finalizeIteration = false;
        last_norm2_of_residual=norm2_of_b;
        last_norm_max_of_residual=norm_max_of_b;
        net_time_start=escript::gettime();

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
                    norm_max_of_residual_local=std::max(std::abs(r[i]),norm_max_of_residual_local);
                }
                #pragma omp critical
                {
                    norm2_of_residual += norm2_of_residual_local;
                    norm_max_of_residual = std::max(norm_max_of_residual_local,norm_max_of_residual);
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
                throw PasoException("Solver: No improvement during iteration. Iterative solver gives up.");

            } else {
                if (norm2_of_residual>tolerance*norm2_of_b ||
                        norm_max_of_residual>tolerance*norm_max_of_b ) {

                    tol=tolerance*std::min(norm2_of_b,0.1*norm2_of_residual/norm_max_of_residual*norm_max_of_b);
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
                        if (options->verbose)
                            std::cout << "Solver: Maximum number of "
                                "iterations reached." << std::endl;
                        break;
                    } else if (errorCode == InputError) {
                        if (options->verbose)
                            std::cout << "Solver: Internal error!\n";
                        break;
                    } else if (errorCode == NegativeNormError) {
                        if (options->verbose)
                            std::cout << "Solver: negative energy norm"
                               " (try other solver or preconditioner)!\n";
                        break;
                    } else if (errorCode == Breakdown) {
                        if (cntIter <= 1) {
                            if (options->verbose)
                                std::cout << "Solver: Uncurable break "
                                    "down!" << std::endl;
                            break;
                        } else {
                            if (options->verbose)
                                std::cout << "Solver: Breakdown at iter "
                                    << totIter << " (residual = "
                                    << tol << "). Restarting ...\n";
                            finalizeIteration = false;
                            errorCode = NoError;
                        }
                    } else {
                        if (options->verbose)
                            std::cout << "Solver: Generic error in solver!\n";
                        break;
                    }
                } else {
                    if (options->verbose)
                        std::cout << " convergence!" << std::endl;
                    options->converged = true;
                }
            }
        } // while
        options->net_time = escript::gettime()-net_time_start;
        options->num_iter = totIter;
        A->applyBalanceInPlace(x, false);
    }
    delete[] r;
    delete[] x0;
    options->time = escript::gettime()-time_iter;
    Performance_stopMonitor(pp, PERFORMANCE_ALL);
    return errorCode;
}

SolverResult Solver(SystemMatrix_ptr<cplx_t> A, cplx_t* x, cplx_t* b, Options* options,
                    Performance* pp)
{
    throw PasoException("Solver(): complex not implemented.");
}

void Solver_free(SystemMatrix<cplx_t>* A)
{
    throw PasoException("Solver_free(): complex not implemented.");
}

} // namespace paso

