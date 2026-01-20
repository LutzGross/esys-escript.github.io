
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


/*
*
*  Purpose
*  =======
*
*  TFQMR solves the linear system A*x = b
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  r       (input) DOUBLE PRECISION array, dimension N.
*
*
*  x       (input/output) DOUBLE PRECISION array, dimension N.
*
*
*  ITER    (input/output) INT
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  INFO    (output) INT
*
*          = SOLVER_NO_ERROR: Successful exit. Iterated approximate solution returned.
*          = SOLVER_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVER_BREAKDOWN: If parameters RHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters RHO or OMEGA become smaller
*
*  ==============================================================
*/

#include "Solver.h"
#include "PasoUtil.h"

namespace paso {

//#define PASO_DYNAMIC_SCHEDULING_MVM
#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined __OPENMP
#define USE_DYNAMIC_SCHEDULING
#endif

SolverResult Solver_TFQMR(SystemMatrix_ptr<double> A, double* r, double* x, dim_t* iter,
                          double* tolerance, Performance* pp)
{
    int m=1;
    int j=0;
    dim_t num_iter=0;
    bool breakFlag=false, maxIterFlag=false, convergeFlag=false;
    SolverResult status = NoError;
    const dim_t n = A->getTotalNumRows();
    double eta,theta,tau,rho,beta,alpha,sigma,rhon,c;
    double norm_of_residual;

    double* u1 = new double[n];
    double* u2 = new double[n];
    double* y1 = new double[n];
    double* y2 = new double[n];
    double* d = new double[n];
    double* w = new double[n];
    double* v = new double[n];
    double* temp_vector = new double[n];
    double* res = new double[n];

    dim_t maxit = *iter;

    if (maxit <= 0) {
        status = InputError;
    }

    util::zeroes(n, x);

    Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER);
    A->solvePreconditioner(res, r);
    Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER);

    Performance_startMonitor(pp, PERFORMANCE_SOLVER);
    util::zeroes(n,u2);
    util::zeroes(n,y2);
    util::copy(n,w,res);
    util::copy(n,y1,res);
    util::zeroes(n,d);
    Performance_stopMonitor(pp, PERFORMANCE_SOLVER);

    Performance_startMonitor(pp, PERFORMANCE_MVM);
    A->MatrixVector_CSR_OFFSET0(PASO_ONE, y1, PASO_ZERO, temp_vector);
    Performance_stopMonitor(pp, PERFORMANCE_MVM);
    Performance_startMonitor(pp, PERFORMANCE_SOLVER);

    Performance_stopMonitor(pp, PERFORMANCE_SOLVER);
    Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER);
    A->solvePreconditioner(v,temp_vector);
    Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER);

    Performance_startMonitor(pp, PERFORMANCE_SOLVER);
    // v = P^{-1} * A y1
    util::copy(n, u1, v);

    theta = 0.0;
    eta = 0.0;
    tau = util::l2(n,res,A->mpi_info);
    rho = tau * tau;
    norm_of_residual=tau;

    while (!(convergeFlag || maxIterFlag || breakFlag || (status!=NoError) )) {
        sigma = util::innerProduct(n,res,v,A->mpi_info);
        if (! (breakFlag = (std::abs(sigma) == 0.))) {
            alpha = rho / sigma;
            for (j=0; j<=1; j=j+1) {
                // Compute y2 and u2 only if you have to
                if (j == 1) {
                    // y2 = y1 - alpha*v
                    util::linearCombination(n, y2, PASO_ONE, y1, -alpha, v);

                    Performance_stopMonitor(pp, PERFORMANCE_SOLVER);
                    Performance_startMonitor(pp, PERFORMANCE_MVM);
                    A->MatrixVector_CSR_OFFSET0(PASO_ONE, y2,PASO_ZERO,temp_vector);
                    Performance_stopMonitor(pp, PERFORMANCE_MVM);
                    Performance_startMonitor(pp, PERFORMANCE_SOLVER);

                    Performance_stopMonitor(pp, PERFORMANCE_SOLVER);
                    Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER);
                    A->solvePreconditioner(u2,temp_vector);
                    Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER);
                    Performance_startMonitor(pp, PERFORMANCE_SOLVER);
                    // u2 = P^{-1} * A y2
                }
                m = 2 * (num_iter+1) - 2 + (j+1);

                if (j==0) {
                    // w = w - alpha * u1
                    util::update(n, 1., w, -alpha, u1);
                    // d = (theta * theta * eta / alpha)*d + y1
                    util::update(n, (theta * theta * eta / alpha), d, 1., y1);
                } else if (j==1) {
                    // w = w - -alpha * u2
                    util::update(n, 1., w, -alpha, u2);
                    // d = (theta * theta * eta / alpha)*d + y2
                    util::update(n, (theta * theta * eta / alpha), d, 1., y2);
                }

                theta =util::l2(n,w,A->mpi_info)/tau;
                //printf("tau = %e, %e %e\n",tau, util::l2(n,w,A->mpi_info)/tau, theta);
                c = PASO_ONE / sqrt(PASO_ONE + theta * theta);
                tau = tau * theta * c;
                eta = c * c * alpha;
                util::update(n,1.,x,eta,d);
            }

            breakFlag = (std::abs(rho) == 0);

            rhon = util::innerProduct(n, res, w, A->mpi_info);
            beta = rhon / rho;
            rho = rhon;

            // y1 = w + beta * y2
            util::linearCombination(n,y1, PASO_ONE,w,beta,y2);

            Performance_stopMonitor(pp, PERFORMANCE_SOLVER);
            Performance_startMonitor(pp, PERFORMANCE_MVM);
            A->MatrixVector_CSR_OFFSET0(PASO_ONE, y1, PASO_ZERO, temp_vector);
            Performance_stopMonitor(pp, PERFORMANCE_MVM);

            Performance_startMonitor(pp, PERFORMANCE_PRECONDITIONER);
            A->solvePreconditioner(u1, temp_vector);
            Performance_stopMonitor(pp, PERFORMANCE_PRECONDITIONER);
            Performance_startMonitor(pp, PERFORMANCE_SOLVER);
            //  u1 = P^{-1} * A y1

            // t = u2 + beta * v
            util::linearCombination(n,temp_vector,PASO_ONE,u2,beta,v);
            // v = u1 + beta * t
            util::linearCombination(n,v,PASO_ONE,u1,beta,temp_vector);
        }
        maxIterFlag = (num_iter > maxit);
        norm_of_residual = tau*sqrt((double)(m + 1));
        convergeFlag = (norm_of_residual<(*tolerance));

        if (maxIterFlag) {
            status = MaxIterReached;
        } else if (breakFlag) {
            status = Breakdown;
        }
        ++num_iter;
    } // end of iterations

    Performance_stopMonitor(pp, PERFORMANCE_SOLVER);
    delete[] u1;
    delete[] u2;
    delete[] y1;
    delete[] y2;
    delete[] d;
    delete[] w;
    delete[] v;
    delete[] temp_vector;
    delete[] res;
    *iter=num_iter;
    *tolerance=norm_of_residual;
    return status;
}

} // namespace paso

