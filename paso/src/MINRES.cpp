
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


/* MINRES iterations */

#include "Solver.h"
#include "PasoUtil.h"
#include "SystemMatrix.h"

namespace paso {

/*
*
*  Purpose
*  =======
*
*  MINRES solves the linear system A*x = b
*
*  Convergence test: norm( b - A*x )< TOL.
*
*  Arguments
*  =========
*
*  R      (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of initial guess x
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess.
*
*  ITER    (input/output) INT
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  INFO    (output) INT
*
*          = NoError: Successful exit. Iterated approximate solution returned.
*          = MaxIterReached
*          = InputError Illegal parameter:
*          = MemoryError :
*          = NegativeNormError :
*
*  ==============================================================
*/

SolverResult Solver_MINRES(SystemMatrix_ptr<double> A, double* R, double* X,
                           dim_t* iter, double* tolerance, Performance* pp)
{
    const dim_t maxit = *iter;
    if (maxit <= 0) {
        return InputError;
    }

    double delta,gamma=0.,gamma_old=0.,eta=0.,dp0=0., c=1.0,c_old=1.0;
    double c_ancient=1.,s=0.0,s_old=0.0,s_ancient, norm_of_residual=0., rnorm_prec=1;
    double tol=1., norm_scal=1.;
    double alpha_0,alpha_1,alpha_2,alpha_3,dp = 0.0;
    dim_t num_iter = 0;
    const dim_t n = A->getTotalNumRows();
    bool convergeFlag=false;
    SolverResult status = NoError;

    double* ZNEW = new double[n];
    double* Z = new double[n];
    double* AZ = new double[n];
    double* W = new double[n];
    double* R_old = new double[n];
    double* W_old = new double[n];
    double* R_ancient = new double[n];
    double* W_ancient = new double[n];

    // z  <- Prec*r
    A->solvePreconditioner(Z, R);
    // gamma <- r'*z
    dp = util::innerProduct(n, R ,Z,A->mpi_info);
    dp0 = dp;
    if (dp < 0) {
        status = NegativeNormError;
    } else if (std::abs(dp) <= 0) {
        // happy break down
        convergeFlag = true;
    } else {
        // gamma <- sqrt(r'*z)
        gamma = sqrt(dp);
        eta = gamma;
        rnorm_prec = gamma;
        norm_of_residual=util::l2(n, R, A->mpi_info);
        norm_scal=rnorm_prec/norm_of_residual;
        tol=(*tolerance)*norm_scal;
    }

    while (!convergeFlag && status == NoError) {
        //  z <- z / gamma
        util::scale(n, Z, 1./gamma);

        //  Az <- A*z
        A->MatrixVector_CSR_OFFSET0(PASO_ONE, Z, PASO_ZERO, AZ);

        //  delta <- Az'.z
        delta = util::innerProduct(n, AZ, Z, A->mpi_info);

        //  r_new <- Az-delta/gamma * r - gamma/gamma_old r_old
        if (num_iter>0)
            util::copy(n, R_ancient, R_old); //  r__ancient <- r_old

        util::copy(n, R_old, R); //  r_old <- r
        util::copy(n, R, AZ);    //  r <- Az
        util::AXPY(n, R, -delta/gamma, R_old); // r <- r - delta/gamma v
        if (num_iter > 0)
            util::AXPY(n, R, -gamma/gamma_old, R_ancient); // r <- r - gamma/gamma_old r__ancient

        //  z <- prec*r
        A->solvePreconditioner(ZNEW, R);

        dp = util::innerProduct(n, R, ZNEW, A->mpi_info);
        if (dp < 0.) {
            status = NegativeNormError;
        } else if (std::abs(dp) == 0.) {
            // happy break down
            convergeFlag = true;
        } else if (std::abs(dp) > 0.e-13 * std::abs(dp0)) {
            //  gamma <- sqrt(r'*z)
            gamma_old = gamma;
            gamma = sqrt(dp);
            // QR factorisation
            c_ancient = c_old; c_old = c;
            s_ancient = s_old; s_old = s;

            alpha_0 = c_old * delta - c_ancient * s_old * gamma_old;
            alpha_1 = sqrt(alpha_0*alpha_0 + gamma*gamma);
            alpha_2 = s_old * delta + c_ancient * c_old * gamma_old;
            alpha_3 = s_ancient * gamma_old;

            // Givens rotation
            c = alpha_0 / alpha_1;
            s = gamma / alpha_1;

            rnorm_prec = rnorm_prec * s;

            // w_new <- (z-alpha_3 w - alpha_2 w_old)/alpha_1

            if (num_iter > 1)
                util::copy(n, W_ancient, W_old); //  w__ancient <- w_old
            if (num_iter > 0)
                util::copy(n, W_old, W); //  w_old <- w

            util::copy(n, W, Z);
            if (num_iter > 1)
                util::AXPY(n, W, -alpha_3, W_ancient); // w <- w - alpha_3 w__ancient
            if (num_iter > 0)
                util::AXPY(n, W, -alpha_2, W_old); // w <- w - alpha_2 w_old
            util::scale(n, W, 1.0 / alpha_1);      // w <- w / alpha_1

            util::AXPY(n, X, c * eta, W);          // x <- x + c eta w
            eta = - s * eta;
            convergeFlag = rnorm_prec <= tol;
        } else {
            status = Breakdown;
        }
        util::copy(n, Z, ZNEW);
        ++num_iter;
        if (!convergeFlag && num_iter >= maxit)
            status = MaxIterReached;
    }
    delete[] Z;
    delete[] ZNEW;
    delete[] AZ;
    delete[] R_old;
    delete[] R_ancient;
    delete[] W;
    delete[] W_old;
    delete[] W_ancient;

    *iter=num_iter;
    *tolerance=rnorm_prec/norm_scal;
    return status;
}

} // namespace paso

