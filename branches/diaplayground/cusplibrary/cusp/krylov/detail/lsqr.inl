/*
 *  Copyright 2011 The Regents of the University of California
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/monitor.h>
#include <cusp/linear_operator.h>

namespace blas = cusp::blas;

namespace cusp {
namespace krylov {

template <class LinearOperator,
          class Vector,
          typename RealType>
lsqr_results<RealType> lsqr(LinearOperator& A,
                            Vector& x,
                            Vector& b,
                            const lsqr_parameters<RealType> & p)
{
    typedef typename LinearOperator::value_type ValueType;
    cusp::default_monitor<ValueType> monitor(b);
    return cusp::krylov::lsqr(A, x, b, p, monitor);
}

template <class LinearOperator,
          class Vector,
          typename RealType,
          class Monitor>
lsqr_results<RealType> lsqr(LinearOperator& A,
                            Vector& x,
                            Vector& b,
                            const lsqr_parameters<RealType> & p,
                            Monitor& monitor)
{
    CUSP_PROFILE_SCOPED();

    typedef typename LinearOperator::value_type   ValueType;
    typedef typename LinearOperator::memory_space MemorySpace;

    RealType damp = p.damp;
    const int n = A.num_cols;

    int istop  = 0;
    int nstop = 0;
    int nconv;
    RealType
        anorm  = 0,
        acond  = 0,
        rnorm  = 0,
        arnorm = 0,
        xnorm  = 0,
        atol   = p.atol,
        btol   = p.btol,
        conlim = p.conlim;
    const bool damped = damp > 0;
    RealType
        alpha, beta, bnorm,
        cs, cs1, cs2, ctol,
        delta, dknorm, dnorm,
        gamma, gambar, phi, phibar, psi,
        res2, rho, rhobar, rhbar1,
        rhs, rtol, sn, sn1, sn2,
        tau, temp, test1, test2, test3,
        theta, t1, t2, t3, xnorm1, z, zbar;

    rtol = monitor.relative_tolerance();
    ctol = (conlim > 0 ? 1.0 / conlim : 0.);
    dnorm  =  0;
    test2  =  0;
    res2   =  0;
    psi    =  0;
    xnorm1 =  0;
    cs2    = -1;
    sn2    =  0;
    z      =  0;

    //  ------------------------------------------------------------------
    //  Set up the first vectors u and v for the bidiagonalization.
    //  These satisfy  beta*u = b,  alpha*v = A^T*u.
    //  ------------------------------------------------------------------

    cusp::array1d<ValueType,MemorySpace> w(n);
    cusp::array1d<ValueType,MemorySpace> v(n, ValueType(0));
    cusp::array1d<ValueType,MemorySpace> tmp(n);
    blas::fill(x, ValueType(0));

    alpha = 0;
    beta  = blas::nrm2(b);

    if (beta > 0) {
        blas::scal(b, ValueType(1.0)/beta);
        //      v = v + At*b;
        cusp::transposed_multiply(A, b, tmp);
        blas::axpy(tmp, v, 1);
        alpha = blas::nrm2(v);
    }

    if (alpha > 0) {
        blas::scal(v, ValueType(1.0)/alpha);
        blas::copy(v, w);
    }
    cusp::array1d<ValueType,cusp::host_memory> residue(1);

    arnorm = alpha * beta;
    rhobar = alpha;
    phibar = beta;
    bnorm  = beta;
    rnorm  = beta;
    residue[0] = rnorm;

    // ==================================================================
    // Main iteration loop.
    // ==================================================================
    if (arnorm) {
        while (!monitor.finished(residue)) {
            if (istop != 0) {
                printf("istop = %d\n",istop);
                break;
            }

            ++monitor;
            // -----------------------------------------------------------
            // Perform the next step of the bidiagonalization to obtain
            // the next  beta, u, alpha, v.  These satisfy the relations
            //        beta*u  =  A*v    -  alpha*u,
            //       alpha*v  =  A^T*u  -   beta*v.
            // -----------------------------------------------------------
            blas::scal(b, -alpha);
            //  b = b + A*v;
            cusp::multiply(A, v, tmp);
            blas::axpy(tmp, b, 1);
            beta = blas::nrm2(b);

            // Accumulate anorm = || Bk ||
            //             =  sqrt( sum of  alpha**2 + beta**2 + damp**2 )

            temp   =   hypot( alpha, beta );
            temp   =   hypot( temp , damp );
            anorm  =   hypot( anorm, temp );

            if (beta > 0) {
                blas::scal(b, ValueType(1.0)/beta);
                blas::scal (v, -beta);
                //    v = v + At*b;
                cusp::transposed_multiply(A, b, tmp);
                blas::axpy(tmp, v, 1);
                //    aprod ( 2, m, n, v, u, UsrWrk );
                alpha = blas::nrm2(v);
                if (alpha > 0) {
                    blas::scal(v, ValueType(1.0)/alpha);
                }
            }

            // -----------------------------------------------------------
            // Use a plane rotation to eliminate the damping parameter.
            // This alters the diagonal (rhobar) of the lower-bidiagonal
            // matrix.
            // -----------------------------------------------------------
            rhbar1 = rhobar;
            if ( damped ) {
                rhbar1 = hypot( rhobar, damp );
                cs1    = rhobar / rhbar1;
                sn1    = damp   / rhbar1;
                psi    = sn1 * phibar;
                phibar = cs1 * phibar;
            }

            // -----------------------------------------------------------
            // Use a plane rotation to eliminate the subdiagonal element
            // (beta) of the lower-bidiagonal matrix, giving an
            // upper-bidiagonal matrix.
            // -----------------------------------------------------------
            rho    =   hypot( rhbar1, beta );
            cs     =   rhbar1 / rho;
            sn     =   beta   / rho;
            theta  =   sn * alpha;
            rhobar = - cs * alpha;
            phi    =   cs * phibar;
            phibar =   sn * phibar;
            tau    =   sn * phi;

            // -----------------------------------------------------------
            // Update  x, w  and (perhaps) the standard error estimates.
            // -----------------------------------------------------------
            t1     =   phi   / rho;
            t2     = - theta / rho;
            t3     =   1.0 / rho;
            dknorm =   0;

            blas::axpy(w, x, t1);
            dknorm = t3 * blas::nrm2(w);
            blas::axpby(v, w, w, 1, t2);
            // -----------------------------------------------------------
            // Monitor the norm of d_k, the update to x.
            // dknorm = norm( d_k )
            // dnorm  = norm( D_k ),  where   D_k = (d_1, d_2, ..., d_k )
            // dxk    = norm( phi_k d_k ),  where new x = x_k + phi_k d_k.
            // -----------------------------------------------------------
            dnorm  = hypot( dnorm, dknorm );
            //dxk    = abs( phi * dknorm );
            //if (dxmax < dxk) {
            //    dxmax = dxk;
            //}

            // -----------------------------------------------------------
            // Use a plane rotation on the right to eliminate the
            // super-diagonal element (theta) of the upper-bidiagonal
            // matrix. Then use the result to estimate  norm(x).
            // -----------------------------------------------------------
            delta  =   sn2 * rho;
            gambar = - cs2 * rho;
            rhs    =   phi    - delta * z;
            zbar   =   rhs    / gambar;
            xnorm  =   hypot( xnorm1, zbar  );
            gamma  =   hypot( gambar, theta );
            cs2    =   gambar / gamma;
            sn2    =   theta  / gamma;
            z      =   rhs    / gamma;
            xnorm1 =   hypot( xnorm1, z     );

            // -----------------------------------------------------------
            // Test for convergence.
            // First, estimate the norm and condition of the matrix  Abar,
            // and the norms of  rbar  and  Abar^T*rbar.
            // -----------------------------------------------------------
            acond  = anorm * dnorm;
            res2   = hypot( res2, psi    );
            rnorm  = hypot( res2, phibar );
            arnorm = alpha * abs( tau );

            // Now use these norms to estimate certain other quantities,
            // some of which will be small near a solution.

            test1 = rnorm / bnorm;
            test2 = 0;
            if (rnorm > 0) {
                test2 = arnorm / (anorm * rnorm);
            }

            test3 = ValueType(1.0) / acond;
            t1    = test1 / (ValueType(1.0) + anorm * xnorm / bnorm);
            if (btol || atol) {
                rtol = btol + atol * anorm * xnorm / bnorm;
            }

            residue[0] = rnorm;

            // The following tests guard against extremely small values of
            // atol, btol  or  ctol.  (The user may have set any or all of
            // the parameters  atol, btol, conlim  to zero.)
            // The effect is equivalent to the normal tests using
            // atol = relpr,  btol = relpr,  conlim = 1/relpr.

            t3 = 1.0 + test3;
            t2 = 1.0 + test2;
            t1 = 1.0 + t1;
            if (monitor.iteration_count() >= monitor.iteration_limit())
                istop = 5;
            if (t3 <= 1.0) istop = 4;
            if (t2 <= 1.0) istop = 2;
            if (t1 <= 1.0) istop = 1;

            // Allow for tolerances set by the user.

            if (test3 <= ctol) istop = 4;
            if (test2 <= atol) istop = 2;
            if (test1 <= rtol) istop = 1;

            if (istop == 0) {
                nstop = 0;
            } else {
                nconv = 1;
                nstop = nstop + 1;
                if (nstop < nconv && monitor.iteration_count() < monitor.iteration_limit()) {
                    istop = 0;
                }
            }
        }
    }
    // ==================================================================
    // End of iteration loop.
    // ==================================================================

    // Decide if istop = 2 or 3.
    // Print the stopping condition.
    if (damped  &&  istop == 2){
        istop = 3;
    }
    // Assign output variables from local copies.
    return lsqr_results<RealType>(istop, anorm, acond, rnorm, test2, xnorm);
}

} // end namespace krylov
} // end namespace cusp

