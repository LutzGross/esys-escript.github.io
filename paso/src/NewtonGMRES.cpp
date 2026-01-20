
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


/*
*  Purpose
*  =======
*
*  NewtonGMRES solves the nonlinear system F(x)=0 using the
*  restarted GMRES method.
*
*  Convergence test: norm(F(x))< rtol*norm(F(x0))+atol
*
*  x input initial guess
*    output solution approximation
*
*  Based on code by C. T. Kelley, July 1, 1994.
*
*/

#include "Paso.h"
#include "Options.h"
#include "PasoUtil.h"
#include "Solver.h"

#include <iostream>

namespace paso {

SolverResult Solver_NewtonGMRES(Function* F, double* x, Options* options,
                                Performance* pp)
{
    const double inner_tolerance_safety=.9;
    dim_t gmres_iter;
    double stop_tol, norm2_f,norm2_fo, normsup_f,reduction_f, gmres_tol, rtmp, quad_tolerance;
    bool convergeFlag=false, maxIterFlag=false, breakFlag=false;
    double *f=NULL, *step=NULL;
    SolverResult status=NoError;
    const bool debug = options->verbose;
    const dim_t n = F->getLen();
    dim_t iteration_count = 0;
    const double atol=options->absolute_tolerance;  /* absolute tolerance */
    const double rtol=options->tolerance;           /* relative tolerance */
    const dim_t maxit=options->iter_max;            /* max iteration counter */
    const dim_t lmaxit=options->inner_iter_max*10;  /* max inner iteration counter */
    const bool adapt_inner_tolerance=options->adapt_inner_tolerance;
    const double max_inner_tolerance=options->inner_tolerance *1.e-11;
    double inner_tolerance=max_inner_tolerance;
    /*
   * max_inner_tolerance = Maximum error tolerance for residual in inner
   *              iteration. The inner iteration terminates
   *              when the relative linear residual is smaller than
   *              inner_tolerance*| F(x_c) |. Inner_tolerance is determined
   *              by the modified Eisenstat-Walker formula, if
   *              adapt_inner_tolerance is set, otherwise
   *              inner_tolerance = max_inner_tolerance iteration.
    */

    f=new double[n];
    step=new double[n];
    /*
     * initial evaluation of F
     */
    F->call(f, x, pp);
    iteration_count++;
    norm2_f=util::l2(n,f,F->mpi_info);
    normsup_f=util::lsup(n,f,F->mpi_info);
    /*
     * stopping criterion:
     */
    stop_tol=atol + rtol*normsup_f;
    if (stop_tol<=0) {
        status=InputError;
        if (debug)
            std::cout << "NewtonGMRES: zero tolerance given." << std::endl;
    } else {
        iteration_count=1;
        if (debug) {
            std::cout << "NewtonGMRES: Start Jacobi-free Newton scheme" << std::endl
                << "NewtonGMRES: lsup tolerance rel/abs= "
                << rtol << "/" << atol << std::endl
                << "NewtonGMRES: lsup stopping tolerance = " << stop_tol
                << std::endl << "NewtonGMRES: max. inner iterations (GMRES) = "
                << lmaxit << std::endl;
            if (adapt_inner_tolerance) {
                std::cout << "NewtonGMRES: inner tolerance is adapted." << std::endl
                    << "NewtonGMRES: max. inner l2 tolerance (GMRES) = "
                    << max_inner_tolerance << std::endl;
            } else {
                std::cout << "NewtonGMRES: inner l2 tolerance (GMRES) = "
                    << inner_tolerance << std::endl;
            }
        }
        /*
         *  main iteration loop
         */
        while (! (convergeFlag || maxIterFlag || breakFlag)) {
            // keep track of the ratio (reduction_f = norm2_f/norm2_fo) of
            // successive residual norms and the iteration counter
            // (iteration_count)
            if (debug)
                std::cout << "NewtonGMRES: iteration step " << iteration_count
                    << ": lsup-norm of F = " << normsup_f << std::endl;

            // call GMRES to get increment
            gmres_iter=lmaxit;
            gmres_tol=inner_tolerance;
            status = Solver_GMRES2(F, f, x, step, &gmres_iter, &gmres_tol, pp);
            inner_tolerance=std::max(inner_tolerance, gmres_tol/norm2_f);
            std::cout << "NewtonGMRES: actual rel. inner tolerance = "
                << inner_tolerance << std::endl;
            iteration_count+=gmres_iter;
            if ((status==NoError) || (status==MaxIterReached)) {
                status=NoError;
                // update x
                norm2_fo=norm2_f;
                util::update(n,1.,x,1.,step);
                F->call(f, x, pp);
                iteration_count++;
                norm2_f=util::l2(n,f,F->mpi_info);
                normsup_f=util::lsup(n,f,F->mpi_info);
                reduction_f=norm2_f/norm2_fo;
                // adjust inner_tolerance
                if (adapt_inner_tolerance) {
                    quad_tolerance = inner_tolerance_safety * reduction_f * reduction_f;
                    rtmp=inner_tolerance_safety * inner_tolerance * inner_tolerance;
                    if (rtmp>.1)  {
                        inner_tolerance=std::min(max_inner_tolerance, std::max(quad_tolerance,rtmp));
                    } else {
                        inner_tolerance=std::min(max_inner_tolerance, quad_tolerance);
                    }
                    inner_tolerance=std::min(max_inner_tolerance, std::max(inner_tolerance, .5*stop_tol/normsup_f));
                }
                convergeFlag = (normsup_f <= stop_tol);
            } else {
                breakFlag = true;
            }
            maxIterFlag = (iteration_count > maxit);
        }
        if (debug) {
            std::cout << "NewtonGMRES: iteration step " << iteration_count
                << ": lsup-norm of F = " << normsup_f << std::endl;

            if (convergeFlag)
                std::cout << "NewtonGMRES: convergence reached after "
                    << iteration_count << " steps." << std::endl;
            if (breakFlag)
                std::cout << "NewtonGMRES: iteration break down after "
                    << iteration_count << " steps." << std::endl;
            if (maxIterFlag)
                std::cout << "NewtonGMRES: maximum number of iteration steps "
                    << maxit << " reached." << std::endl;
        }
        if (breakFlag) status = Breakdown;
        if (maxIterFlag) status = MaxIterReached;
    }
    delete[] f;
    delete[] step;
    if (debug)
        std::cout << "NewtonGMRES: STATUS return = " << status << std::endl;
    return status;
}

} // namespace paso

