
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

#include "Common.h"
#include "Solver.h"
#include "PasoUtil.h"

namespace paso {

err_t Solver_NewtonGMRES(Function* F, double* x, Options* options,
                         Paso_Performance* pp) 
{
    const double inner_tolerance_safety=.9;
    dim_t gmres_iter;
    double stop_tol, norm2_f,norm2_fo, normsup_f,reduction_f, gmres_tol, rtmp, quad_tolerance;
    bool convergeFlag=FALSE, maxIterFlag=FALSE, breakFlag=FALSE;
    double *f=NULL, *step=NULL;
    err_t Status=SOLVER_NO_ERROR;
    const bool debug=options->verbose;
    const dim_t n = F->getLen(); 
    dim_t iteration_count=0;
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
        Status=SOLVER_INPUT_ERROR;
        if (debug) printf("NewtonGMRES: zero tolerance given.\n");
    } else {
        iteration_count=1;
        if (debug) {
            printf("NewtonGMRES: Start Jacobi-free Newton scheme\n");
            printf("NewtonGMRES: lsup tolerance rel/abs= %e/%e\n",rtol,atol);
            printf("NewtonGMRES: lsup stopping tolerance = %e\n",stop_tol);
            printf("NewtonGMRES: max. inner iterations (GMRES) = %d\n",lmaxit);
            if (adapt_inner_tolerance) {
                printf("NewtonGMRES: inner tolerance is adapted.\n");
                printf("NewtonGMRES: max. inner l2 tolerance (GMRES) = %e\n",max_inner_tolerance);
            } else {
                printf("NewtonGMRES: inner l2 tolerance (GMRES) = %e\n",inner_tolerance);
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
                printf("NewtonGMRES: iteration step %d: lsup-norm of F =%e\n",iteration_count,normsup_f);

            // call GMRES to get increment
            gmres_iter=lmaxit;
            gmres_tol=inner_tolerance;
            Status = Solver_GMRES2(F,f,x,step,&gmres_iter,&gmres_tol,pp);
            inner_tolerance=MAX(inner_tolerance, gmres_tol/norm2_f);
            printf("NewtonGMRES: actual rel. inner tolerance = %e\n",inner_tolerance);
            iteration_count+=gmres_iter;
            if ((Status==SOLVER_NO_ERROR) || (Status==SOLVER_MAXITER_REACHED)) {
                Status=SOLVER_NO_ERROR;
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
                        inner_tolerance=MIN(max_inner_tolerance, MAX(quad_tolerance,rtmp));
                    } else {
                        inner_tolerance=MIN(max_inner_tolerance, quad_tolerance); 
                    }
                    inner_tolerance=MIN(max_inner_tolerance, MAX(inner_tolerance, .5*stop_tol/normsup_f));
                }
                convergeFlag = (normsup_f <= stop_tol);
            } else {
                breakFlag = true;
            }
            maxIterFlag = (iteration_count > maxit);
        }
        if (debug) {
            printf("NewtonGMRES: iteration step %d: lsup-norm of F =%e\n",iteration_count,normsup_f);

            if (convergeFlag) printf("NewtonGMRES: convergence reached after %d steps.\n",iteration_count);
            if (breakFlag) printf("NewtonGMRES: iteration break down after %d steps.\n",iteration_count);
            if (maxIterFlag) printf("NewtonGMRES: maximum number of iteration steps %d reached.\n",maxit);
        }
        if (breakFlag) Status = SOLVER_BREAKDOWN;
        if (maxIterFlag) Status = SOLVER_MAXITER_REACHED;
    }
    delete[] f;
    delete[] step;
    if (debug) printf("NewtonGMRES: STATUS return = %d\n",Status);
    return Status;
}

} // namespace paso

