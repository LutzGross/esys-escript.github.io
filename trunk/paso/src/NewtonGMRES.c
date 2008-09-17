/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/*
*  Purpose
*  =======
*
*  NewtonGMRES solves the non-linear system F(x)=0 using the
*  restarted GMRES method.
*
*  Convergence test: norm(F(x))< rtol*norm(F(x0))+atol
*
*  x input initial guess
*    output solution approximation
*
*  based on code by C. T. Kelley, July 1, 1994.
*
*/

#include "Common.h"
#include "Solver.h"
#include "PasoUtil.h"
#ifdef _OPENMP
#include <omp.h>
#endif
err_t Paso_Solver_NewtonGMRES(
    Paso_Function *F,    /* function evaluation */
    double *x,           /* in: initial guess, out: new approximation */
    Paso_Options* options,
    Paso_Performance* pp) 

{
   const double inner_tolerance_safety=.9;
   dim_t gmres_iter;
   double stop_tol, norm_f,norm_fo, reduction_f,old_inner_tolerance, gmres_tol, rtmp;
   bool_t convergeFlag=FALSE, maxIterFlag=FALSE, breakFlag=FALSE;
   double *f=NULL, *step=NULL;
   err_t Status=SOLVER_NO_ERROR;
   const bool_t debug=options->verbose;
   const dim_t n = F->n; 
   dim_t iteration_count=0;
   const double atol=options->absolute_tolerance;  /* absolute tolerance */
   const double rtol=options->tolerance;           /* relative tolerance */
   const dim_t maxit=options->iter_max;            /* max iteration counter */
   const dim_t lmaxit=options->inner_iter_max;     /* max inner iteration counter */
   const bool_t adapt_inner_tolerance=options->adapt_inner_tolerance;
   const double max_inner_tolerance=options->inner_tolerance;  /*inner tolerance counter */
   double inner_tolerance=max_inner_tolerance;
  /*
   * max_inner_tolerance = Maximum error tolerance for residual in inner
   *              iteration. The inner iteration terminates
   *              when the relative linear residual is smaller than inner_tolerance*| F(x_c) |. inner_tolerance is determined
   *              by the modified Eisenstat-Walker formula, if adapt_inner_tolerance is set, otherwise 
   *              inner_tolerance = max_inner_tolerance iteration.
   */
  
  f=TMPMEMALLOC(n,double);
  step=TMPMEMALLOC(n,double);
  if (f==NULL || step ==NULL) {
      Status=SOLVER_MEMORY_ERROR;
  } else {
     /* 
      * initial evaluation of F
      */
     Paso_FunctionCall(F,f,x);
     norm_f=Paso_l2(n,f,F->mpi_info);
     /*
      * stoping criterium:
      */
     stop_tol=atol + rtol*norm_f;
     iteration_count=1;
     if (debug) printf("Start Jacobi-free Newton scheme:\n\ttolerance = %e\n\tstopping tolerance = %e\n",rtol,stop_tol);
     /* 
      *  main iteration loop
      */
     while (! (convergeFlag || maxIterFlag || breakFlag)) {
         /* 
          * keep track of the ratio (reduction_f = norm_f/frnmo) of successive residual norms and 
          * the iteration counter (iteration_count)
          */
         if (debug) printf("iteration step %d: norm of F =%lg\n",iteration_count,norm_f);
         /*
          * call GMRES to get increment
          */
         gmres_iter=lmaxit;
         gmres_tol=inner_tolerance;
         if (debug) printf("GMRES called with tolerance = %e (max iter=%d)\n",inner_tolerance,gmres_iter);
         Status=Paso_Solver_GMRES2(F,f,x,step,&gmres_iter,&gmres_tol,pp);
         if (debug) printf("GMRES finalized after %d steps (residual = %e)\n",gmres_iter,gmres_tol);
         iteration_count+=gmres_iter;
         if ((Status==SOLVER_NO_ERROR) || (Status==SOLVER_MAXITER_REACHED)) {
            Status=SOLVER_NO_ERROR;
            /* 
             * update x:
             */
            norm_fo=norm_f; 
            Paso_Update(n,1.,x,1.,step);
            Paso_FunctionCall(F,f,x);
            iteration_count++;
            norm_f=Paso_l2(n,f,F->mpi_info);
            reduction_f=norm_f/norm_fo;
            /*
             *   adjust inner_tolerance 
             */
            if (adapt_inner_tolerance) {
                 old_inner_tolerance=inner_tolerance;
                 inner_tolerance=inner_tolerance_safety * reduction_f * reduction_f;
                 rtmp=inner_tolerance_safety * old_inner_tolerance * old_inner_tolerance;
                 if (rtmp>.1) inner_tolerance=MAX(inner_tolerance,rtmp);
                 inner_tolerance=MAX(MIN(inner_tolerance,max_inner_tolerance), .5*stop_tol/norm_f);
            }
            convergeFlag = (norm_f <= stop_tol);
         } else {
            breakFlag=TRUE;
         }
         maxIterFlag = (iteration_count > maxit);
      }
      if (debug) {
           if (convergeFlag) printf("convergence reached after %d steps with residual %e.\n",iteration_count,norm_f);
           if (breakFlag)  printf("iteration break down after %d steps.\n",iteration_count);
           if (maxIterFlag)  printf("maximum number of iteration step %s is reached.\n",maxit);
      }
      if (breakFlag) Status=SOLVER_BREAKDOWN;
      if (maxIterFlag) Status=SOLVER_MAXITER_REACHED;
  }
  TMPMEMFREE(f);
  TMPMEMFREE(step);
printf("STATUS return = %d\n",Status);
  return Status;
}
