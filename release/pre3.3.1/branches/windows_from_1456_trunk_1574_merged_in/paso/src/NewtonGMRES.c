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
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
err_t Paso_Solver_NewtonGMRES(
    Paso_Function *F,
    double *x,
    Paso_Options* options,
    Paso_Performance* pp) 

{
   const double gamma=.9;
   dim_t gmres_iter;
   double stop_tol, fnrm,fnrmo, rat,etaold,etanew, gmres_tol;
   bool_t convergeFlag=FALSE, maxIterFlag=FALSE, breakFlag=FALSE;
   double *f=NULL, *step=NULL;
   err_t Status=SOLVER_NO_ERROR;
   bool_t debug=options->verbose;
   dim_t n = n=F->local_n; 
   dim_t itc=0;
   double atol=options->absolute_tolerance;
   double rtol=options->tolerance;
   dim_t maxit=options->iter_max;
   dim_t lmaxit=options->inner_iter_max;
   bool_t adapt_eta=options->adapt_inner_tolerance;
   double etamax=options->inner_tolerance;
   double eta=etamax;
  /*
   * etamax = Maximum error tolerance for residual in inner
   *              iteration. The inner iteration terminates
   *              when the relative linear residual is
   *              smaller than eta*| F(x_c) |. eta is determined
   *              by the modified Eisenstat-Walker formula if etamax > 0.
   *              If etamax < 0, then eta = |etamax| for the entire
   *              iteration.
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
     fnrm=Paso_l2(n,f,F->mpi_info);
     /*
      * stoping criterium:
      */
     stop_tol=atol + rtol*fnrm;
     itc=0;
     if (debug) printf("Start Jacobi-free Newton scheme:\n\ttolerance = %e\n\tstoping tolerance = %e\n",rtol,stop_tol);
     /* 
      *  main iteration loop
      */
     while (! (convergeFlag || maxIterFlag || breakFlag)) {
         /* 
          * keep track of the ratio (rat = fnrm/frnmo) of successive residual norms and 
          * the iteration counter (itc)
          */
         itc++;
         if (debug) printf("iteration step %d: norm of F =%d\n",itc,fnrm);
         /*
          * call GMRES to get increment
          */
         gmres_iter=lmaxit;
         gmres_tol=eta;
         if (debug) printf("GMRES called with tolerance = %d\n",eta);
         Status=Paso_Solver_NLGMRES(F,f,x,step,&gmres_iter,&gmres_tol,pp);
         itc+=gmres_iter;
         if ((Status==SOLVER_NO_ERROR) || (Status==SOLVER_MAXITER_REACHED)) {
            Status=SOLVER_NO_ERROR;
            /* 
             * update x:
             */
            Paso_Update(n,1.,x,1.,step);
            Paso_FunctionCall(F,f,x);
            fnrm=Paso_l2(n,f,F->mpi_info);
            fnrmo=fnrm; 
            rat=fnrm/fnrmo;
            /*
             *   adjust eta 
             */
             if (adapt_eta) {
                 etaold=eta;
                 etanew=gamma*rat*rat;
                 if (gamma*etaold*etaold >1.) etanew=MAX(etanew,gamma*etaold*etaold);
                 eta=MAX(MIN(etanew,etamax), .5*stop_tol/fnrm);
                }
         } else {
           breakFlag=TRUE;
         }
         maxIterFlag = (itc > maxit);
         convergeFlag = (fnrm <= stop_tol);
      }
  }
  TMPMEMFREE(f);
  TMPMEMFREE(step);
  return Status;
}
