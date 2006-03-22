/* $Id$ */

/*
   Crude modifications and translations for Paso by Matt Davies and Lutz Gross
*/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with
*  preconditioning.
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  A       (input) 
*
*  R       (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of inital guess X
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. 
*
*  ITER    (input/output) INT
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) 
*          On output, the final value of this measure.
*
*  return value
*
*          = SOLVER_NO_ERROR: Successful exit. Iterated approximate solution returned.
*          = SOLVER_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVER_BREAKDOWN: If parameters RHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters RHO or OMEGA become smaller
*
*  ==============================================================
*/

err_t Paso_Solver_BiCGStab(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double * tolerance,
    Paso_Performance* pp) {


  /* Local variables */
  double *rtld=NULL,*p=NULL,*v=NULL,*t=NULL,*phat=NULL,*shat=NULL,*s=NULL;
  double beta,norm_of_residual,sum_1,sum_2,sum_3,sum_4,norm_of_residual_global;
  double alpha, omega, omegaNumtr, omegaDenumtr, rho, tol, rho1;
  dim_t num_iter=0,maxit,num_iter_global;
  dim_t i0;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  dim_t status = SOLVER_NO_ERROR;

  /* adapt original routine parameters */
  dim_t n = A->num_cols * A-> col_block_size;;
  double * resid = tolerance;

  /* Executable Statements */

  /*     allocate memory: */
  rtld=TMPMEMALLOC(n,double);
  p=TMPMEMALLOC(n,double);
  v=TMPMEMALLOC(n,double);
  t=TMPMEMALLOC(n,double);
  phat=TMPMEMALLOC(n,double);
  shat=TMPMEMALLOC(n,double);
  s=TMPMEMALLOC(n,double);

  /*     Test the input parameters. */

  if (n < 0) {
    status = SOLVER_INPUT_ERROR;
  } else if (rtld==NULL || p==NULL || v==NULL || t==NULL || phat==NULL || shat==NULL || s==NULL) {
    status = SOLVER_MEMORY_ERROR;
  } else {

    /* now bicgstab starts : */
    maxit = *iter;
    tol = *resid;
   
#pragma omp parallel firstprivate(maxit,tol,convergeFlag,maxIterFlag,breakFlag) \
       private(rho,omega,num_iter,norm_of_residual,beta,alpha,rho1)
    {
      num_iter =0;

      /* initialize arrays */
   
      #pragma omp for private(i0) schedule(static)
      for (i0 = 0; i0 < n; i0++) {
	rtld[i0]=0;
	p[i0]=0;
	v[i0]=0;
	t[i0]=0;
	phat[i0]=0;
	shat[i0]=0;
      }
      #pragma omp for private(i0) schedule(static)
      for (i0 = 0; i0 < n; i0++) rtld[i0] = r[i0];
   
      /*     Perform BiConjugate Gradient Stabilized iteration. */
   
    L10:
      ++(num_iter);
      #pragma omp barrier
      #pragma omp master
      {
	sum_1 = 0;
	sum_2 = 0;
	sum_3 = 0;
	sum_4 = 0;
	omegaNumtr = 0.0;
	omegaDenumtr = 0.0;
      }
      #pragma omp barrier
      #pragma omp for private(i0) reduction(+:sum_1) schedule(static)
      for (i0 = 0; i0 < n; i0++) sum_1 += rtld[i0] * r[i0];
      rho = sum_1;
      
      if (! (breakFlag = (ABS(rho) <= TOLERANCE_FOR_SCALARS))) {
	/*        Compute vector P. */
      
	if (num_iter > 1) {
	  beta = rho / rho1 * (alpha / omega);
          #pragma omp for private(i0) schedule(static)
	  for (i0 = 0; i0 < n; i0++) p[i0] = r[i0] + beta * (p[i0] - omega * v[i0]);
	} else {
          #pragma omp for private(i0) schedule(static)
	  for (i0 = 0; i0 < n; i0++) p[i0] = r[i0];
	}
   
	/*        Compute direction adjusting vector PHAT and scalar ALPHA. */
   
        Paso_Solver_solvePreconditioner(A,&phat[0], &p[0]);
	Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, &phat[0],ZERO, &v[0]);
   
        #pragma omp for private(i0) reduction(+:sum_2) schedule(static)
	for (i0 = 0; i0 < n; i0++) sum_2 += rtld[i0] * v[i0];
        if (! (breakFlag = (ABS(sum_2) <= TOLERANCE_FOR_SCALARS))) {
	   alpha = rho / sum_2;

           #pragma omp for private(i0) reduction(+:sum_3) schedule(static) 
	   for (i0 = 0; i0 < n; i0++) {
	     r[i0] -= alpha * v[i0];
	     s[i0] = r[i0];
	     sum_3 += s[i0] * s[i0];
	   }
	   norm_of_residual = sqrt(sum_3);
        
	   /*        Early check for tolerance. */
	   if ( (convergeFlag = (norm_of_residual <= tol)) ) {
             #pragma omp for  private(i0) schedule(static)
	     for (i0 = 0; i0 < n; i0++) x[i0] += alpha * phat[i0];
	     maxIterFlag = FALSE;
	     breakFlag = FALSE;
	   } else {
	     /*           Compute stabilizer vector SHAT and scalar OMEGA. */
             Paso_Solver_solvePreconditioner(A,&shat[0], &s[0]);
	     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, &shat[0],ZERO,&t[0]);
   
             #pragma omp for private(i0) reduction(+:omegaNumtr,omegaDenumtr) schedule(static)
	     for (i0 = 0; i0 < n; i0++) {
	       omegaNumtr +=t[i0] * s[i0];
	       omegaDenumtr += t[i0] * t[i0];
	     }
             if (! (breakFlag = (ABS(omegaDenumtr) <= TOLERANCE_FOR_SCALARS))) {
	        omega = omegaNumtr / omegaDenumtr;
   
                #pragma omp for private(i0) reduction(+:sum_4) schedule(static)
	        for (i0 = 0; i0 < n; i0++) {
	          x[i0] += alpha * phat[i0] + omega * shat[i0];
	          r[i0] -= omega * t[i0];
	          sum_4 += r[i0] * r[i0];
	        }
	        norm_of_residual = sqrt(sum_4);
	        convergeFlag = norm_of_residual <= tol;
	        maxIterFlag = num_iter == maxit;
	        breakFlag = (ABS(omega) <= TOLERANCE_FOR_SCALARS);
	      }
	   }
	}
	if (!(convergeFlag || maxIterFlag || breakFlag)) {
	  rho1 = rho;
	  goto L10;
	}
      }
      /* end of iteration */
      #pragma omp master 
      {
	num_iter_global=num_iter;
	norm_of_residual_global=norm_of_residual;
	if (maxIterFlag) { 
	    status = SOLVER_MAXITER_REACHED;
	} else if (breakFlag) {
	    status = SOLVER_BREAKDOWN;
	}
      }
    }  /* end of parallel region */
  }
  TMPMEMFREE(rtld);
  TMPMEMFREE(p);
  TMPMEMFREE(v);
  TMPMEMFREE(t);
  TMPMEMFREE(phat);
  TMPMEMFREE(shat);
  TMPMEMFREE(s);
  *iter=num_iter_global;
  *resid=norm_of_residual_global;

  /*     End of BICGSTAB */
  return status;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:40  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:49  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
