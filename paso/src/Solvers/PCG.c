/* $Id$ */


/* PCG iterations */

#include "../SystemMatrix.h"
#include "../Paso.h"
#include "Solver.h"
/* #include <math.h> */
#ifdef _OPENMP
#include <omp.h>
#endif

/*
*
*  Purpose
*  =======
*
*  PCG solves the linear system A*x = b using the
*  preconditioned conjugate gradient method plus a smoother
*  A has to be symmetric.
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  r       (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of inital guess x
*
*  x       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess.
*
*  ITER    (input/output) INT
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  INFO    (output) INT
*
*          = SOLVER_NO_ERROR: Successful exit. Iterated approximate solution returned.
*          = SOLVEr_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVEr_BREAKDOWN: If parameters rHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters rHO or OMEGA become smaller
*
*  ==============================================================
*/

err_t Paso_Solver_PCG(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double * tolerance,
    Paso_Performance* pp) {


  /* Local variables */
  dim_t num_iter=0,maxit,num_iter_global;
  dim_t i0;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  err_t status = SOLVER_NO_ERROR;
  dim_t n = A->num_cols * A-> col_block_size;
  double *resid = tolerance, *rs=NULL, *p=NULL, *v=NULL, *x2=NULL ;
  double tau_old,tau,beta,delta,gamma_1,gamma_2,alpha,sum_1,sum_2,sum_3,sum_4,sum_5,tol;
  double norm_of_residual,norm_of_residual_global;
  register double r_tmp,d,rs_tmp,x2_tmp,x_tmp;

/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Start of Calculation :                                        */
/*   ---------------------                                         */
/*                                                                 */
/*                                                                 */
  rs=TMPMEMALLOC(n,double);
  p=TMPMEMALLOC(n,double);
  v=TMPMEMALLOC(n,double);
  x2=TMPMEMALLOC(n,double);

  /*     Test the input parameters. */

  if (n < 0) {
    status = SOLVER_INPUT_ERROR;
  } else if (rs==NULL || p==NULL || v==NULL || x2==NULL) {
    status = SOLVER_MEMORY_ERROR;
  } else {
    maxit = *iter;
    tol = *resid;
    #pragma omp parallel firstprivate(maxit,tol,convergeFlag,maxIterFlag,breakFlag) \
                                           private(tau_old,tau,beta,delta,gamma_1,gamma_2,alpha,norm_of_residual,num_iter)
    {
       Performance_startMonitor(pp,PERFORMANCE_SOLVER);
       /* initialize data */
       #pragma omp for private(i0) schedule(static)
       for (i0=0;i0<n;i0++) {
          rs[i0]=r[i0];
          x2[i0]=x[i0];
       } 
       #pragma omp for private(i0) schedule(static)
       for (i0=0;i0<n;i0++) {
          p[i0]=0;
          v[i0]=0;
       } 
       num_iter=0;
       /* start of iteration */
       while (!(convergeFlag || maxIterFlag || breakFlag)) {
           ++(num_iter);
           #pragma omp barrier
           #pragma omp master
           {
	       sum_1 = 0;
	       sum_2 = 0;
	       sum_3 = 0;
	       sum_4 = 0;
	       sum_5 = 0;
           }
           /* v=prec(r)  */
           Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
           Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
           Paso_Solver_solvePreconditioner(A,v,r);
           Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
           Performance_startMonitor(pp,PERFORMANCE_SOLVER);
           /* tau=v*r    */
           #pragma omp for private(i0) reduction(+:sum_1) schedule(static)
           for (i0=0;i0<n;i0++) sum_1+=v[i0]*r[i0];
           tau_old=tau;
           tau=sum_1;
           /* p=v+beta*p */
           if (num_iter==1) {
               #pragma omp for private(i0)  schedule(static)
               for (i0=0;i0<n;i0++) p[i0]=v[i0];
           } else {
               beta=tau/tau_old;
               #pragma omp for private(i0)  schedule(static)
               for (i0=0;i0<n;i0++) p[i0]=v[i0]+beta*p[i0];
           }
           /* v=A*p */
           Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
           Performance_startMonitor(pp,PERFORMANCE_MVM);
	   Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, p,ZERO,v);
           Performance_stopMonitor(pp,PERFORMANCE_MVM);
           Performance_startMonitor(pp,PERFORMANCE_SOLVER);
           /* delta=p*v */
           #pragma omp for private(i0) reduction(+:sum_2) schedule(static)
           for (i0=0;i0<n;i0++) sum_2+=v[i0]*p[i0];
           delta=sum_2;

   
           if (! (breakFlag = (ABS(delta) <= TOLERANCE_FOR_SCALARS))) {
               alpha=tau/delta;
               /* smoother */
               #pragma omp for private(i0) schedule(static)
               for (i0=0;i0<n;i0++) r[i0]-=alpha*v[i0];
               #pragma omp for private(i0,d) reduction(+:sum_3,sum_4) schedule(static)
               for (i0=0;i0<n;i0++) {
                     d=r[i0]-rs[i0];
                     sum_3+=d*d;
                     sum_4+=d*rs[i0];
                }
                gamma_1= ( (ABS(sum_3)<= ZERO) ? 0 : -sum_4/sum_3) ;
                gamma_2= ONE-gamma_1;
                #pragma omp for private(i0,x2_tmp,x_tmp,rs_tmp) schedule(static)
                for (i0=0;i0<n;++i0) {
                  rs[i0]=gamma_2*rs[i0]+gamma_1*r[i0];
                  x2[i0]+=alpha*p[i0];
                  x[i0]=gamma_2*x[i0]+gamma_1*x2[i0];
                }
                #pragma omp for private(i0) reduction(+:sum_5) schedule(static)
                for (i0=0;i0<n;++i0) sum_5+=rs[i0]*rs[i0];
                norm_of_residual=sqrt(sum_5);
                convergeFlag = norm_of_residual <= tol;
                maxIterFlag = num_iter == maxit;
                breakFlag = (ABS(tau) <= TOLERANCE_FOR_SCALARS);
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
       Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
    }  /* end of parallel region */
    TMPMEMFREE(rs);
    TMPMEMFREE(x2);
    TMPMEMFREE(v);
    TMPMEMFREE(p);
    *iter=num_iter_global;
    *resid=norm_of_residual_global;
  }
  /*     End of PCG */
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
