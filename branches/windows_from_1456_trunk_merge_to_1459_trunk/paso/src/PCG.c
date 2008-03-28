
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/* PCG iterations */

#include "SystemMatrix.h"
#include "Paso.h"
#include "Solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PASO_MPI
#include <mpi.h>
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
  dim_t n = Paso_SystemMatrix_getTotalNumRows(A);
  double *resid = tolerance, *rs=NULL, *p=NULL, *v=NULL, *x2=NULL ;
  double tau_old,tau,beta,delta,gamma_1,gamma_2,alpha,sum_1,sum_2,sum_3,sum_4,sum_5,tol;
#ifdef PASO_MPI
  double loc_sum[2], sum[2];
#endif
  double norm_of_residual,norm_of_residual_global;
  register double d;

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
           for (i0=0;i0<n;i0++) sum_1+=v[i0]*r[i0]; /* Limit to local values of v[] and r[] */
           #ifdef PASO_MPI
	        /* In case we have many MPI processes, each of which may have several OMP threads:
	           OMP master participates in an MPI reduction to get global sum_1 */
                #pragma omp master
	        {
	          loc_sum[0] = sum_1;
	          MPI_Allreduce(loc_sum, &sum_1, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
	        }
           #endif
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
	   Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, p,ZERO,v);
           Performance_stopMonitor(pp,PERFORMANCE_MVM);
           Performance_startMonitor(pp,PERFORMANCE_SOLVER);
           /* delta=p*v */
           #pragma omp for private(i0) reduction(+:sum_2) schedule(static)
           for (i0=0;i0<n;i0++) sum_2+=v[i0]*p[i0];
           #ifdef PASO_MPI
               #pragma omp master
	       {
	         loc_sum[0] = sum_2;
	         MPI_Allreduce(loc_sum, &sum_2, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
	       }
           #endif
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
               #ifdef PASO_MPI
                   #pragma omp master
	           {
	             loc_sum[0] = sum_3;
	             loc_sum[1] = sum_4;
	             MPI_Allreduce(loc_sum, sum, 2, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
	             sum_3=sum[0];
	             sum_4=sum[1];
	           }
                #endif
                gamma_1= ( (ABS(sum_3)<= ZERO) ? 0 : -sum_4/sum_3) ;
                gamma_2= ONE-gamma_1;
                #pragma omp for private(i0) schedule(static)
                for (i0=0;i0<n;++i0) {
                  rs[i0]=gamma_2*rs[i0]+gamma_1*r[i0];
                  x2[i0]+=alpha*p[i0];
                  x[i0]=gamma_2*x[i0]+gamma_1*x2[i0];
                }
                #pragma omp for private(i0) reduction(+:sum_5) schedule(static)
                for (i0=0;i0<n;++i0) sum_5+=rs[i0]*rs[i0];
                #ifdef PASO_MPI
                   #pragma omp master
	           {
	              loc_sum[0] = sum_5;
	              MPI_Allreduce(loc_sum, &sum_5, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
	           }
                #endif
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
