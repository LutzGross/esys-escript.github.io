
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/* TFQMR iterations */

#include "SystemMatrix.h"
#include "Paso.h"
#include "Solver.h"
#include "PasoUtil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ESYS_MPI
#include <mpi.h>
#endif

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

/* #define PASO_DYNAMIC_SCHEDULING_MVM */

#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined __OPENMP 
#define USE_DYNAMIC_SCHEDULING
#endif

err_t Paso_Solver_TFQMR(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double * tolerance,
    Paso_Performance* pp) {

  /* Local variables */
  
  int m=1;  
  int j=0;
  
  dim_t num_iter=0,maxit;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  err_t status = SOLVER_NO_ERROR;
  dim_t n = Paso_SystemMatrix_getTotalNumRows(A);
  double  *u1=NULL, *u2=NULL, *y1=NULL, *y2=NULL, *d=NULL, *w=NULL, *v=NULL, *temp_vector=NULL,*res=NULL;

  double eta,theta,tau,rho,beta,alpha,sigma,rhon,c;

  double norm_of_residual;
  
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Start of Calculation :                                        */
/*   ---------------------                                         */
/*                                                                 */
/*                                                                 */
  u1=TMPMEMALLOC(n,double);
  u2=TMPMEMALLOC(n,double);
  y1=TMPMEMALLOC(n,double);
  y2=TMPMEMALLOC(n,double);
  d=TMPMEMALLOC(n,double);
  w=TMPMEMALLOC(n,double);
  v=TMPMEMALLOC(n,double);
  temp_vector=TMPMEMALLOC(n,double);
  res=TMPMEMALLOC(n,double);
 
 if (u1 ==NULL || u2== NULL || y1 == NULL || y2== NULL || d==NULL || w==NULL || v==NULL ) {
     status=SOLVER_MEMORY_ERROR;
  }
 
    maxit = *iter;

 /*     Test the input parameters. */
  if (n < 0 || maxit<=0 ) {
    status=SOLVER_INPUT_ERROR;
  }
  
  Paso_zeroes(n,x);
  
  Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Paso_SystemMatrix_solvePreconditioner(A,res,r);
  Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
  
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Paso_zeroes(n,u2);
  Paso_zeroes(n,y2);
  
  Paso_Copy(n,w,res);
  Paso_Copy(n,y1,res);
      
  Paso_zeroes(n,d);
  
  Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
  Performance_startMonitor(pp,PERFORMANCE_MVM);
  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y1,PASO_ZERO,temp_vector);
  Performance_stopMonitor(pp,PERFORMANCE_MVM);
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
  Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Paso_SystemMatrix_solvePreconditioner(A,v,temp_vector);
  Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  /* v = P^{-1} * A y1 */
  
  Paso_Copy(n,u1,v);
    
  theta = 0.0;
  eta = 0.0;
  
  tau = Paso_l2(n,res,A->mpi_info);
  
  rho = tau * tau;
      
  norm_of_residual=tau;
  
  while (!(convergeFlag || maxIterFlag || breakFlag || (status !=SOLVER_NO_ERROR) ))
  {
          
 
     sigma=Paso_InnerProduct(n,res,v,A->mpi_info);
     
     if (! (breakFlag = (ABS(sigma) == 0.))) {
     
     alpha = rho / sigma;
     
     for (j=0; j<=1; j=j+1)
       {
         /*  Compute y2 and u2 only if you have to */
         if ( j == 1 ){
          Paso_LinearCombination(n,y2,PASO_ONE,y1,-alpha,v); /* y2 = y1 - alpha*v */
          
          Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
          Performance_startMonitor(pp,PERFORMANCE_MVM);
          Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y2,PASO_ZERO,temp_vector);
          Performance_stopMonitor(pp,PERFORMANCE_MVM);
          Performance_startMonitor(pp,PERFORMANCE_SOLVER);
          
          Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
          Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
	  Paso_SystemMatrix_solvePreconditioner(A,u2,temp_vector);  
          Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
          Performance_startMonitor(pp,PERFORMANCE_SOLVER);
          /* u2 = P^{-1} * A y2 */
         } 
         m = 2 * (num_iter+1) - 2 + (j+1);

          if (j==0) { 
            Paso_Update(n,1.,w,-alpha,u1); /* w = w - alpha * u1 */
            Paso_Update(n,( theta * theta * eta / alpha ),d,1.,y1); /* d = ( theta * theta * eta / alpha )*d + y1 */
          }
          if (j==1) {
            Paso_Update(n,1.,w,-alpha,u2);  /* w = w - -alpha * u2 */
            Paso_Update(n,( theta * theta * eta / alpha ),d,1.,y2); /* d = ( theta * theta * eta / alpha )*d + y2 */
          } 
         
         theta =Paso_l2(n,w,A->mpi_info)/tau;
     printf("tau = %e, %e %e\n",tau, Paso_l2(n,w,A->mpi_info)/tau, theta);
         c = PASO_ONE / sqrt ( PASO_ONE + theta * theta );
         tau = tau * theta * c;
         eta = c * c * alpha;
         Paso_Update(n,1.,x,eta,d);       
       }

     breakFlag = (ABS(rho) == 0);

     rhon = Paso_InnerProduct(n,res,w,A->mpi_info);
     beta = rhon / rho;
     rho = rhon;
  
     Paso_LinearCombination(n,y1, PASO_ONE,w,beta,y2); /* y1 = w + beta * y2 */

     Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
     Performance_startMonitor(pp,PERFORMANCE_MVM);
     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y1,PASO_ZERO,temp_vector);
     Performance_stopMonitor(pp,PERFORMANCE_MVM);
     
     Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
     Paso_SystemMatrix_solvePreconditioner(A,u1,temp_vector);
     Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
     Performance_startMonitor(pp,PERFORMANCE_SOLVER);
     /*  u1 = P^{-1} * A y1 */
     
     Paso_LinearCombination(n,temp_vector,PASO_ONE,u2,beta,v); /* t = u2 + beta * v */
     Paso_LinearCombination(n,v,PASO_ONE,u1,beta,temp_vector); /* v = u1 + beta * t */
     }
     maxIterFlag = (num_iter > maxit);
     norm_of_residual=tau*sqrt ( (double) (m + 1 ) );
     convergeFlag=(norm_of_residual<(*tolerance));
    
    
     if (maxIterFlag) {
         status = SOLVER_MAXITER_REACHED;
     } else if (breakFlag) {
         status = SOLVER_BREAKDOWN;
     }
    ++(num_iter);
  }
    /* end of iteration */
    
    Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
    TMPMEMFREE(u1);
    TMPMEMFREE(u2);
    TMPMEMFREE(y1); 
    TMPMEMFREE(y2);
    TMPMEMFREE(d);
    TMPMEMFREE(w); 
    TMPMEMFREE(v); 
    TMPMEMFREE(temp_vector);
    TMPMEMFREE(res);
    *iter=num_iter;
    *tolerance=norm_of_residual;
    
  /*     End of TFQMR */
  return status;
}
