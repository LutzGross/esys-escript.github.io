
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/* TFQMR iterations */

#include "SystemMatrix.h"
#include "Paso.h"
#include "Solver.h"
#include "PasoUtil.h"

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
*          = SOLVEr_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVEr_BREAKDOWN: If parameters rHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters rHO or OMEGA become smaller
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
  double  *u1=NULL, *u2=NULL, *y1=NULL, *y2=NULL, *d=NULL, *w=NULL, *v=NULL, *v_old=NULL;

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
  v_old=TMPMEMALLOC(n,double);
  
 
 if (u1 ==NULL || u2== NULL || y1 == NULL || y2== NULL || d==NULL || w==NULL || v==NULL || v_old==NULL) {
     status=SOLVER_MEMORY_ERROR;
  }
 
    maxit = *iter;

 /*     Test the input parameters. */
  if (n < 0 || maxit<=0 ) {
    status=SOLVER_INPUT_ERROR;
  }
  
  Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Paso_Solver_solvePreconditioner(A,r,r);
  Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
  
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Paso_zeroes(n,u2);
  Paso_zeroes(n,y2);
  
  Paso_Copy(n,w,r);
  Paso_Copy(n,y1,r);
      
  Paso_zeroes(n,d);
  
  Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
  Performance_startMonitor(pp,PERFORMANCE_MVM);
  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y1,PASO_ZERO,v);
  Performance_stopMonitor(pp,PERFORMANCE_MVM);
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
  Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Paso_Solver_solvePreconditioner(A,v,v);
  Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Paso_Copy(n,u1,v);
    
  theta = 0.0;
  eta = 0.0;
  
  tau = Paso_l2(n,r,A->mpi_info);
  
  rho = tau * tau;
      
  norm_of_residual=tau*sqrt ( m + 1 );
  
  while (!(convergeFlag || maxIterFlag || breakFlag || (status !=SOLVER_NO_ERROR) ))
  {
          
 
     sigma=Paso_InnerProduct(n,r,v,A->mpi_info);
     
     if (! (breakFlag = (ABS(sigma) == 0.))) {
     
     alpha = rho / sigma;
     
     for (j=0; j<=1; j=j+1)
       {
         /*  Compute y2 and u2 only if you have to */
         if ( j == 1 ){
          Paso_LinearCombination(n,y2,1.,y1,-alpha,v);
          
          Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
          Performance_startMonitor(pp,PERFORMANCE_MVM);
          Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y2,PASO_ZERO,u2);
          Performance_stopMonitor(pp,PERFORMANCE_MVM);
          Performance_startMonitor(pp,PERFORMANCE_SOLVER);
          
          Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
          Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
          Paso_Solver_solvePreconditioner(A,u2,u2);
          Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
          Performance_startMonitor(pp,PERFORMANCE_SOLVER);
         } 
         m = 2 * (num_iter+1) - 2 + (j+1);
          if (j==0) { 
            Paso_Update(n,1.,w,-alpha,u1);
            Paso_Update(n,( theta * theta * eta / alpha ),d,1.,y1);
          }
          if (j==1) {
            Paso_Update(n,1.,w,-alpha,u2);
            Paso_Update(n,( theta * theta * eta / alpha ),d,1.,y2);
          } 
         
         theta =Paso_l2(n,w,A->mpi_info)/tau;
         c = 1.0 / sqrt ( 1.0 + theta * theta );
         tau = tau * theta * c;
         eta = c * c * alpha;
         Paso_Update(n,1.,x,eta,d);       
       }

     breakFlag = (ABS(rho) == 0);

     rhon = Paso_InnerProduct(n,r,w,A->mpi_info);
     beta = rhon / rho;
     rho = rhon;
  
     Paso_LinearCombination(n,y1,1.,w,beta,y2);

     Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
     
     Performance_startMonitor(pp,PERFORMANCE_MVM);
     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, y1,PASO_ZERO,u1);
     Performance_stopMonitor(pp,PERFORMANCE_MVM);
     
     Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
     Paso_Solver_solvePreconditioner(A,u1,u1);
     Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
     
     Performance_startMonitor(pp,PERFORMANCE_SOLVER);
     
     Paso_Copy(n,v_old,v);
          
     Paso_Update(n,beta,v_old,1,u2);
     Paso_LinearCombination(n,v,1.,u1,beta,v_old);
     }
     maxIterFlag = (num_iter > maxit);
     norm_of_residual=tau*sqrt ( m + 1 );
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
    TMPMEMFREE(v_old);
  
    *iter=num_iter;
    *tolerance=norm_of_residual;
    
  /*     End of TFQMR */
  return status;
}
