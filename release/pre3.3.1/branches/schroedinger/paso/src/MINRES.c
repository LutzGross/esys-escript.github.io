
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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
*  MINRES solves the linear system A*x = b
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
*          = SOLVER_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVER_MEMORY_ERROR :
*          = SOLVER_NEGATIVE_NORM_ERROR :
*
*  ==============================================================
*/

/* #define PASO_DYNAMIC_SCHEDULING_MVM */

#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined __OPENMP 
#define USE_DYNAMIC_SCHEDULING
#endif

err_t Paso_Solver_MINRES(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double *tol,
    double * tolerance,
    Paso_Performance* pp) {

  /* Local variables */
  
  dim_t num_iter=0,maxit;
  bool_t breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  err_t status = SOLVER_NO_ERROR;
  dim_t n = Paso_SystemMatrix_getTotalNumRows(A);
  double  *w=NULL, *w1=NULL, *w2=NULL, *r1=NULL, *r2=NULL, *y=NULL, *v=NULL;

  double Anorm,ynorm,oldb,dbar,epsln,phibar,rhs1,rhs2,rnorm,tnorm2,ynorm2,cs,sn,eps,s,alfa,denom,z,beta1,beta;
  double gmax,gmin,oldeps,delta,gbar,gamma,phi;
 
  double norm_of_residual;
  
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Start of Calculation :                                        */
/*   ---------------------                                         */
/*                                                                 */
/*                                                                 */
  w=TMPMEMALLOC(n,double);
  w1=TMPMEMALLOC(n,double);
  w2=TMPMEMALLOC(n,double);
  r1=TMPMEMALLOC(n,double);
  r2=TMPMEMALLOC(n,double);
  y=TMPMEMALLOC(n,double);
  v=TMPMEMALLOC(n,double);
  
 
 if (w ==NULL || w1== NULL || w2== NULL || r1 == NULL || r2== NULL || y==NULL || v==NULL ) {
     status=SOLVER_MEMORY_ERROR;
  }
 
    maxit = *iter;

 /*     Test the input parameters. */
  if (n < 0 || maxit<=0 ) {
    status=SOLVER_INPUT_ERROR;
  }
  
  Paso_Copy(n,r1,r);
  
  Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
  Paso_Solver_solvePreconditioner(A,y,r1);
  Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
  
  beta1=Paso_InnerProduct(n,r1,y,A->mpi_info);
  if (beta1<0) {
    status=SOLVER_NEGATIVE_NORM_ERROR;
  }
  
  if (beta1>0) {
     beta1=sqrt(beta1);
  }
  
  Performance_startMonitor(pp,PERFORMANCE_SOLVER);
  
  Paso_zeroes(n,w);
  Paso_zeroes(n,w2);
  
  Paso_Copy(n,r2,r1);
  
  Anorm = 0;
  ynorm = 0;
  oldb   = 0;
  beta   = beta1;
  dbar   = 0;
  epsln  = 0;
  phibar = beta1;
  rhs1   = beta1;
  rhs2   = 0;
  rnorm  = phibar;
  tnorm2 = 0;
  ynorm2 = 0;
  cs     = -1;
  sn     = 0;
  eps    = 0.0001;
 
  while (!(convergeFlag || maxIterFlag || breakFlag || (status !=SOLVER_NO_ERROR) ))
  {
          
     s=1/beta;
     Paso_Update(n, 0., v, s, y);
     
     Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
     Performance_startMonitor(pp,PERFORMANCE_MVM);
     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, v,ZERO,y);
     Performance_stopMonitor(pp,PERFORMANCE_MVM);
     Performance_startMonitor(pp,PERFORMANCE_SOLVER);
    
     if (num_iter >= 1) {
        Paso_Update(n, 1., y, -(beta/oldb), r1);
     }

     alfa = Paso_InnerProduct(n,v,y,A->mpi_info);
     Paso_Update(n, 1., y, -(alfa/beta), r2);
     Paso_Copy(n,r1,r2);
     Paso_Copy(n,r2,y);

     Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
     Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
     Paso_Solver_solvePreconditioner(A,y,r2);
     Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
     Performance_startMonitor(pp,PERFORMANCE_SOLVER);

     oldb   = beta;                         
     beta   = Paso_InnerProduct(n,y,r2,A->mpi_info);           
     if (beta<0) {
        status=SOLVER_NEGATIVE_NORM_ERROR;
     }

    beta   = sqrt( beta );
    tnorm2 = tnorm2 + alfa*alfa + oldb*oldb + beta*beta;
        
    if (num_iter==0) {
        gmax   = ABS(alfa);      
        gmin   = gmax;
    }

     /* Apply previous rotation Qk-1 to get
       [deltak epslnk+1] = [cs  sn][dbark    0   ]
       [gbar k dbar k+1]   [sn -cs][alfak betak+1]. */
    
     oldeps = epsln;
     delta  = cs * dbar  +  sn * alfa ; 
     gbar   = sn * dbar  -  cs * alfa ; 
     epsln  =               sn * beta ; 
     dbar   =            -  cs * beta;
     
     gamma  = sqrt(gbar*gbar+beta*beta) ;
     gamma  = MAX(gamma,eps) ;
     cs     = gbar / gamma ;            
     sn     = beta / gamma ;            
     phi    = cs * phibar ;             
     phibar = sn * phibar ;             

     /* Update  x. */

     denom = 1/gamma ;
     Paso_Copy(n,w1,w2);
     Paso_Copy(n,w2,w);
     
     Paso_LinearCombination(n,w,denom,v,-(denom*oldeps),w1);
     Paso_Update(n, 1., w, -(delta*denom), w2) ;
     Paso_Update(n, 1., x, phi, w) ;

     /* Go round again. */

     gmax   = MAX(gmax,gamma);
     gmin   = MIN(gmin,gamma);
     z      = rhs1 / gamma;
     ynorm2 = z*z  +  ynorm2;
     rhs1   = rhs2 -  delta*z;
     rhs2   =      -  epsln*z;

     Anorm  = sqrt( tnorm2 ) ;
     ynorm  = sqrt( ynorm2 ) ;

     rnorm  = phibar;

     maxIterFlag = (num_iter > maxit);
     norm_of_residual=rnorm;
     convergeFlag=(norm_of_residual<Anorm*ynorm*(*tolerance));
    
    
     if (maxIterFlag) {
         status = SOLVER_MAXITER_REACHED;
     } else if (breakFlag) {
         status = SOLVER_BREAKDOWN;
     }
    ++(num_iter);
    /*printf("residual norm %.10f < %.10f %.10f %.10f \n",rnorm,Anorm*ynorm*(*tolerance), Anorm*ynorm, (*tolerance));*/
  }
    /* end of iteration */
    
    Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
    TMPMEMFREE(w);
    TMPMEMFREE(w1);
    TMPMEMFREE(w2); 
    TMPMEMFREE(r1);
    TMPMEMFREE(r2);
    TMPMEMFREE(y); 
    TMPMEMFREE(v); 
  
    *iter=num_iter;
    *tol=norm_of_residual;
    
  /*     End of MINRES */
  return status;
}
