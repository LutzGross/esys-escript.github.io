
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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


/* MINRES iterations */

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
*  MINRES solves the linear system A*x = b
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  R      (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of initial guess x
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
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



err_t Paso_Solver_MINRES(
    Paso_SystemMatrix * A,
    double * R,
    double * X,
    dim_t *iter,
    double * tolerance,
    Paso_Performance* pp) 
{

   double    delta,gamma=0.,gamma_old=0.,eta=0.,dp0=0., c=1.0,c_old=1.0,c_ancient=1.,s=0.0,s_old=0.0,s_ancient, norm_of_residual=0., rnorm_prec=1;
   double tol=1., norm_scal=1.;
    const dim_t maxit = *iter;
    double    alpha_0,alpha_1,alpha_2,alpha_3,dp = 0.0;
    dim_t     num_iter = 0;
    double    *Z=NULL, *W=NULL, *AZ=NULL, *R_old=NULL, *R_ancient=NULL, *W_old=NULL, *W_ancient=NULL, *ZNEW=NULL;
    const dim_t n = Paso_SystemMatrix_getTotalNumRows(A);
    bool_t convergeFlag=FALSE;
    err_t status = SOLVER_NO_ERROR;
/*                                                                 
 *                                                                 
 *   Start of Calculation :                                        
 *   ---------------------                                         
 *                                                                 
 *                                                                 */
   

   /*     Test the input parameters. */
   if (n < 0 || maxit<=0 ) {
      status=SOLVER_INPUT_ERROR;
   }
   
   ZNEW       = new double[n];
   Z       = new double[n];
   AZ    = new double[n];
   W       = new double[n];
   R_old    = new double[n];
   W_old    = new double[n];
   R_ancient   = new double[n];
   W_ancient   = new double[n];
   
   if (R_ancient==NULL || Z==NULL || W==NULL || AZ==NULL || R_old==NULL || W_old==NULL || W_ancient==NULL || ZNEW==NULL) {
      status=SOLVER_MEMORY_ERROR;
   }
      
   if (status ==SOLVER_NO_ERROR) { 
      
      Paso_SystemMatrix_solvePreconditioner(A, Z, R); /*     z  <- Prec*r       */
      /* gamma <- r'*z */
          dp=Paso_InnerProduct(n, R ,Z,A->mpi_info); /* gamma <- r'*z */
	  dp0=dp;
      if (dp<0) {
	 status=SOLVER_NEGATIVE_NORM_ERROR;
      } else if (! ABS(dp)>0) {
	    convergeFlag = TRUE;            /* happy break down */
      } else {
            gamma   = sqrt(dp); /*  gamma <- sqrt(r'*z)  */
            eta  = gamma;
            rnorm_prec = gamma;
            norm_of_residual=Paso_l2(n, R, A->mpi_info);
            norm_scal=rnorm_prec/norm_of_residual;
            tol=(*tolerance)*norm_scal;
      }
   }
   while (!(convergeFlag || (status !=SOLVER_NO_ERROR) ))
   {
        /*    z <- z / gamma     */
           Paso_Scale(n, Z,1./gamma);        

        /*      Az <- A*z   */
           Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, Z,PASO_ZERO,AZ); 

	/*  delta <- Az'.z */
	    delta=Paso_InnerProduct(n,AZ,Z,A->mpi_info); 

       /*  r_new <- Az-delta/gamma * r - gamma/gamma_old r_old */
          if (num_iter>0) Paso_Copy(n, R_ancient, R_old);   /*  r__ancient <- r_old */
          Paso_Copy(n, R_old, R);       /*  r_old <- r */
          Paso_Copy(n, R, AZ);       /*  r <- Az */
	  Paso_AXPY(n, R, -delta/gamma, R_old);     /*  r <- r - delta/gamma v     */
	  if (num_iter>0) Paso_AXPY(n, R, -gamma/gamma_old, R_ancient);   /*  r <- r - gamma/gamma_old r__ancient  */

	/*  z <- prec*r   */
	  Paso_SystemMatrix_solvePreconditioner(A, ZNEW, R); 
        
	
	 dp=Paso_InnerProduct(n,R,ZNEW,A->mpi_info);
	 if (dp < 0.) {
		  status=SOLVER_NEGATIVE_NORM_ERROR;
	 } else if (ABS(dp) == 0.) {
		  convergeFlag = TRUE;            /* happy break down */
	 } else if (ABS(dp) > 0.e-13 * ABS(dp0) ) {
	       /*  gamma <- sqrt(r'*z)   */
		gamma_old=gamma;
	        gamma = sqrt(dp);                             
	       /*    QR factorisation    */
	 
	       c_ancient = c_old; c_old = c; 
               s_ancient = s_old; s_old = s;
	 
	       alpha_0 = c_old * delta - c_ancient * s_old * gamma_old;
	       alpha_1 = sqrt(alpha_0*alpha_0 + gamma*gamma);
	       alpha_2 = s_old * delta + c_ancient * c_old * gamma_old;
	       alpha_3 = s_ancient * gamma_old;
	 
	       /*     Givens rotation    */
	 
	       c = alpha_0 / alpha_1;
	       s = gamma / alpha_1;

               rnorm_prec = rnorm_prec * s;

               /* w_new <- (z-alpha_3 w - alpha_2 w_old)/alpha_1 */
	 
	             if (num_iter>1) Paso_Copy(n, W_ancient, W_old);     /*  w__ancient <- w_old      */
	             if (num_iter>0) Paso_Copy(n, W_old, W);         /*  w_old  <- w          */
	 
	             Paso_Copy(n, W, Z);
	             if (num_iter>1) Paso_AXPY(n, W,- alpha_3,W_ancient); /*  w <- w - alpha_3 w__ancient */
	             if (num_iter>0) Paso_AXPY(n, W,- alpha_2,W_old);  /*  w <- w - alpha_2 w_old  */
   	             Paso_Scale(n, W, 1.0 / alpha_1);      /*  w <- w / alpha_1        */
               /*                                                        */
	       Paso_AXPY(n, X,c * eta,W);      /*  x <- x + c eta w     */ 
	       eta = - s * eta;
	       convergeFlag = rnorm_prec <= tol;
	 } else {
		   status=SOLVER_BREAKDOWN;
	 }
         Paso_Copy(n, Z, ZNEW);       
	 ++(num_iter);
	 if ( !convergeFlag && (num_iter>=maxit)) status = SOLVER_MAXITER_REACHED;
   }
   delete[] Z;
   delete[] ZNEW;
   delete[] AZ;
   delete[] R_old;
   delete[] R_ancient;
   delete[] W;
   delete[] W_old;
   delete[] W_ancient;
      
   *iter=num_iter;
   *tolerance=rnorm_prec/norm_scal;
      
   /*     End of MINRES */
   return status;
}

