
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: Transport solver
 *
 * solves Mu_t=Ku+q
 *        
 *  using an operator splitting approach (K=L+D where L is row sum zero and D is a diagonal matrix):
 *
 *   - Mu_t=Du+q u(0)=u             (reactive part)
 *   - Mv_t=Lv   v(0)=u(dt/2)       (transport part uses flux correction (FCT))
 *   - Mu_t=Du+q u(dt/2)=v(dt/2)    (reactive part)
 *
 *  to return u(dt)
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Transport.h"
#include "FCTSolver.h"
#include "ReactiveSolver.h"
#include "Solver.h"
#include "PasoUtil.h"


void Paso_TransportProblem_solve(Paso_TransportProblem* fctp, double* u, double dt, double* u0, double* q, Paso_Options* options) {
   const double reduction_after_divergence_factor = 0.5;
   const dim_t num_failures_max=50;
   const double increase_after_convergence_factor = 1.1; /* has not much of an impact if substepping is beeing used */
   
   Paso_Performance pp;  
   Paso_Function *fctfunction=NULL;
   Paso_ReactiveSolver *rsolver=NULL;

   dim_t i_substeps, n_substeps, num_failures=0;
   double *u_save=NULL;
   double dt_max=LARGE_POSITIVE_FLOAT, dt2,t=0, dt3;
   err_t errorCode=SOLVER_NO_ERROR;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const dim_t blockSize=Paso_TransportProblem_getBlockSize(fctp);


   if (dt<=0.) {
       Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: dt must be positive.");
   }
   if (blockSize>1) {
       Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: block size >0 is not supported.");
    } 
   if (options->verbose) {
      if (fctp->useBackwardEuler) {
         printf("Paso_TransportProblem_solve: Backward Euler is used (dt = %e)\n",dt);
      } else {
         printf("Paso_TransportProblem_solve: Crank-Nicolson is used (dt = %e).\n",dt);
      }
   }
   if (Paso_noError()) {
        dt_max=Paso_TransportProblem_getSafeTimeStepSize(fctp);
        /* 
         *  allocate memory
         *
         */
          fctfunction=Paso_FCTSolver_Function_alloc(fctp,options);
          rsolver=Paso_ReactiveSolver_alloc(fctp); 
	  u_save=TMPMEMALLOC(n,double);
          Paso_checkPtr(u_save);
    }
   if (Paso_noError()) {  
       /*
        * let the show begin!!!!
        *
        */
	
        /* while(i_substeps<n_substeps && Paso_noError()) */
        if (dt_max < LARGE_POSITIVE_FLOAT) {
            dt2=MIN(dt_max,dt);
        } else {
            dt2=dt;
        }
	num_failures=0;
	
	Paso_Copy(n,u,u0); /* copy initial value */
	Paso_Copy(n,u_save, u); /* create copy for restart in case of failure */
	
	
	while( (dt-t)>dt*sqrt(EPSILON) && Paso_noError()) {

	    n_substeps=ceil((dt-t)/dt2);
	    dt3=(dt-t)/n_substeps;
	    if (options->verbose) printf("Paso_TransportProblem_solve: number of substeps = %d with dt = %e.\n",n_substeps,dt3);
            i_substeps=0;
	    /* initialize the iteration matrix */
            Paso_FCTSolver_Function_initialize(dt3, fctp,options, &pp);
            Paso_ReactiveSolver_initialize(dt3, fctp, options);
	    /* start iteration */
	    for (i_substeps =0; (i_substeps<n_substeps) && (errorCode==SOLVER_NO_ERROR) && Paso_noError(); i_substeps++) {
	        if (options->verbose) printf("Paso_TransportProblem_solve: substep step %d of %d at t = %e (dt = %e)\n",i_substeps,n_substeps,t+dt3,dt3);
   	        /* updates u */
	        errorCode=Paso_ReactiveSolver_solve(rsolver,fctp,u,dt3/2,q, options, &pp); /* Mu_t=Du+q u(0)=u */ 
	        if (errorCode == SOLVER_NO_ERROR) { 
	             errorCode=Paso_FCTSolver_solve(fctfunction, u,dt3,options, &pp); /* Mv_t=Lv   v(0)=u(dt/2) */
	        }
	        if (errorCode == SOLVER_NO_ERROR) {   
                   errorCode=Paso_ReactiveSolver_solve(rsolver,fctp,u,dt3/2,q, options, &pp); /*  Mu_t=Du+q u(dt/2)=v(dt/2) */
	        }
		if (errorCode == SOLVER_NO_ERROR) {
                    t+=dt3;
                    Paso_Copy(n,u_save, u); /* copy u to u_save */
		}
            }
            /*
	     * error handeling:
	     */
             if (errorCode == SOLVER_NO_ERROR) {
                    num_failures=0;
                    i_substeps++;
                    if (fctp->dt_max < LARGE_POSITIVE_FLOAT) {
                       fctp->dt_max=MIN3(dt_max, dt3*increase_after_convergence_factor ,fctp->dt_failed/increase_after_convergence_factor);
                    } else {
                       fctp->dt_max=MIN(dt3*increase_after_convergence_factor ,fctp->dt_failed/increase_after_convergence_factor);
                    }
	     } else if ( (errorCode == SOLVER_MAXITER_REACHED) || (errorCode == SOLVER_DIVERGENCE) ) {
                    /* if num_failures_max failures in a row: give up */
		    fctp->dt_failed=MIN(fctp->dt_failed, dt3);
                    if (num_failures >= num_failures_max) {
                       Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: No convergence after time step reductions.");
                    } else {
                       if (options->verbose) printf("Paso_TransportProblem_solve: no convergence. Time step size is reduced.\n");
                       dt2=dt3*reduction_after_divergence_factor;
                       num_failures++;
		       Paso_Copy(n,u, u_save); /* reset initial value */
		       errorCode = SOLVER_NO_ERROR;
                    }

            } else if (errorCode == SOLVER_INPUT_ERROR) {
	        Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: input error for solver.");
	    } else if (errorCode == SOLVER_MEMORY_ERROR) {
	        Paso_setError(MEMORY_ERROR,"Paso_TransportProblem_solve: memory allication failed.");
	    } else if (errorCode == SOLVER_BREAKDOWN) {
	        Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: solver break down.");
	    } else if (errorCode == SOLVER_NEGATIVE_NORM_ERROR) {
	        Paso_setError(VALUE_ERROR,"Paso_TransportProblem_solve: negative norm.");
	    } else {
	        Paso_setError(SYSTEM_ERROR,"Paso_TransportProblem_solve: general error.");
	    }
	}
     }
    /* 
     *  clean-up:
     *
     */
    Paso_FCTSolver_Function_free(fctfunction);
    Paso_ReactiveSolver_free(rsolver);
    TMPMEMFREE(u_save);
   
}


double Paso_TransportProblem_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
   double dt_max, dt1=LARGE_POSITIVE_FLOAT, dt2=LARGE_POSITIVE_FLOAT;
   if ( ! fctp->valid_matrices) {
     fctp->dt_failed=LARGE_POSITIVE_FLOAT;
     /* set row-sum of mass_matrix */
     Paso_SystemMatrix_rowSum(fctp->mass_matrix,fctp->lumped_mass_matrix);
     /* split off row-sum from transport_matrix */
     Paso_SystemMatrix_makeZeroRowSums(fctp->transport_matrix,fctp->reactive_matrix);
     /* get a copy of the main diagonal of the mass matrix */
     Paso_SystemMatrix_copyFromMainDiagonal(fctp->mass_matrix,fctp->main_diagonal_mass_matrix);

     if (Paso_noError()) {
          dt1=Paso_ReactiveSolver_getSafeTimeStepSize(fctp); 
          if (dt1<LARGE_POSITIVE_FLOAT) dt1*=2;
     }
     if (Paso_noError()) dt2=Paso_FCTSolver_getSafeTimeStepSize(fctp);
     printf("Paso_TransportProblem_getSafeTimeStepSize: dt_max from reactive part = %e\n",dt1);
     printf("Paso_TransportProblem_getSafeTimeStepSize: dt_max from transport part = %e\n",dt2);
     dt_max=MIN(dt1,dt2);

     if ( (dt_max <= 0.) &&  Paso_noError() ) {
            Paso_setError(TYPE_ERROR,"Paso_TransportProblem_solve: dt must be positive.");
      } 
     fctp->dt_max=dt_max;
     fctp->valid_matrices=Paso_noError();
   }
   return fctp->dt_max;
}
