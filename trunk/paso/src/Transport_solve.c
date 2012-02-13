
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
#include "FCT_Solver.h"
#include "ReactiveSolver.h"
#include "Solver.h"
#include "PasoUtil.h"


void Paso_TransportProblem_solve(Paso_TransportProblem* fctp, double* u, double dt, double* u0, double* q, Paso_Options* options) {
   const double reduction_after_divergence_factor = 0.5;
   const dim_t num_failures_max=50;
   
   Paso_Performance pp;
   Paso_ReactiveSolver *rsolver=NULL;
   Paso_FCT_Solver *fctsolver=NULL;
      

   dim_t i_substeps=0, n_substeps=1, num_failures=0;
   double *u_save=NULL, *u2=NULL;
   double  dt2,t=0, dt3;
   err_t errorCode=SOLVER_NO_ERROR;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const dim_t blockSize=Paso_TransportProblem_getBlockSize(fctp);
   options->time_step_backtracking_used=FALSE;
   options->num_iter=0;   

   if (dt<=0.) {
       Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: dt must be positive.");
   }
   if (blockSize>1) {
       Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: block size >1 is not supported.");
    } 
   if (options->verbose) {
      if (options->ode_solver == PASO_BACKWARD_EULER) {
         printf("Paso_TransportProblem_solve: Backward Euler is used (dt = %e)\n",dt);
      } else  if (options->ode_solver == PASO_LINEAR_CRANK_NICOLSON) {
         printf("Paso_TransportProblem_solve: linear Crank-Nicolson is used (dt = %e).\n",dt);
      } else  if (options->ode_solver == PASO_CRANK_NICOLSON) {
         printf("Paso_TransportProblem_solve: Crank-Nicolson is used (dt = %e).\n",dt);
      } else {
	 Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: unknown ODE solver."); 
      }
   }
   if (Esys_noError()) {
         Paso_TransportProblem_getSafeTimeStepSize(fctp); 
        /* 
         *  allocate memory
         *
         */
          fctsolver=Paso_FCT_Solver_alloc(fctp,options);
          rsolver=Paso_ReactiveSolver_alloc(fctp); 
	  u_save=TMPMEMALLOC(n,double);
	  u2=TMPMEMALLOC(n,double);
          Esys_checkPtr(u_save);
	  Esys_checkPtr(u2);
    }
   if (Esys_noError()) {  
       /*
        * let the show begin!!!!
        *
        */
	const double dt_R=fctp->dt_max_R;
        const double dt_T=fctp->dt_max_T;
	dt2=dt;
	if (dt_R < LARGE_POSITIVE_FLOAT) dt2 = MIN(dt_R *2, dt); /* as we half the steps size for the RT bit */
	if (dt_T < LARGE_POSITIVE_FLOAT) {
	      if ( ( options->ode_solver == PASO_LINEAR_CRANK_NICOLSON) || ( options->ode_solver == PASO_CRANK_NICOLSON ) ) {
	         dt2 = MIN(dt_T, dt); 
	      } /* PASO_BACKWARD_EULER does not require a restriction */
	}
	
        num_failures=0;
	Paso_Copy(n,u,u0); /* copy initial value to return */
		
	while( (dt-t)>dt*sqrt(EPSILON) && Esys_noError()) {

	    n_substeps=ceil((dt-t)/dt2);
	    if ( n_substeps <= 0) {
	       Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: time stepping break down.");
	    } else { 
	      dt3=(dt-t)/n_substeps;
	      if (options->verbose) printf("Paso_TransportProblem_solve: number of substeps = %d with dt = %e.\n",n_substeps,dt3);
	      	      /* initialize the iteration matrix */
	      Paso_FCT_Solver_initialize(dt3, fctsolver, options, &pp);
	      Paso_ReactiveSolver_initialize(dt3/2, rsolver, options);
	      errorCode = SOLVER_NO_ERROR;

	      /* start iteration */
	      for (i_substeps =0; (i_substeps<n_substeps) && (errorCode==SOLVER_NO_ERROR) && Esys_noError(); i_substeps++) {
	        if (options->verbose) printf("Paso_TransportProblem_solve: substep %d of %d at t = %e (dt = %e)\n",i_substeps,n_substeps,t+dt3,dt3);
   	        Paso_Copy(n,u_save, u); /* create copy for restart in case of failure */
		/* updates u */
	        errorCode=Paso_ReactiveSolver_solve(rsolver,fctp,u2, u ,q, options, &pp); /* Mu_t=Du+q u(0)=u */ 
	        if (errorCode == SOLVER_NO_ERROR) { 
                  errorCode=Paso_FCT_Solver_update(fctsolver, u, u2,options, &pp); /* Mv_t=Lv   v(0)=u(dt/2) */
	             
	        }
	        if (errorCode == SOLVER_NO_ERROR) {   
                   errorCode=Paso_ReactiveSolver_solve(rsolver,fctp,u2, u, q, options, &pp); /*  Mu_t=Du+q u(dt/2)=v(dt/2) */
	        }
                

                /*
	         * error handling:
	         */
                if (errorCode == SOLVER_NO_ERROR) {
                    num_failures=0;
		    t+=dt3;
		    Paso_Copy(n,u, u2);
		}
	      }
	      if (errorCode == SOLVER_NO_ERROR) {
	      } else if ( (errorCode == SOLVER_MAXITER_REACHED) || (errorCode == SOLVER_DIVERGENCE) ) {
                    /* if num_failures_max failures in a row: give up */
                    if (num_failures >= num_failures_max) {
                       Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: No convergence after time step reductions.");
                    } else {
                       options->time_step_backtracking_used=TRUE;
                       if (options->verbose) printf("Paso_TransportProblem_solve: no convergence. Time step size is reduced.\n");
                       dt2=dt3*reduction_after_divergence_factor;
                       num_failures++;
		       Paso_Copy(n,u, u_save); /* reset initial value */
                    }
	      } else if (errorCode == SOLVER_INPUT_ERROR) {
	        Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: input error for solver.");
	      } else if (errorCode == SOLVER_MEMORY_ERROR) {
	        Esys_setError(MEMORY_ERROR,"Paso_TransportProblem_solve: memory allocation failed.");
	      } else if (errorCode == SOLVER_BREAKDOWN) {
	        Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: solver break down.");
	      } else if (errorCode == SOLVER_NEGATIVE_NORM_ERROR) {
	        Esys_setError(VALUE_ERROR,"Paso_TransportProblem_solve: negative norm.");
	      } else {
	        Esys_setError(SYSTEM_ERROR,"Paso_TransportProblem_solve: general error.");
	      }
	    }
	}
     } /* end of time loop */ 
    /* 
     *  clean-up:
     *
     */

    Paso_FCT_Solver_free(fctsolver);
    Paso_ReactiveSolver_free(rsolver);
    TMPMEMFREE(u_save);
    TMPMEMFREE(u2);
}


double Paso_TransportProblem_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
   double dt_max=0., dt_R=LARGE_POSITIVE_FLOAT, dt_T=LARGE_POSITIVE_FLOAT;
   index_t fail_loc=0;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   
   if ( ! fctp->valid_matrices) {
     /* set row-sum of mass_matrix */
     Paso_SystemMatrix_rowSum(fctp->mass_matrix,fctp->lumped_mass_matrix);
     {  
        /* check for positive entries in lumped_mass_matrix and set negative value for constraints */
	index_t fail=0, i;
        #pragma omp parallel private(i,fail_loc )
        {
	       
               #pragma omp for schedule(static)
               for (i=0;i<n;++i) {
                  const double m_i=fctp->lumped_mass_matrix[i];
		  if ( m_i > 0 ) {
		     if (fctp->constraint_mask[i] > 0 ) fctp->lumped_mass_matrix[i]=-1.;
		  } else {
		      fail_loc=1;
		  }
               }
               #pragma omp critical 
               {
		  fail=MIN(fail, fail_loc);
               }
        }
        #ifdef ESYS_MPI
        { 
	       fail_loc=fail;
               MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MIN, fctp->mpi_info->comm);
	}
        #endif
        if (fail < 0 )
	   Esys_setError(VALUE_ERROR, "Paso_FCTSolver_getSafeTimeStepSize: negative mass matrix entries detected.");
     }
     /* split off row-sum from transport_matrix */
     Paso_SystemMatrix_makeZeroRowSums(fctp->transport_matrix,fctp->reactive_matrix);
     /* get a copy of the main diagonal of the mass matrix */
     Paso_SystemMatrix_copyFromMainDiagonal(fctp->mass_matrix,fctp->main_diagonal_mass_matrix);

     if (Esys_noError()) {
          dt_R=Paso_ReactiveSolver_getSafeTimeStepSize(fctp);
	  dt_T=Paso_FCT_Solver_getSafeTimeStepSize(fctp);
     }
     if (Esys_noError()) {
          fctp->dt_max_R=dt_R;
          fctp->dt_max_T=dt_T;
          fctp->valid_matrices=TRUE;
	  dt_max=MIN( 2 * dt_R,dt_T);
       }
   } else {
        dt_max=MIN(2 * fctp->dt_max_R,fctp->dt_max_T);   /* factor 2 as we use operator splitting */
   }
   return dt_max;
}
