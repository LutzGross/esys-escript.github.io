
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

/* Paso: Transport solver with flux correction (L is row sum zero) 
 *
 *   - Mv_t=Lv   v(0)=u
 *
 *  to return v(dt)
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "FCT_Solver.h"
#include "Preconditioner.h"
#include "PasoUtil.h"


Paso_FCT_Solver* Paso_FCT_Solver_alloc(Paso_TransportProblem *fctp, Paso_Options* options)
{
    Paso_FCT_Solver* out=NULL;
    const dim_t blockSize=Paso_TransportProblem_getBlockSize(fctp);
    const dim_t n =  Paso_TransportProblem_getTotalNumRows(fctp);
   
    out=MEMALLOC(1,Paso_FCT_Solver);
    if (! Esys_checkPtr(out)) {
      	out->transportproblem  = Paso_TransportProblem_getReference(fctp);
        out->mpi_info          = Esys_MPIInfo_getReference(fctp->mpi_info);
	out->flux_limiter      = Paso_FCT_FluxLimiter_alloc(fctp);
	out->b                 = MEMALLOC(n, double);
	if ( (options->ode_solver == PASO_CRANK_NICOLSON) || (options->ode_solver == PASO_BACKWARD_EULER) ) {
	    out->du = MEMALLOC(n, double);
	    out->z = MEMALLOC(n, double);
        } else {
	     out->du = NULL;
	     out->z=NULL;
        }
	out->u_coupler = Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp), blockSize);
        out->u_old_coupler = Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp), blockSize);
	out->omega=0;
 	
	if ( options->ode_solver == PASO_LINEAR_CRANK_NICOLSON ) {
	     out->method   = PASO_LINEAR_CRANK_NICOLSON;
	} else if ( options->ode_solver == PASO_CRANK_NICOLSON ) {
	     out->method = PASO_CRANK_NICOLSON;
	} else if ( options->ode_solver == PASO_BACKWARD_EULER ) {
	     out->method = PASO_BACKWARD_EULER;
	} else {
	    Esys_setError(VALUE_ERROR, "Paso_FCT_Solver_alloc: unknown integration scheme."); 
	    out->method = UNKNOWN;
	}
	  
    }
    
    if (Esys_noError()) {
        return out;
    } else {
        Paso_FCT_Solver_free(out);
        return NULL;
    }    
  
}

void Paso_FCT_Solver_free(Paso_FCT_Solver *in)
{
    if (in != NULL) {
          Paso_TransportProblem_free(in->transportproblem);
	  Paso_FCT_FluxLimiter_free(in->flux_limiter);
          Esys_MPIInfo_free(in->mpi_info);
	  Paso_Coupler_free(in->u_old_coupler);
	  Paso_Coupler_free(in->u_coupler);

	  MEMFREE(in->b);
	  MEMFREE(in->z);
	  MEMFREE(in->du);
          MEMFREE(in);
      
    }
}

double Paso_FCT_Solver_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
   dim_t i, n;
   double dt_max=LARGE_POSITIVE_FLOAT;
   index_t fail=0;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   /* set low order transport operator */
   Paso_FCT_setLowOrderOperator(fctp);
          
   if (Esys_noError()) {
        /*
         *  calculate time step size:                                           
        */
        dt_max=LARGE_POSITIVE_FLOAT;
	fail=0;
        #pragma omp parallel private(i)
        {
               double dt_max_loc=LARGE_POSITIVE_FLOAT;
	       index_t fail_loc=0;
               #pragma omp for schedule(static)
               for (i=0;i<n;++i) {
                  const double l_ii=fctp->main_diagonal_low_order_transport_matrix[i];
                  const double m_i=fctp->lumped_mass_matrix[i];
		  if ( (m_i > 0) ) {
		      if (l_ii<0) dt_max_loc=MIN(dt_max_loc,m_i/(-l_ii));
		  } else {
		      fail_loc=-1;
		  }
               }
               #pragma omp critical 
               {
                  dt_max=MIN(dt_max,dt_max_loc);
		  fail=MIN(fail, fail_loc);
               }
        }
        #ifdef ESYS_MPI
        {
	       double rtmp_loc[2], rtmp[2];
               rtmp_loc[0]=dt_max;
	       rtmp_loc[1]= (double) fail;
               MPI_Allreduce(rtmp_loc, rtmp, 2, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
	       dt_max=rtmp[0];
	       fail = rtmp[1] < 0 ? -1 : 0;
	}
        #endif
        if (fail < 0 ) {
	   Esys_setError(VALUE_ERROR, "Paso_FCTSolver_getSafeTimeStepSize: negative mass matrix entries detected.");
	   return -1;
	} else {
	    if (dt_max<LARGE_POSITIVE_FLOAT) dt_max*=2.;
	}
   }
   return dt_max;
}

/* modifies the main diagonal of the iteration matrix to introduce new dt */
void Paso_FCT_Solver_initialize(const double dt, Paso_FCT_Solver *fct_solver, Paso_Options* options, Paso_Performance* pp) 
{
   Paso_TransportProblem* fctp = fct_solver->transportproblem;
   const index_t* main_iptr=Paso_TransportProblem_borrowMainDiagonalPointer(fctp);
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const double theta = Paso_FCT_Solver_getTheta(fct_solver);
   const double omega=1./(dt* theta);
   dim_t i;
   Paso_Options options2;
   


   Paso_solve_free(fctp->iteration_matrix);
   /*
    *   fctp->iteration_matrix[i,i]=m[i]/(dt theta) -l[i,i]
    *
   */
    fct_solver->omega=omega; 
    fct_solver->dt = dt;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; ++i) {
           const double m=fctp->lumped_mass_matrix[i];
	   const double l_ii = fctp->main_diagonal_low_order_transport_matrix[i];
           fctp->iteration_matrix->mainBlock->val[main_iptr[i]] = m * omega - l_ii;
    }
    
    /* allocate preconditioner/solver */
    Paso_Options_setDefaults(&options2);
    options2.verbose = options->verbose;
    if (fct_solver->method == PASO_LINEAR_CRANK_NICOLSON  ) {
        options2.preconditioner = PASO_GS;
    } else  { 
	options2.preconditioner = PASO_JACOBI;
/* options2.preconditioner = PASO_GS; */
    }
    options2.use_local_preconditioner = FALSE;
    options2.sweeps=-1;
    
    Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
    Paso_SystemMatrix_setPreconditioner(fctp->iteration_matrix, &options2);
    Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
}

/* entry point for update proceedures */
err_t Paso_FCT_Solver_update(Paso_FCT_Solver *fct_solver, double* u, double *u_old,  Paso_Options* options, Paso_Performance *pp) 
{    
    const index_t method=fct_solver->method;
    err_t err_out = SOLVER_NO_ERROR;
    
   
    if (method == PASO_LINEAR_CRANK_NICOLSON) {
        err_out=Paso_FCT_Solver_update_LCN(fct_solver, u, u_old, options, pp);
      
    } else if (method == PASO_CRANK_NICOLSON) {
        err_out=Paso_FCT_Solver_updateNL(fct_solver, u, u_old, options, pp);
      
    } else if (method == PASO_BACKWARD_EULER) {
        err_out=Paso_FCT_Solver_updateNL(fct_solver, u, u_old, options, pp);
    } else {
        err_out = SOLVER_INPUT_ERROR;
    }
    return err_out;
  
}

/* linear crank-nicolson update */
err_t Paso_FCT_Solver_update_LCN(Paso_FCT_Solver *fct_solver, double * u, double *u_old, Paso_Options* options, Paso_Performance *pp) 
{
    double const dt = fct_solver->dt;
    dim_t sweep_max;
    double *b = fct_solver->b;
    double const RTOL = options->tolerance;
    const dim_t n=Paso_TransportProblem_getTotalNumRows(fct_solver->transportproblem); 
    Paso_SystemMatrix * iteration_matrix = fct_solver->transportproblem->iteration_matrix;
    err_t errorCode = SOLVER_NO_ERROR;
  
    Paso_Coupler_startCollect(fct_solver->u_old_coupler,u_old);
    Paso_Coupler_finishCollect(fct_solver->u_old_coupler);
    
    /* b[i]=m*u_tilde[i] = m u_old[i] + dt/2 sum_{j <> i} l_{ij}*(u_old[j]-u_old[i]) 
        * note that iteration_matrix stores the negative values of the 
        * low order transport matrix l. Therefore a=-dt*0.5 is used. */
	  
    Paso_FCT_Solver_setMuPaLu(b, fct_solver->transportproblem->lumped_mass_matrix, 
			     fct_solver->u_old_coupler, -dt*0.5, iteration_matrix); 
    /* solve for u_tilde : u_tilda = m^{-1} * b   */ 
    Paso_FCT_FluxLimiter_setU_tilda(fct_solver->flux_limiter, b);
    /* u_tilda_connector is completed */
    
    /* calculate anti-diffusive fluxes for u_tilda */
    Paso_FCT_setAntiDiffusionFlux_linearCN(fct_solver->flux_limiter->antidiffusive_fluxes, 
					   fct_solver->transportproblem, dt, 
					   fct_solver->flux_limiter->u_tilde_coupler, 
					   fct_solver->u_old_coupler);


   /* b_i += sum_{j} limitation factor_{ij} * antidiffusive_flux_{ij} */
   Paso_FCT_FluxLimiter_addLimitedFluxes_Start(fct_solver->flux_limiter);
   Paso_FCT_FluxLimiter_addLimitedFluxes_Complete(fct_solver->flux_limiter, b); 
   
   
   /* solve (m-dt/2*L) u = b in the form (omega*m-L) u = b * omega with omega*dt/2=1 */
   
   /* initial guess is u<- -u + 2*u_tilde */
   Paso_Update(n, -1., u, 2., fct_solver->flux_limiter->u_tilde);
   Paso_Scale(n, b,fct_solver->omega ); 
   sweep_max = MAX((int) (- 2 * log(RTOL)/log(2.)-0.5),1);

   if (options->verbose) {
       const double norm_u_tilde=Paso_lsup(n, fct_solver->flux_limiter->u_tilde, fct_solver->flux_limiter->mpi_info);
       printf("Paso_FCT_Solver_update_LCN: u_tilda lsup = %e (rtol = %e, max. sweeps = %d)\n",norm_u_tilde,RTOL,sweep_max); 
   }
   errorCode = Paso_Preconditioner_Smoother_solve_byTolerance( iteration_matrix,  ((Paso_Preconditioner*) (iteration_matrix->solver_p))->gs, 
						   u, b, RTOL, &sweep_max, TRUE);
   if (errorCode == PRECONDITIONER_NO_ERROR) {
      if (options->verbose) printf("Paso_FCT_Solver_update_LCN: convergence after %d Gauss-Seidel steps.\n",sweep_max);
      errorCode=SOLVER_NO_ERROR;
   } else {
      if (options->verbose) printf("Paso_FCT_Solver_update_LCN: Gauss-Seidel failed within %d stesp (rel. tolerance %e).\n",sweep_max,RTOL);
      errorCode= SOLVER_MAXITER_REACHED;
   }
   return errorCode;
   
}   

err_t Paso_FCT_Solver_updateNL(Paso_FCT_Solver *fct_solver, double* u, double *u_old, Paso_Options* options, Paso_Performance *pp) 
{
   const dim_t num_critical_rates_max=3; /* number of rates >=critical_rate accepted before divergence is triggered */
   const double critical_rate=0.95;   /* expected value of convergence rate */

   double *b = fct_solver->b;
   double *z = fct_solver->z;
   double *du = fct_solver->du;
   double const dt = fct_solver->dt; 
   Paso_TransportProblem* fctp = fct_solver->transportproblem;
   Paso_FCT_FluxLimiter* flux_limiter = fct_solver->flux_limiter;
   dim_t i;
   double norm_u_tilde, ATOL, norm_du=LARGE_POSITIVE_FLOAT, norm_du_old, rate=1.;
   err_t errorCode=SOLVER_NO_ERROR;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const double atol=options->absolute_tolerance;  
   const double rtol=options->tolerance;
   const dim_t max_m=options->iter_max;
   dim_t m=0, num_critical_rates=0 ;
  /* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
   
   bool_t converged=FALSE, max_m_reached=FALSE,diverged=FALSE;   
   options->num_iter=0;
   
   Paso_Coupler_startCollect(fct_solver->u_old_coupler,u_old);
   Paso_Coupler_finishCollect(fct_solver->u_old_coupler);
    /* prepare u_tilda and flux limiter */
    if ( fct_solver->method == PASO_BACKWARD_EULER ) {
          /* b[i]=m_i* u_old[i] */
          #pragma omp for private(i) schedule(static) 
          for (i = 0; i < n; ++i) {
               b[i]=u_old[i]* fctp->lumped_mass_matrix[i];
          } 
    } else {
       /* b[i]=m_i* u_old[i] + dt/2 sum_{j <> i} l_{ij}*(u_old[j]-u_old[i]) = m_i * u_tilde_i 
        * note that iteration_matrix stores the negative values of the 
        * low order transport matrix l. Therefore a=-dt*0.5 is used. */
       Paso_FCT_Solver_setMuPaLu(b,fctp->lumped_mass_matrix,fct_solver->u_old_coupler,-dt*0.5,fctp->iteration_matrix);
    }        
    Paso_FCT_FluxLimiter_setU_tilda(flux_limiter, b); /* u_tilda = m^{-1} b */
    /* u_tilda_connector is completed */
   
    /**********************************************************************************************************************/   
    /* calculate stopping criterium */
    norm_u_tilde=Paso_lsup(n, flux_limiter->u_tilde, flux_limiter->mpi_info);
    ATOL= rtol * norm_u_tilde + atol ;
    if (options->verbose) printf("Paso_FCT_Solver_updateNL: iteration starts u_tilda lsup = %e (abs. tol = %e)\n",norm_u_tilde,ATOL);
    /* ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */

    /* u_old is an initial guess for u*/
    Paso_Copy(n,u,u_old);
    
    while ( (!converged) && (!diverged) && (! max_m_reached) && Esys_noError()) {
         Paso_Coupler_startCollect(fct_solver->u_coupler,u);
         Paso_Coupler_finishCollect(fct_solver->u_coupler); 

         /*  set antidiffusive_flux_{ij} for u */
         if (fct_solver->method == PASO_BACKWARD_EULER) {
             Paso_FCT_setAntiDiffusionFlux_BE(fct_solver->flux_limiter->antidiffusive_fluxes, fctp, dt, fct_solver->u_coupler, fct_solver->u_old_coupler);   
         } else {         
  	     Paso_FCT_setAntiDiffusionFlux_CN(fct_solver->flux_limiter->antidiffusive_fluxes, fctp, dt, fct_solver->u_coupler, fct_solver->u_old_coupler);
         }
         /* start the calculation of the limitation factors_{fct_solver->ij} */ 
         Paso_FCT_FluxLimiter_addLimitedFluxes_Start(flux_limiter); /* uses u_tilde */  

         /*
	  * z_m[i]=b[i] - (m_i*u[i] - omega*sum_{j<>i} l_{ij} (u[j]-u[i]) ) omega = dt/2 or dt .
	  *
	  * note that iteration_matrix stores the negative values of the 
	  * low order transport matrix l. Therefore a=dt*theta is used.
	  */
	  if (fct_solver-> method == PASO_BACKWARD_EULER) {
	      Paso_FCT_Solver_setMuPaLu(z, fctp->lumped_mass_matrix, fct_solver->u_coupler, dt, fctp->iteration_matrix);
	  } else {
	      Paso_FCT_Solver_setMuPaLu(z, fctp->lumped_mass_matrix, fct_solver->u_coupler, dt/2, fctp->iteration_matrix);
	  }
   
 
	  Paso_Update(n,-1.,z,1.,b);  /* z=b-z */


	  /* z_i += sum_{j} limitation factor_{ij} * antidiffusive_flux_{ij} */
	  Paso_FCT_FluxLimiter_addLimitedFluxes_Complete(flux_limiter, z); 
   
	  /* we solve (m/omega - L ) * du = z */
	  if (fct_solver->method == PASO_BACKWARD_EULER) {
	      dim_t cntIter = options->iter_max;
	      double tol= Paso_l2(n, z, fctp->mpi_info) ;

	      if ( m ==0) {
		  tol *=0.5;
	      } else {
		  tol *= MIN(MAX(rate*rate, 1e-2), 0.5);
	      } 
	      /* use BiCGSTab with jacobi preconditioner ( m - omega * L ) */
    	      Paso_zeroes(n,du); 
	      errorCode = Paso_Solver_BiCGStab(fctp->iteration_matrix, z, du, &cntIter, &tol, pp);

	      /* errorCode =  Paso_Solver_GMRES(fctp->iteration_matrix, z, du, &cntIter, &tol, 10, 2000, pp); */
	      if (options->verbose) printf("Paso_FCT_Solver_updateNL: BiCGStab is completed after %d steps (residual =%e).\n",cntIter, tol);
	      options->num_iter+=cntIter;
	      if ( errorCode !=  SOLVER_NO_ERROR) break;
	   } else {
	      /* just use the main diagonal of (m/omega - L ) */

	      Paso_Preconditioner_Smoother_solve(fctp->iteration_matrix, ((Paso_Preconditioner*) (fctp->iteration_matrix->solver_p))->jacobi, 
					  du, z, 1, FALSE); 
	      options->num_iter++;
	   }           
	   
 
	   Paso_Update(n,1.,u,fct_solver->omega,du);
	   norm_du_old=norm_du;
           norm_du=Paso_lsup(n,du, fctp->mpi_info);
	   if (m ==0) {
                if (options->verbose) printf("Paso_FCT_Solver_updateNL: step %d: increment= %e\n",m+1, norm_du * fct_solver->omega);
	   } else {
	        if (norm_du_old > 0.) {
	            rate=norm_du/norm_du_old;
	        } else if (norm_du <= 0.) {
		     rate=0.;
		} else {
		     rate=LARGE_POSITIVE_FLOAT;
		}
                if (options->verbose) printf("Paso_FCT_Solver_updateNL: step %d: increment= %e (rate = %e)\n",m+1, norm_du * fct_solver->omega, rate);
		num_critical_rates+=( rate<critical_rate ? 0 : 1);
                max_m_reached=(m>max_m);
		diverged = (num_critical_rates >= num_critical_rates_max);
                converged=(norm_du * fct_solver->omega  <= ATOL) ; 
           }
           m++;
      } /* end of while loop */
      if (errorCode == SOLVER_NO_ERROR) {
	        if (converged) {
	            if (options->verbose) printf("Paso_FCT_Solver_updateNL: iteration is completed.\n"); 
                    errorCode=SOLVER_NO_ERROR;
                } else if (diverged) {
		    if (options->verbose) printf("Paso_FCT_Solver_updateNL: divergence.\n");
		    errorCode=SOLVER_DIVERGENCE;
		} else if (max_m_reached) {
		     if (options->verbose) printf("Paso_FCT_Solver_updateNL: maximum number of iteration steps reached.\n");
		    errorCode=SOLVER_MAXITER_REACHED;
		}
	    
    }
    return errorCode;
}


/*   
 *  AntiDiffusionFlux: 
 *
 *        f_{ij} = (m_{ij} - dt (1-theta) d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt theta d_{ij}) (u[j]-u[i])
 *         
 *     m=fc->mass matrix
 *     d=artifical diffusion matrix = L - K = - fc->iteration matrix - fc->transport matrix (away from main diagonal)
 *
 *   for CN : theta =0.5
 *   for BE : theta = 1.
 */

void Paso_FCT_setAntiDiffusionFlux_CN(Paso_SystemMatrix *flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			              const Paso_Coupler* u_coupler,  
				      const Paso_Coupler* u_old_coupler)
{
  dim_t i;
  index_t iptr_ij;
  
  const double *u    = Paso_Coupler_borrowLocalData(u_coupler);
  const double *u_old= Paso_Coupler_borrowLocalData(u_old_coupler);
  const double *remote_u=Paso_Coupler_borrowRemoteData(u_coupler);
  const double *remote_u_old=Paso_Coupler_borrowRemoteData(u_old_coupler);
  const double dt_half= dt/2;
  const Paso_SystemMatrixPattern *pattern=fct->iteration_matrix->pattern;
  const dim_t n=Paso_SystemMatrix_getTotalNumRows(fct->iteration_matrix);

  #pragma omp parallel for schedule(static) private(i, iptr_ij)
  for (i = 0; i < n; ++i) {
      const double u_i     = u[i];
      const double u_old_i = u_old[i];

      #pragma ivdep
      for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
	const index_t j      = pattern->mainPattern->index[iptr_ij];
	const double m_ij    = fct->mass_matrix->mainBlock->val[iptr_ij];
	const double d_ij    = fct->transport_matrix->mainBlock->val[iptr_ij]+fct->iteration_matrix->mainBlock->val[iptr_ij]; /* this is in fact -d_ij */
	const double u_old_j = u_old[j];
	const double u_j     = u[j];
	
	/* (m_{ij} - dt (1-theta) d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt theta d_{ij}) (u[j]-u[i]) */
	flux_matrix->mainBlock->val[iptr_ij]=(m_ij+dt_half*d_ij)*(u_old_j-u_old_i) - (m_ij-dt_half*d_ij)*(u_j-u_i);

      }
      #pragma ivdep
      for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
	const index_t j      = pattern->col_couplePattern->index[iptr_ij];
	const double m_ij    = fct->mass_matrix->col_coupleBlock->val[iptr_ij];
	const double d_ij    = fct->transport_matrix->col_coupleBlock->val[iptr_ij]+fct->iteration_matrix->col_coupleBlock->val[iptr_ij]; /* this is in fact -d_ij */
        const double u_old_j = remote_u_old[j];
	const double u_j     = remote_u[j];
	flux_matrix->col_coupleBlock->val[iptr_ij]=(m_ij+dt_half*d_ij)*(u_old_j-u_old_i)- (m_ij-dt_half*d_ij)*(u_j-u_i);
      }
  }
}

void Paso_FCT_setAntiDiffusionFlux_BE(Paso_SystemMatrix *flux_matrix,
                                      const Paso_TransportProblem* fct, 
				      const double dt,
			              const Paso_Coupler* u_coupler,  
				      const Paso_Coupler* u_old_coupler)
{
  dim_t i;
  index_t iptr_ij;
  
  const double *u=Paso_Coupler_borrowLocalData(u_coupler);
  const double *u_old= Paso_Coupler_borrowLocalData(u_old_coupler);
  const double *remote_u=Paso_Coupler_borrowRemoteData(u_coupler);
  const double *remote_u_old=Paso_Coupler_borrowRemoteData(u_old_coupler);
  const Paso_SystemMatrixPattern *pattern=fct->iteration_matrix->pattern;
  const dim_t  n=Paso_SystemMatrix_getTotalNumRows(fct->iteration_matrix);

  #pragma omp parallel for schedule(static) private(i, iptr_ij)
  for (i = 0; i < n; ++i) {
      const double u_i     = u[i];
      const double u_old_i = u_old[i];
      #pragma ivdep
      for (iptr_ij=pattern->mainPattern->ptr[i]; iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {

	const index_t j      = pattern->mainPattern->index[iptr_ij];
	const double m_ij    = fct->mass_matrix->mainBlock->val[iptr_ij];
	const double d_ij    = fct->transport_matrix->mainBlock->val[iptr_ij]+fct->iteration_matrix->mainBlock->val[iptr_ij]; /* this is in fact -d_ij */
	const double u_old_j = u_old[j];
	const double u_j     = u[j];
	
	flux_matrix->mainBlock->val[iptr_ij]=m_ij*(u_old_j-u_old_i)- (m_ij-dt*d_ij)*(u_j-u_i);
      }
      #pragma ivdep
      for (iptr_ij=pattern->col_couplePattern->ptr[i]; iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
	const index_t j      = pattern->col_couplePattern->index[iptr_ij];
	const double m_ij    = fct->mass_matrix->col_coupleBlock->val[iptr_ij]; /* this is in fact -d_ij */
	const double d_ij    = fct->transport_matrix->col_coupleBlock->val[iptr_ij]+fct->iteration_matrix->col_coupleBlock->val[iptr_ij];
        const double u_old_j = remote_u_old[j];
	const double u_j     = remote_u[j];
	
	flux_matrix->col_coupleBlock->val[iptr_ij]=m_ij*(u_old_j-u_old_i)- (m_ij-dt*d_ij)*(u_j-u_i);
      }
  }
}

/* special version of the ant-diffusive fluxes for the linear Crank-Nicolson scheme 
 * in fact this is evaluated for u = 2*u_tilde - u_old which is the predictor
 * of the solution of the the stabilized problem at time dt using the forward Euler scheme 
 * 
 *   f_{ij} = (m_{ij} - dt/2 d_{ij}) (u_old[j]-u_old[i]) - (m_{ij} + dt/2 d_{ij}) (u[j]-u[i])
 *    =  (m_{ij} - dt/2 d_{ij}) * (u_old[j]-u_old[i]) - (m_{ij} + dt/2 d_{ij}) * ( 2*(u_tilde[j]-u_tilde[i]) - (u_old[j] -u_old [i]) )
 *    =  2*  m_{ij} * ( (u_old[j]-u_tilde[j] - (u_old[i]) - u_tilde[i]) ) - dt d_{ij} * (u_tilde[j]-u_tilde[i]) 
 * 
 */
  
void Paso_FCT_setAntiDiffusionFlux_linearCN(Paso_SystemMatrix *flux_matrix,
                                            const Paso_TransportProblem* fct, 
				            const double dt,
			                    const Paso_Coupler* u_tilde_coupler,  
				            const Paso_Coupler* u_old_coupler)
{
  dim_t i;
  index_t iptr_ij;
  
  const double *u_tilde=Paso_Coupler_borrowLocalData(u_tilde_coupler);
  const double *u_old= Paso_Coupler_borrowLocalData(u_old_coupler);
  const double *remote_u_tilde=Paso_Coupler_borrowRemoteData(u_tilde_coupler);
  const double *remote_u_old=Paso_Coupler_borrowRemoteData(u_old_coupler);
  const Paso_SystemMatrixPattern *pattern=fct->iteration_matrix->pattern;
  const dim_t n=Paso_SystemMatrix_getTotalNumRows(fct->iteration_matrix);

  #pragma omp parallel for schedule(static) private(i, iptr_ij)
  for (i = 0; i < n; ++i) {
      const double u_tilde_i = u_tilde[i];
      const double u_old_i   = u_old[i];
      const double du_i      = u_tilde_i - u_old_i;
      #pragma ivdep
      for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
	
	  const index_t j      = pattern->mainPattern->index[iptr_ij];
	  const double m_ij    = fct->mass_matrix->mainBlock->val[iptr_ij];
	  const double d_ij    = fct->transport_matrix->mainBlock->val[iptr_ij]+fct->iteration_matrix->mainBlock->val[iptr_ij]; /* this is in fact -d_ij */
          const double u_tilde_j = u_tilde[j];
	  const double u_old_j = u_old[j];
	  const double du_j    = u_tilde_j - u_old_j;
	  
	  flux_matrix->mainBlock->val[iptr_ij]= 2 * m_ij * ( du_i - du_j ) - dt * d_ij * ( u_tilde_i - u_tilde_j);
      }
      #pragma ivdep
      for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
	
	const index_t j      = pattern->col_couplePattern->index[iptr_ij];
	const double m_ij    = fct->mass_matrix->col_coupleBlock->val[iptr_ij];
	const double d_ij    = fct->transport_matrix->col_coupleBlock->val[iptr_ij]+fct->iteration_matrix->col_coupleBlock->val[iptr_ij];/* this is in fact -d_ij */
        const double u_tilde_j = remote_u_tilde[j];
	const double u_old_j = remote_u_old[j];
	const double du_j    = u_tilde_j - u_old_j;	

	flux_matrix->col_coupleBlock->val[iptr_ij]= 2 * m_ij * ( du_i - du_j ) - dt * d_ij *  ( u_tilde_i - u_tilde_j); 
	
      }
  }

}

/**************************************************************/

/* Creates the low order transport matrix and stores its negative values 
 * into the iteration_matrix except for the main diagonal which is stored
 * separately.
 * If fc->iteration_matrix==NULL, fc->iteration_matrix is allocated 
 *
 * a=transport_matrix 
 * b= low_order_transport_matrix = - iteration_matrix
 * c=main diagonal low_order_transport_matrix 
 * initialize c[i] mit a[i,i] 
 *
 *    d_ij=max(0,-a[i,j],-a[j,i])  
 *    b[i,j]=-(a[i,j]+d_ij)         
 *    c[i]-=d_ij                  
 */

void Paso_FCT_setLowOrderOperator(Paso_TransportProblem * fc) {
  
  dim_t i;
  index_t iptr_ij, iptr_ji;
  const index_t*  main_iptr=Paso_TransportProblem_borrowMainDiagonalPointer(fc);

  if (fc->iteration_matrix==NULL) {
      fc->iteration_matrix=Paso_SystemMatrix_alloc(fc->transport_matrix->type,
                                                   fc->transport_matrix->pattern,
                                                   fc->transport_matrix->row_block_size,
                                                   fc->transport_matrix->col_block_size, TRUE);
  }

  if (Esys_noError()) {
      const Paso_SystemMatrixPattern *pattern=fc->iteration_matrix->pattern;
      const dim_t n=Paso_SystemMatrix_getTotalNumRows(fc->iteration_matrix);
      #pragma omp parallel for private(i, iptr_ij, iptr_ji)  schedule(static)
      for (i = 0; i < n; ++i) {
          double sum=fc->transport_matrix->mainBlock->val[main_iptr[i]];
	  
/* printf("sum[%d] = %e -> ", i, sum); */
          /* look at a[i,j] */
          for (iptr_ij=pattern->mainPattern->ptr[i];iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
              const index_t j    = pattern->mainPattern->index[iptr_ij];
              const double rtmp1 = fc->transport_matrix->mainBlock->val[iptr_ij];
	      if (j!=i) {
                 /* find entry a[j,i] */
                 #pragma ivdep
                 for (iptr_ji=pattern->mainPattern->ptr[j]; iptr_ji<pattern->mainPattern->ptr[j+1]; ++iptr_ji) {
		   
                    if ( pattern->mainPattern->index[iptr_ji] == i) {
		        const double rtmp2=fc->transport_matrix->mainBlock->val[iptr_ji];
/*
printf("a[%d,%d]=%e\n",i,j,rtmp1);
printf("a[%d,%d]=%e\n",j,i,rtmp2);
*/

                        const double d_ij=-MIN3(0.,rtmp1,rtmp2);
                        fc->iteration_matrix->mainBlock->val[iptr_ij]=-(rtmp1+d_ij);
/* printf("l[%d,%d]=%e\n",i,j,fc->iteration_matrix->mainBlock->val[iptr_ij]); */
                        sum-=d_ij;
                        break;
                    }
                 }
             }
          }
          for (iptr_ij=pattern->col_couplePattern->ptr[i];iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
              const index_t    j = pattern->col_couplePattern->index[iptr_ij];
              const double  rtmp1 = fc->transport_matrix->col_coupleBlock->val[iptr_ij];
	      /* find entry a[j,i] */
              #pragma ivdep
              for (iptr_ji=pattern->row_couplePattern->ptr[j]; iptr_ji<pattern->row_couplePattern->ptr[j+1]; ++iptr_ji) {
                    if (pattern->row_couplePattern->index[iptr_ji]==i) {
                        const double rtmp2=fc->transport_matrix->row_coupleBlock->val[iptr_ji];
                        const double d_ij=-MIN3(0.,rtmp1,rtmp2);
                        fc->iteration_matrix->col_coupleBlock->val[iptr_ij]=-(rtmp1+d_ij);
                        fc->iteration_matrix->row_coupleBlock->val[iptr_ji]=-(rtmp2+d_ij);
                        sum-=d_ij;
                        break;
                    }
              }
         }
         /* set main diagonal entry */
         fc->main_diagonal_low_order_transport_matrix[i]=sum;
/*	 printf("%e \n", sum); */
      }

  }
}

/*
 * out_i=m_i u_i + a * \sum_{j <> i} l_{ij} (u_j-u_i)
 *
 */
void Paso_FCT_Solver_setMuPaLu(double* out,
                              const double* M, 
                              const Paso_Coupler* u_coupler,
                              const double a,
                              const Paso_SystemMatrix *L)
{
  dim_t i;
  const Paso_SystemMatrixPattern *pattern = L->pattern;
  const double *u=Paso_Coupler_borrowLocalData(u_coupler);
  const double *remote_u=Paso_Coupler_borrowRemoteData(u_coupler);
  index_t iptr_ij;
  const dim_t n=Paso_SystemMatrix_getTotalNumRows(L);

  #pragma omp parallel for private(i) schedule(static)
  for (i = 0; i < n; ++i) {
      out[i]=M[i]*u[i];
  } 
  if (ABS(a)>0) {
      #pragma omp parallel for schedule(static) private(i, iptr_ij) 
      for (i = 0; i < n; ++i) {
          double sum=0;
          const double u_i=u[i];
          #pragma ivdep
  	  for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
               const index_t j=pattern->mainPattern->index[iptr_ij];
               const double l_ij=L->mainBlock->val[iptr_ij];
               sum+=l_ij*(u[j]-u_i);
	       
          }
          #pragma ivdep
  	  for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
               const index_t j=pattern->col_couplePattern->index[iptr_ij];
               const double l_ij=L->col_coupleBlock->val[iptr_ij];
               sum+=l_ij*(remote_u[j]-u_i);
          }
          out[i]+=a*sum;
      }
  }
}

/* *************************************************************************************************************************** */







