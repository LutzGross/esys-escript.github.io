
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

/* Paso: Transport solver with flux corretion (L is row sum zero) 
 *        
 *   - Mv_t=Lv   v(0)=u       
 *
 *  to return v(dt)
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "FCTSolver.h"
#include "Solver.h"
#include "PasoUtil.h"

Paso_Function* Paso_FCTSolver_Function_alloc(Paso_TransportProblem *fctp, Paso_Options* options)
{
    const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
    const dim_t blockSize=Paso_TransportProblem_getBlockSize(fctp);
    Paso_Function * out=NULL;
    Paso_FCTSolver *more=NULL;
    out=MEMALLOC(1,Paso_Function);
    if (! Paso_checkPtr(out)) {
        out->kind=FCT;
        out->mpi_info=Paso_MPIInfo_getReference(fctp->mpi_info);
        out->n=n;
	out->b=TMPMEMALLOC(n,double); /*b_m */
	Paso_checkPtr(out->b);        
	out->tmp=TMPMEMALLOC(n,double); /* z_m */
	Paso_checkPtr(out->tmp); 
	
	more=MEMALLOC(1,Paso_FCTSolver);
	out->more=(void*) more;
	if (! Paso_checkPtr(more)) {
	    more->transportproblem=Paso_TransportProblem_getReference(fctp);
	    more->uTilde_n=TMPMEMALLOC(n,double);
	    Paso_checkPtr(more->uTilde_n);
	    more->QN_n=TMPMEMALLOC(n,double);
	    Paso_checkPtr(more->QN_n);
	    more->QP_n=TMPMEMALLOC(n,double);
	    Paso_checkPtr(more->QP_n);
	    more->RN_m=TMPMEMALLOC(n,double);
	    Paso_checkPtr(more->RN_m);
	    more->RP_m=TMPMEMALLOC(n,double);
	    Paso_checkPtr(more->RP_m);
	    more->QN_n_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->QP_n_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->RN_m_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->RP_m_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->uTilde_n_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->u_m_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    more->flux_matrix_m=Paso_SystemMatrix_alloc(fctp->transport_matrix->type,
						    fctp->transport_matrix->pattern,
						    fctp->transport_matrix->row_block_size,
						    fctp->transport_matrix->col_block_size, TRUE);
						    
        }
        Paso_Solver_setPreconditioner(fctp->iteration_matrix,options);
    }
    if (Paso_noError()) {
        return out;
    } else {
        Paso_FCTSolver_Function_free(out);
        return NULL;
    }
}

void Paso_FCTSolver_Function_free(Paso_Function * in) 
{
   Paso_FCTSolver *more=NULL;
   if (in!=NULL) {
       Paso_MPIInfo_free(in->mpi_info);
       MEMFREE(in->tmp);
       MEMFREE(in->b);
       more=(Paso_FCTSolver *) (in->more);
       if (NULL !=more) {
          Paso_TransportProblem_free(more->transportproblem);
          Paso_SystemMatrix_free(more->flux_matrix_m);
          TMPMEMFREE(more->uTilde_n);
          TMPMEMFREE(more->QN_n);
          TMPMEMFREE(more->QP_n);
          TMPMEMFREE(more->RN_m);
          TMPMEMFREE(more->RP_m);
          Paso_Coupler_free(more->QN_n_coupler);
          Paso_Coupler_free(more->QP_n_coupler);
          Paso_Coupler_free(more->RN_m_coupler);
          Paso_Coupler_free(more->RP_m_coupler);
          Paso_Coupler_free(more->uTilde_n_coupler);
          Paso_Coupler_free(more->u_m_coupler);
       }
       MEMFREE(in);
   }
}
double Paso_FCTSolver_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
   dim_t i, n;
   double dt_max=LARGE_POSITIVE_FLOAT, dt_max_loc;
   register double l_ii,m;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   /* extract the row sum of the advective part */
   Paso_SystemMatrix_rowSum(fctp->mass_matrix,fctp->lumped_mass_matrix);

   /* set low order transport operator */
   Paso_FCTSolver_setLowOrderOperator(fctp);
          
   if (Paso_noError()) {
        /*
         *  calculate time step size:                                           
        */
        dt_max=LARGE_POSITIVE_FLOAT;
        #pragma omp parallel private(dt_max_loc)
        {
               dt_max_loc=LARGE_POSITIVE_FLOAT;
               #pragma omp for schedule(static) private(i,l_ii,m) 
               for (i=0;i<n;++i) {
                  l_ii=fctp->main_diagonal_low_order_transport_matrix[i];
                  m=fctp->lumped_mass_matrix[i];
                  if ( (l_ii<0 && m>0.) || (l_ii>0 && m<0) ) {
                     dt_max_loc=MIN(dt_max_loc,-m/l_ii);
                  }
               }
               #pragma omp critical 
               {
                  dt_max=MIN(dt_max,dt_max_loc);
               }
        }
        #ifdef PASO_MPI
               dt_max_loc = dt_max;
               MPI_Allreduce(&dt_max_loc, &dt_max, 1, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
        #endif
	if (dt_max<LARGE_POSITIVE_FLOAT) dt_max*=fctp->dt_factor;
   }
   return dt_max;
}
/*
 * evaluates F as 
 * (m_i*omega - l_ij)*(F_j/omega)
 * =(m u^n[i] + dt*(1-theta) sum_{j <> i} l_{ij}*(u^n[j]-u^n[i]))-(m_i*u[i] - dt*theta*sum_{j<>i} l_{ij} (u[j]-u[i]) )- antidiffusion
 *
 */
err_t Paso_FCTSolver_Function_call(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp)
{
     Paso_FCTSolver *more=(Paso_FCTSolver *) (F->more);


     /* set z_m in tmp */
      Paso_FCTSolver_setUpRightHandSide(more->transportproblem,more->dt,arg,more->u_m_coupler,F->tmp,more->flux_matrix_m,
				     	more->uTilde_n_coupler,F->b, 
                                        more->QN_n_coupler,more->QP_n_coupler,more->RN_m,
					more->RN_m_coupler,more->RP_m,more->RP_m_coupler,pp);
      /* apply preconditoner do get du = */
      Paso_Solver_solvePreconditioner(more->transportproblem->iteration_matrix,value,F->tmp);
      return NO_ERROR;
}

err_t Paso_FCTSolver_solve(Paso_Function* F, double* u, double dt, Paso_Options* options, Paso_Performance *pp) 
{
   double norm_u_tilde, ATOL, norm_u, norm_du, *du=NULL;
   err_t errorCode=SOLVER_NO_ERROR;
   Paso_FCTSolver *more=(Paso_FCTSolver *) (F->more);
   const double omega=1./dt*( more->transportproblem->useBackwardEuler ? 1. : 2.);
   Paso_TransportProblem* fctp=more->transportproblem;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const double atol=options->absolute_tolerance;  
   const double rtol=options->tolerance;
   const dim_t max_m=options->iter_max;
   dim_t m=0;
   bool_t converged=FALSE, max_m_reached=FALSE;

   /* tolerance? */
    more->dt=dt;
    Paso_FCTSolver_setUp(fctp,dt,u, F->b,more->uTilde_n,more->uTilde_n_coupler,more->QN_n,more->QN_n_coupler,more->QP_n,more->QP_n_coupler,
			 options,pp);
			 
			 
    norm_u_tilde=Paso_lsup(n,more->uTilde_n,F->mpi_info);
    ATOL= rtol * norm_u_tilde + atol ;
    if (options->verbose) printf("Paso_FCTSolver_solve: iteration starts u_tilda lsup = %e (abs. tol = %e)\n",norm_u_tilde,ATOL);
    if  (more->transportproblem->useBackwardEuler && FALSE) {
        options->absolute_tolerance = ATOL/omega;
        options->tolerance=0;
        errorCode=Paso_Solver_NewtonGMRES(F,u,options,pp);
        options->absolute_tolerance=atol;
        options->tolerance=rtol;
    } else {
            m=0;
            converged=FALSE;
            max_m_reached=FALSE;
	    du=MEMALLOC(n,double);
	    if (Paso_checkPtr(du)) {
                   errorCode=SOLVER_MEMORY_ERROR;
	    } else {
	         /* tolerance? */
	         
                 while ( (!converged) && (! max_m_reached) && Paso_noError()) {
	            errorCode=Paso_FCTSolver_Function_call(F,du, u, pp);
                    Paso_Update(n,1.,u,omega,du);
                     norm_u=Paso_lsup(n,u, fctp->mpi_info);
                     norm_du=Paso_lsup(n,du, fctp->mpi_info);
                     if (options->verbose) printf("Paso_FCTSolver_solve: step %d: increment= %e (tolerance = %e)\n",m+1, norm_du*omega, ATOL);
                     max_m_reached=(m>max_m);
                     converged=(norm_du*omega <= ATOL) ;
                     m++;
                 }
	    }
	    MEMFREE(du);
	    if (errorCode == SOLVER_NO_ERROR) {
	        if (converged) {
	            if (options->verbose) printf("Paso_FCTSolver_solve: iteration is completed.\n"); 
                    errorCode=SOLVER_NO_ERROR;
                } else if (max_m_reached) {
		    errorCode=SOLVER_MAXITER_REACHED;
		}
	    }
    }
    return errorCode;
}



err_t Paso_FCTSolver_setUpRightHandSide(Paso_TransportProblem* fctp, const double dt, const double *u_m, 
					Paso_Coupler* u_m_coupler,  double * z_m,
					Paso_SystemMatrix* flux_matrix, Paso_Coupler* uTilde_coupler, const double *b, 
					Paso_Coupler* QN_coupler, Paso_Coupler* QP_coupler,
					double *RN_m, Paso_Coupler* RN_m_coupler, double* RP_m, Paso_Coupler* RP_m_coupler, 
					Paso_Performance* pp)
{
   const double omega=1./dt* (fctp->useBackwardEuler ? 1. : 2.);
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   /* distribute u */
   Paso_Coupler_startCollect(u_m_coupler,u_m);
   Paso_Coupler_finishCollect(u_m_coupler);
   /*
    *  set the ant diffusion fluxes:
    *
    */
   Paso_FCTSolver_setAntiDiffusionFlux(dt,fctp,flux_matrix,u_m_coupler);
   /*
    *  apply pre flux-correction: f_{ij}:=0 if f_{ij}*(\tilde{u}[i]- \tilde{u}[j])<=0
    *
    */
   Paso_FCTSolver_applyPreAntiDiffusionCorrection(flux_matrix,uTilde_coupler);
   /*
    *  set flux limiters RN_m,RP_m
    *
    */
   Paso_FCTSolver_setRs(flux_matrix,fctp->lumped_mass_matrix,QN_coupler,QP_coupler,RN_m,RP_m);
   Paso_Coupler_startCollect(RN_m_coupler,RN_m);
   Paso_Coupler_startCollect(RP_m_coupler,RP_m);
    /*
     * z_m[i]=b[i] - (m_i*u_m[i] - dt*theta*sum_{j<>i} l_{ij} (u_m[j]-u_m[i]) )
     *
     * note that iteration_matrix stores the negative values of the 
     * low order transport matrix l therefore a=dt*theta is used.
     */
     Paso_FCTSolver_setMuPaLu(z_m,fctp->lumped_mass_matrix, u_m_coupler,1./omega,fctp->iteration_matrix);
 /*{
 int kk;
 for (kk=0;kk<n;kk++) printf("z_m %d : %e %e -> %e\n",kk,z_m[kk], b[kk],z_m[kk]-b[kk]);
 }
*/
   /* z_m=b-z_m */
   Paso_Update(n,-1.,z_m,1.,b);

   Paso_Coupler_finishCollect(RN_m_coupler);
   Paso_Coupler_finishCollect(RP_m_coupler);
   /* add corrected fluxes into z_m */
   Paso_FCTSolver_addCorrectedFluxes(z_m,flux_matrix,RN_m_coupler,RP_m_coupler); 
   return SOLVER_NO_ERROR;
}

void Paso_FCTSolver_setUp(Paso_TransportProblem* fctp, const double dt, const double *u, double* b, double* uTilde, 
                          Paso_Coupler* uTilde_coupler, double *QN, Paso_Coupler* QN_coupler, double *QP, Paso_Coupler* QP_coupler,
                          Paso_Options* options, Paso_Performance* pp)
{
   dim_t i;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   const double omega=1./dt* (fctp->useBackwardEuler ? 1. : 2.);
   register double m, u_tilde_i, rtmp4; 
   /* distribute u */
   Paso_Coupler_startCollect(fctp->u_coupler,u);
   Paso_Coupler_finishCollect(fctp->u_coupler);
   /*
    * b^n[i]=m u^n[i] + dt*(1-theta) sum_{j <> i} l_{ij}*(u^n[j]-u^n[i])
    *
    * note that iteration_matrix stores the negative values of the 
    * low order transport matrix l therefore a=-dt*(1-fctp->theta) is used.
    *
    */
    if (fctp->useBackwardEuler) { 
       Paso_FCTSolver_setMuPaLu(b,fctp->lumped_mass_matrix,fctp->u_coupler,0.,fctp->iteration_matrix);
    } else {
       Paso_FCTSolver_setMuPaLu(b,fctp->lumped_mass_matrix,fctp->u_coupler,-dt*0.5,fctp->iteration_matrix);
    }      
   /*
    *   uTilde[i]=b[i]/m[i]
    *
    *   fctp->iteration_matrix[i,i]=m[i]/(dt theta) -l[i,i]
    *
    */

    Paso_solve_free(fctp->iteration_matrix);     
    #pragma omp parallel for private(i,m,u_tilde_i,rtmp4)
    for (i = 0; i < n; ++i) {
           m=fctp->lumped_mass_matrix[i];
           u_tilde_i=b[i]/m;
           rtmp4=m*omega-fctp->main_diagonal_low_order_transport_matrix[i];
           fctp->iteration_matrix->mainBlock->val[fctp->main_iptr[i]]=rtmp4;
           uTilde[i]=u_tilde_i;
/* printf("uTilde %d %e and %e : %e : %e %e\n",i,u_tilde_i, b[i], rtmp4, m, fctp->main_diagonal_low_order_transport_matrix[i]); */
    }
    Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
    Paso_Solver_setPreconditioner(fctp->iteration_matrix,options);
    Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
        
    /* distribute uTilde: */
    Paso_Coupler_startCollect(uTilde_coupler,uTilde);
    Paso_Coupler_finishCollect(uTilde_coupler);
    /* 
     * calculate QP[i] max_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *           QN[i] min_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *
     */
     Paso_FCTSolver_setQs(uTilde_coupler,QN,QP,fctp->iteration_matrix);
     Paso_Coupler_startCollect(QN_coupler,QN);
     Paso_Coupler_startCollect(QP_coupler,QP);
     Paso_Coupler_finishCollect(QN_coupler);
     Paso_Coupler_finishCollect(QP_coupler);
}