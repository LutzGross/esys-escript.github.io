
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: CFT: fluxlimiter for a transport problem
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "FluxLimiter.h"


Paso_FCT_FluxLimiter* Paso_FCT_FluxLimiter_alloc(Paso_TransportProblem *fctp) 
{
     Paso_FCT_FluxLimiter* out=NULL;
     const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
     const dim_t blockSize=Paso_TransportProblem_getBlockSize(fctp);
     
     out=MEMALLOC(1, Paso_FCT_FluxLimiter);
     if (! Esys_checkPtr(out)) {
      
            out->mpi_info = Esys_MPIInfo_getReference(fctp->mpi_info);
	    out->u_tilde=MEMALLOC(n,double);
	    Esys_checkPtr(out->u_tilde);

	    out->MQ=MEMALLOC(2*n,double);
	    Esys_checkPtr(out->MQ);

	    out->R=MEMALLOC(2*n,double);
	    Esys_checkPtr(out->R);

	    out->R_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),2*blockSize);
	    out->u_tilde_coupler=Paso_Coupler_alloc(Paso_TransportProblem_borrowConnector(fctp),blockSize);
	    out->antidiffusive_fluxes=Paso_SystemMatrix_alloc(fctp->transport_matrix->type,
						    fctp->transport_matrix->pattern,
						    fctp->transport_matrix->row_block_size,
						    fctp->transport_matrix->col_block_size, TRUE);
	    out->borrowed_lumped_mass_matrix=fctp->lumped_mass_matrix;
    }
    if (Esys_noError()) {
        return out;
    } else {
        Paso_FCT_FluxLimiter_free(out);
        return NULL;
    }
}

void Paso_FCT_FluxLimiter_free(Paso_FCT_FluxLimiter * in) 
{
   if (in!=NULL) {
       Esys_MPIInfo_free(in->mpi_info);
       Paso_SystemMatrix_free(in->antidiffusive_fluxes);
       MEMFREE(in->u_tilde);
       MEMFREE(in->MQ);
       MEMFREE(in->R);
       Paso_Coupler_free(in->R_coupler);
       Paso_Coupler_free(in->u_tilde_coupler);
       MEMFREE(in);
   }
}

/* sets the predictor u_tilda from Mu_tilda by solving M_C * u_tilda = Mu_tilda
 * and calculates the limiters QP and QN: : */
void Paso_FCT_FluxLimiter_setU_tilda(Paso_FCT_FluxLimiter* flux_limiter, const double *Mu_tilda) {
  
  const dim_t n=Paso_FCT_FluxLimiter_getTotalNumRows(flux_limiter);
  double * remote_u_tilde =NULL;
  const Paso_SystemMatrixPattern *pattern = Paso_FCT_FluxLimiter_getFluxPattern(flux_limiter);
  const double *lumped_mass_matrix = flux_limiter->borrowed_lumped_mass_matrix;
  index_t i, iptr_ij;

  #pragma omp parallel private(i)
  {
    #pragma omp for schedule(static) 
    for (i = 0; i < n; ++i) {
                const double m=lumped_mass_matrix[i];
		if (m > 0 ) {
                    flux_limiter->u_tilde[i]=Mu_tilda[i]/m;
		} else {
		    flux_limiter->u_tilde[i]=Mu_tilda[i];
		}
	   
    }        
  }
  /* distribute u_tilde: */
  Paso_Coupler_startCollect(flux_limiter->u_tilde_coupler,flux_limiter->u_tilde);
  /* 
     * calculate MQ_P[i] = lumped_mass_matrix[i] * max_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *           MQ_N[i] = lumped_mass_matrix[i] * min_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *
     */
  /* first we calculate the min and max of u_tilda in the main block 
       QP, QN are used to hold the result */ 
  #pragma omp parallel private(i)
  {
     #pragma omp  for schedule(static)
     for (i = 0; i < n; ++i) {
        if (  flux_limiter->borrowed_lumped_mass_matrix[i]> 0) { /* no constraint */    
		double u_min_i=LARGE_POSITIVE_FLOAT;
		double u_max_i=-LARGE_POSITIVE_FLOAT;
		#pragma ivdep
		for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
		    const index_t j=pattern->mainPattern->index[iptr_ij];
		    const double u_j=flux_limiter->u_tilde[j];
		    u_min_i=MIN(u_min_i,u_j);
		    u_max_i=MAX(u_max_i,u_j);
		}
		flux_limiter->MQ[2*i] = u_min_i;   
		flux_limiter->MQ[2*i+1] = u_max_i;
	    
	} else {
		    flux_limiter->MQ[2*i  ] =LARGE_POSITIVE_FLOAT;
		    flux_limiter->MQ[2*i+1] =LARGE_POSITIVE_FLOAT;
	}
     }
  } /* parallel region */
  /* complete distribute u_tilde: */
  Paso_Coupler_finishCollect(flux_limiter->u_tilde_coupler);
  remote_u_tilde=Paso_Coupler_borrowRemoteData(flux_limiter->u_tilde_coupler);

  /* now we look at the couple matrix and set the final value for QP, QN */
  #pragma omp parallel private(i, iptr_ij)
  {
      #pragma omp  for schedule(static)
      for (i = 0; i < n; ++i) {
	    if ( flux_limiter->borrowed_lumped_mass_matrix[i]> 0) { /* no constraint */ 
		const double u_i = flux_limiter->u_tilde[i];
		double u_min_i = flux_limiter->MQ[2*i];
		double u_max_i = flux_limiter->MQ[2*i+1];
		#pragma ivdep
		for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
		    const index_t j  = pattern->col_couplePattern->index[iptr_ij];
		    const double u_j = remote_u_tilde[j];
		    u_min_i=MIN(u_min_i,u_j);
		    u_max_i=MAX(u_max_i,u_j);
		}
		flux_limiter->MQ[2*i  ] = ( u_min_i-u_i ) * lumped_mass_matrix[i];  /* M_C*Q_min*/
		flux_limiter->MQ[2*i+1] = ( u_max_i-u_i ) * lumped_mass_matrix[i] ;  /* M_C*Q_max */
	    } 

     }
  } /* parallel region */
}

/* starts to update a vector (not given yet) from the antidiffusion fluxes in flux_limiter->antidiffusive_fluxes
 (needes u_tilde and Q */
void Paso_FCT_FluxLimiter_addLimitedFluxes_Start(Paso_FCT_FluxLimiter* flux_limiter) {
  
  const dim_t n=Paso_FCT_FluxLimiter_getTotalNumRows(flux_limiter);
  const Paso_SystemMatrixPattern *pattern = Paso_FCT_FluxLimiter_getFluxPattern(flux_limiter);
  const double* u_tilde = flux_limiter->u_tilde;
  const double* remote_u_tilde=Paso_Coupler_borrowRemoteData(flux_limiter->u_tilde_coupler);
  Paso_SystemMatrix * adf=flux_limiter->antidiffusive_fluxes;
  index_t iptr_ij;
  dim_t i;
 
  #pragma omp parallel private(i, iptr_ij)
  {
       #pragma omp for schedule(static)   
       for (i = 0; i < n; ++i) {
	  double R_N_i =1;
	  double R_P_i =1;
	  if ( flux_limiter->borrowed_lumped_mass_matrix[i]> 0) { /* no constraint */
		const double u_tilde_i=u_tilde[i];
		double P_P_i=0.;
		double P_N_i=0.;
		const double MQ_min=flux_limiter->MQ[2*i];
		const double MQ_max=flux_limiter->MQ[2*i+1]; 
		#pragma ivdep
		for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
		  const index_t j=pattern->mainPattern->index[iptr_ij];
		  if (i != j ) {
			  const double f_ij=adf->mainBlock->val[iptr_ij];
			  const double u_tilde_j=u_tilde[j];
			  /* pre-limiter */
			  if (f_ij * (u_tilde_j-u_tilde_i) >= 0) {
			    adf->mainBlock->val[iptr_ij]=0;
			  } else {
			      if (f_ij <=0) {
				  P_N_i+=f_ij;
			      } else {
				  P_P_i+=f_ij;
			      }
			  }
		  }
		}
		
		/* now the couple matrix: */
		#pragma ivdep
		for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
		  const index_t j=pattern->col_couplePattern->index[iptr_ij];
		  const double f_ij=adf->col_coupleBlock->val[iptr_ij];
		  const double u_tilde_j=remote_u_tilde[j];
		    /* pre-limiter */
		    if (f_ij * (u_tilde_j-u_tilde_i) >= 0) {
			adf->col_coupleBlock->val[iptr_ij]=0;
		    } else {
			if (f_ij <=0) {
			    P_N_i+=f_ij;
			} else {
			    P_P_i+=f_ij;
			}
		    }
	    }
	    /* finally the R+ and R- are calculated */

	    if (P_N_i<0) R_N_i=MIN(1,MQ_min/P_N_i);
	    if (P_P_i>0) R_P_i=MIN(1,MQ_max/P_P_i);
	  }
          flux_limiter->R[2*i]   = R_N_i;
	  flux_limiter->R[2*i+1] = R_P_i;
      }
  } /* end parallel region */
  
  /* now we kick off the distrribution of the R's */
  Paso_Coupler_startCollect(flux_limiter->R_coupler,flux_limiter->R);
} 


/* completes the exchange of the R factors and adds the weighted antdiffusion fluxxes to the residual b */
void Paso_FCT_FluxLimiter_addLimitedFluxes_Complete(Paso_FCT_FluxLimiter* flux_limiter, double* b) 
{
  const dim_t n=Paso_FCT_FluxLimiter_getTotalNumRows(flux_limiter);
  const Paso_SystemMatrixPattern *pattern = Paso_FCT_FluxLimiter_getFluxPattern(flux_limiter);
  const Paso_SystemMatrix * adf=flux_limiter->antidiffusive_fluxes;
  double *R = flux_limiter->R;
  double *remote_R=NULL;
  dim_t i;
  index_t iptr_ij;

  remote_R=Paso_Coupler_finishCollect(flux_limiter->R_coupler);

  #pragma omp parallel for schedule(static) private(i, iptr_ij)
  for (i = 0; i < n; ++i) {
    
     const double R_N_i=R[2*i];
     const double R_P_i=R[2*i+1];
     double f_i=b[i];
     
     #pragma ivdep
     for (iptr_ij=(pattern->mainPattern->ptr[i]);iptr_ij<pattern->mainPattern->ptr[i+1]; ++iptr_ij) {
         const index_t j    = pattern->mainPattern->index[iptr_ij];
         const double  f_ij = adf->mainBlock->val[iptr_ij];
	 const double  R_P_j = R[2*j+1];
	 const double  R_N_j = R[2*j];
         const double  rtmp=((f_ij >=0 ) ? MIN(R_P_i,R_N_j) : MIN(R_N_i,R_P_j) ) ;
	 f_i+=f_ij*rtmp;
     }
     #pragma ivdep
     for (iptr_ij=(pattern->col_couplePattern->ptr[i]);iptr_ij<pattern->col_couplePattern->ptr[i+1]; ++iptr_ij) {
          const index_t j    = pattern->col_couplePattern->index[iptr_ij];
          const double  f_ij = adf->col_coupleBlock->val[iptr_ij];
	  const double  R_P_j =  remote_R[2*j+1];
	  const double  R_N_j =  remote_R[2*j];
          const double rtmp=(f_ij >=0 ) ? MIN(R_P_i, R_N_j) : MIN(R_N_i, R_P_j) ;
          f_i+=f_ij*rtmp;
     }
     b[i]=f_i;
  } /* and of i loop */
} 
