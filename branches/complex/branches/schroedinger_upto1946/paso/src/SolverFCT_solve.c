
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


/**************************************************************/

/* Paso: Flux correction transport solver
 *
 * solves Mu_t=Ku+q
 *        u(0) >=0
 *
 * Warning: the program assums sum_{j} k_{ij}=0!!!
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "SolverFCT.h"
#include "PasoUtil.h"
#include "escript/blocktimer.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif

/*
 * inserts the source term into the problem
 */
void Paso_FCT_setSource(const dim_t n,const double *source, double* sourceN, double* sourceP) 
{
   dim_t i;
   register double rtmp;
   /* 
    * seperate source into positive and negative part:
    */
   #pragma omp parallel for private(i,rtmp)
   for (i = 0; i < n; ++i) {
       rtmp=source[i];
       if (rtmp <0) {
          sourceN[i]=-rtmp;
          sourceP[i]=0;
       } else {
          sourceN[i]= 0;
          sourceP[i]= rtmp;
       }
   }
}

err_t Paso_FCT_setUpRightHandSide(Paso_FCTransportProblem* fctp, const double dt, const double *u_m, Paso_Coupler* u_m_coupler,  double * z_m,
                                  Paso_SystemMatrix* flux_matrix, Paso_Coupler* uTilde_coupler, const double *b, 
                                  Paso_Coupler* QN_coupler, Paso_Coupler* QP_coupler,
                                  double *RN_m, Paso_Coupler* RN_m_coupler, double* RP_m, Paso_Coupler* RP_m_coupler, const double *sourceN, 
                                  Paso_Performance* pp)
{
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   /* distribute u */
   Paso_Coupler_startCollect(u_m_coupler,u_m);
   Paso_Coupler_finishCollect(u_m_coupler);
   /*
    *  set the ant diffusion fluxes:
    *
    */
   Paso_FCTransportProblem_setAntiDiffusionFlux(dt,fctp,flux_matrix,u_m_coupler);
   /*
    *  apply pre flux-correction: f_{ij}:=0 if f_{ij}*(\tilde{u}[i]- \tilde{u}[j])<=0
    *
    */
   Paso_FCTransportProblem_applyPreAntiDiffusionCorrection(flux_matrix,uTilde_coupler);
   /*
    *  set flux limiters RN_m,RP_m
    *
    */
   Paso_FCTransportProblem_setRs(flux_matrix,fctp->lumped_mass_matrix,QN_coupler,QP_coupler,RN_m,RP_m);
   Paso_Coupler_startCollect(RN_m_coupler,RN_m);
   Paso_Coupler_startCollect(RP_m_coupler,RP_m);
    /*
     * z_m[i]=b[i] - (m_i*u_m[i] - dt*theta*sum_{j<>i} l_{ij} (u_m[j]-u_m[i]) + dt q^-[i])
     *
     * note that iteration_matrix stores the negative values of the 
     * low order transport matrix l therefore a=dt*fctp->theta is used.
     */
   Paso_SolverFCT_setMuPaLuPbQ(z_m,fctp->lumped_mass_matrix, u_m_coupler,dt*fctp->theta,fctp->iteration_matrix,dt,sourceN);
   /* z_m=b-z_m */
   Paso_Update(n,-1.,z_m,1.,b);

   Paso_Coupler_finishCollect(RN_m_coupler);
   Paso_Coupler_finishCollect(RP_m_coupler);
   /* add corrected fluxes into z_m */
   Paso_FCTransportProblem_addCorrectedFluxes(z_m,flux_matrix,RN_m_coupler,RP_m_coupler); 
   return NO_ERROR;
}

void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options) {
   const dim_t FAILURES_MAX=5;
   dim_t m, i_substeps, Failed, i;
   double *z_m=NULL, *b_n=NULL, *sourceP=NULL, *sourceN=NULL, *uTilde_n=NULL, *QN_n=NULL, *QP_n=NULL, *RN_m=NULL, *RP_m=NULL, *du_m=NULL;
   Paso_Coupler *QN_n_coupler=NULL, *QP_n_coupler=NULL, *RN_m_coupler=NULL, *RP_m_coupler=NULL, *uTilde_n_coupler=NULL, *u_m_coupler=NULL;
   Paso_SystemMatrix *flux_matrix_m=NULL;
   double dt_max, dt2,t, norm_u_m, omega, norm_du_m;
   register double mass, rtmp;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   dim_t blockSize=Paso_FCTransportProblem_getBlockSize(fctp);
   const double atol=options->absolute_tolerance;  
   const double rtol=options->tolerance;
   const dim_t max_m=options->iter_max;
   Paso_Performance pp;
   bool_t converged=FALSE, max_m_reached=FALSE;
   if (dt<=0.) {
       Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
   }
   dt_max=Paso_FCTransportProblem_getSafeTimeStepSize(fctp);
   /* 
    *  allocate memory
    *
    */
   z_m=TMPMEMALLOC(n,double);
   Paso_checkPtr(z_m);
   du_m=TMPMEMALLOC(n,double);
   Paso_checkPtr(du_m);
   b_n=TMPMEMALLOC(n,double);
   Paso_checkPtr(b_n);
   sourceP=TMPMEMALLOC(n,double);
   Paso_checkPtr(sourceP);
   sourceN=TMPMEMALLOC(n,double);
   Paso_checkPtr(sourceN);
   uTilde_n=TMPMEMALLOC(n,double);
   Paso_checkPtr(uTilde_n);
   QN_n=TMPMEMALLOC(n,double);
   Paso_checkPtr(QN_n);
   QP_n=TMPMEMALLOC(n,double);
   Paso_checkPtr(QP_n);
   RN_m=TMPMEMALLOC(n,double);
   Paso_checkPtr(RN_m);
   RP_m=TMPMEMALLOC(n,double);
   Paso_checkPtr(RP_m);
   QN_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   QP_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   RN_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   RP_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   uTilde_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   u_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   flux_matrix_m=Paso_SystemMatrix_alloc(fctp->transport_matrix->type,
                                                 fctp->transport_matrix->pattern,
                                                 fctp->transport_matrix->row_block_size,
                                                 fctp->transport_matrix->col_block_size);

   if (Paso_noError()) {
       /*
        *    Preparation:
        *
        */
       Paso_FCT_setSource(n, source, sourceN, sourceP);
       /*
        * let the show begin!!!!
        *
        */
        t=0;
        i_substeps=0;
        Paso_Copy(n,u,fctp->u);
        norm_u_m=Paso_lsup(n,u, fctp->mpi_info);
        /* while(i_substeps<n_substeps && Paso_noError()) { */
        if (fctp->dt_max < LARGE_POSITIVE_FLOAT) {
            dt2=MIN(dt_max,dt);
        } else {
             dt2=dt;
        }
        while(t<dt && Paso_noError()) {
            printf("substep step %d at t=%e (step size= %e)\n",i_substeps+1,t+dt2,dt2);
            Paso_FCT_setUp(fctp,dt2,sourceN,sourceP,b_n,uTilde_n,uTilde_n_coupler,QN_n,QN_n_coupler,QP_n,QP_n_coupler,
                           options,&pp);
            /* now the iteration starts */
            m=0;
            converged=FALSE;
            max_m_reached=FALSE;
            /* tolerance? */
            while ( (!converged) && (! max_m_reached) && Paso_noError()) {
                    /* set z_m */
                    Paso_FCT_setUpRightHandSide(fctp,dt2,u,u_m_coupler,z_m,flux_matrix_m,uTilde_n_coupler,b_n, 
                                  QN_n_coupler,QP_n_coupler,RN_m,RN_m_coupler,RP_m,RP_m_coupler,sourceN,&pp);
                    /* 
                     * now we solve the linear system to get the correction dt:
                     *
                     */
                     if (fctp->theta > 0) {
                          omega=1./(dt2*fctp->theta);
                          Paso_Solver_solvePreconditioner(fctp->iteration_matrix,du_m,z_m);  
                          Paso_Update(n,1.,u,omega,du_m);
                     } else {
                          omega=1;
                          #pragma omp parallel for private(i,mass,rtmp)
                          for (i = 0; i < n; ++i) {
                              mass=fctp->lumped_mass_matrix[i];
                              if (ABS(mass)>0.) {
                                  rtmp=z_m[i]/mass;
                              } else {
                                  rtmp=0;
                              }
                              du_m[i]=rtmp;
                              u[i]+=rtmp;
                          }
                   }
                   norm_u_m=Paso_lsup(n,u, fctp->mpi_info);
                   norm_du_m=Paso_lsup(n,du_m, fctp->mpi_info)*omega;
                   printf("iteration step %d completed: norm increment= %e (tolerance = %e)\n",m+1, norm_du_m, rtol * norm_u_m + atol);

                   max_m_reached=(m>max_m);
                   converged=(norm_du_m <= rtol * norm_u_m + atol);
                   m++;
            }
            if (converged) {
                    Failed=0;
                    /* #pragma omp parallel for schedule(static) private(i) */
                    Paso_Copy(n,fctp->u,u);
                    i_substeps++;
                    t+=dt2;
                    if (fctp->dt_max < LARGE_POSITIVE_FLOAT) {
                       dt2=MIN3(dt_max,dt2*1.5,dt-t);
                    } else {
                       dt2=MIN(dt2*1.5,dt-t);
                    }
            } else if (max_m_reached) {
                    /* if FAILURES_MAX failures in a row: give up */
                    if (Failed > FAILURES_MAX) {
                       Paso_setError(VALUE_ERROR,"Paso_SolverFCT_solve: no convergence after time step reduction.");
                    } else {
                       printf("no convergence in Paso_Solver_NewtonGMRES: Trying smaller time step size.");
                       dt2=dt*0.5;
                       Failed++;
                    }
            }
        
      }
      /* 
       *  clean-up:
       *
       */
      MEMFREE(z_m);
      MEMFREE(du_m);
      TMPMEMFREE(b_n);
      Paso_SystemMatrix_free(flux_matrix_m);
      TMPMEMFREE(sourceP);
      TMPMEMFREE(sourceN);
      TMPMEMFREE(uTilde_n);
      TMPMEMFREE(QN_n);
      TMPMEMFREE(QP_n);
      TMPMEMFREE(RN_m);
      TMPMEMFREE(RP_m);
      Paso_Coupler_free(QN_n_coupler);
      Paso_Coupler_free(QP_n_coupler);
      Paso_Coupler_free(RN_m_coupler);
      Paso_Coupler_free(RP_m_coupler);
      Paso_Coupler_free(uTilde_n_coupler);
      Paso_Coupler_free(u_m_coupler);
      options->absolute_tolerance=atol;
      options->tolerance=rtol;
   }
}
double Paso_FCTransportProblem_getSafeTimeStepSize(Paso_FCTransportProblem* fctp)
{
   dim_t i, n;
   double dt_max, dt_max_loc;
   register double l_ii,m;
   n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   if (! fctp->valid_matrices) {
        fctp->dt_max=LARGE_POSITIVE_FLOAT;
        /* extract the row sum of the advective part */
        Paso_SystemMatrix_rowSum(fctp->mass_matrix,fctp->lumped_mass_matrix);

        /* set low order transport operator */
        Paso_FCTransportProblem_setLowOrderOperator(fctp);
        /*
         *  calculate time step size:                                           
        */
        dt_max=LARGE_POSITIVE_FLOAT;
        if (fctp->theta < 1.) {
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
            if (dt_max<LARGE_POSITIVE_FLOAT) dt_max*=1./(1.-fctp->theta);
         }
         if (dt_max <= 0.)  {
            Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
         } else {
            if (dt_max<LARGE_POSITIVE_FLOAT)
               printf("maximum time step size is %e (theta = %e).\n",dt_max,fctp->theta);   
        }
        fctp->dt_max=dt_max;
        fctp->valid_matrices=Paso_noError();
   }
   return fctp->dt_max;
}
void Paso_FCT_setUp(Paso_FCTransportProblem* fctp, const double dt, const double *sourceN, const double *sourceP, double* b, double* uTilde, 
                     Paso_Coupler* uTilde_coupler, double *QN, Paso_Coupler* QN_coupler, double *QP, Paso_Coupler* QP_coupler,
                     Paso_Options* options, Paso_Performance* pp)
{
   dim_t i;
   const dim_t n=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   double omega, factor;
   register double m, u_tilde_i, rtmp4; 
   /* distribute u */
   Paso_Coupler_startCollect(fctp->u_coupler,fctp->u);
   Paso_Coupler_finishCollect(fctp->u_coupler);
   /*
    * b^n[i]=m u^n[i] + dt*(1-theta) sum_{j <> i} l_{ij}*(u^n[j]-u^n[i]) + dt*sourceP[i]
    *
    * note that iteration_matrix stores the negative values of the 
    * low order transport matrix l therefore a=-dt*(1-fctp->theta) is used.
    *
    */
    Paso_SolverFCT_setMuPaLuPbQ(b,fctp->lumped_mass_matrix,fctp->u_coupler,
                               -dt*(1-fctp->theta),fctp->iteration_matrix,dt,sourceP);
   /*
    *   uTilde[i]=b[i]/m[i]
    *
    *   fctp->iteration_matrix[i,i]=m[i]/(dt theta) + \frac{1}{\theta} \frac{q^-[i]}-l[i,i]
    *
    */
    if (fctp->theta > 0) {
         Paso_solve_free(fctp->iteration_matrix);
         omega=1./(dt*fctp->theta);
         factor=dt*omega;
         #pragma omp parallel for private(i,m,u_tilde_i,rtmp4)
         for (i = 0; i < n; ++i) {
              m=fctp->lumped_mass_matrix[i];
              if (ABS(m)>0.) {
                 u_tilde_i=b[i]/m;
              } else {
                 u_tilde_i=fctp->u[i];
              }
              rtmp4=m*omega-fctp->main_diagonal_low_order_transport_matrix[i];
              if (ABS(u_tilde_i)>0) rtmp4+=sourceN[i]*factor/u_tilde_i;
              fctp->iteration_matrix->mainBlock->val[fctp->main_iptr[i]]=rtmp4;
              uTilde[i]=u_tilde_i;
         }
         Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
         Paso_Solver_setPreconditioner(fctp->iteration_matrix,options);
         Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
    } else {
        #pragma omp parallel for private(i,m,u_tilde_i)
        for (i = 0; i < n; ++i) {
            m=fctp->lumped_mass_matrix[i];
            if (ABS(m)>0.) {
                u_tilde_i=b[i]/m;
            } else {
                u_tilde_i=fctp->u[i];
            }
            uTilde[i]=u_tilde_i;
        }
        /* no update of iteration_matrix required! */
    } /* end (fctp->theta > 0) */

    /* distribute uTilde: */
    Paso_Coupler_startCollect(uTilde_coupler,uTilde);
    Paso_Coupler_finishCollect(uTilde_coupler);
    /* 
     * calculate QP[i] max_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *           QN[i] min_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
     *
     */
     Paso_SolverFCT_setQs(uTilde_coupler,QN,QP,fctp->iteration_matrix);
     Paso_Coupler_startCollect(QN_coupler,QN);
     Paso_Coupler_startCollect(QP_coupler,QP);
     Paso_Coupler_finishCollect(QN_coupler);
     Paso_Coupler_finishCollect(QP_coupler);
}
