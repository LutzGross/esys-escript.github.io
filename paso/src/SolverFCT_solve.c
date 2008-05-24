/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007,2008 by University of Queensland
 *
 *                http://esscc.uq.edu_m.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
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

/* Author: l.gross@uq.edu_m.au */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "SolverFCT.h"
#include "escript/blocktimer.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif

double Paso_FCTransportProblem_getSafeTimeStepSize(Paso_FCTransportProblem* fctp)
{
   dim_t i, n_rows;
   double dt_max, dt_max_loc;
   register double rtmp1,rtmp2;
   n_rows=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
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
               #pragma omp for schedule(static) private(i,rtmp1,rtmp2) 
               for (i=0;i<n_rows;++i) {
                    rtmp1=fctp->main_diagonal_low_order_transport_matrix[i];
                    rtmp2=fctp->lumped_mass_matrix[i];
                    if ( (rtmp1<0 && rtmp2>0.) || (rtmp1>0 && rtmp2<0.) ) {
                        dt_max_loc=MIN(dt_max_loc,-rtmp2/rtmp1);
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
           if (dt_max<LARGE_POSITIVE_FLOAT) printf("maximum time step size is %e (theta = %e).\n",dt_max,fctp->theta);   
        }
        fctp->dt_max=dt_max;
        fctp->valid_matrices=Paso_noError();
   }
   return fctp->dt_max;
}




void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options) {
   index_t i, j;
   int n_substeps,n, m;
   double dt_max, omega, dt2,t;
   double local_norm[2],norm[2],local_norm_u,local_norm_du,norm_u,norm_du, tolerance;
   register double rtmp1,rtmp2,rtmp3,rtmp4, rtmp;
   double *b_n=NULL, *sourceP=NULL, *sourceN=NULL, *uTilde_n=NULL, *QN_n=NULL, *QP_n=NULL, *RN_m=NULL, *RP_m=NULL, *z_m=NULL, *du_m=NULL,*u_m=NULL;
   Paso_Coupler* QN_n_coupler=NULL, *QP_n_coupler=NULL, *uTilde_n_coupler=NULL, *RN_m_coupler=NULL, *RP_m_coupler=NULL, *u_m_coupler=NULL;
   Paso_SystemMatrix *flux_matrix=NULL;
   dim_t n_rows=Paso_SystemMatrix_getTotalNumRows(fctp->transport_matrix);
   bool_t converged;
   dim_t max_m=500;
   double inner_tol=1.e-2;
   dim_t blockSize=Paso_FCTransportProblem_getBlockSize(fctp);
   if (dt<=0.) {
       Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
   }
   dt_max=Paso_FCTransportProblem_getSafeTimeStepSize(fctp);
   /* 
    *  allocate memory
    *
    */
   u_m=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(u_m);
   b_n=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(b_n);
   sourceP=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(sourceP);
   sourceN=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(sourceN);
   uTilde_n=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(uTilde_n);
   QN_n=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(QN_n);
   QP_n=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(QP_n);
   RN_m=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(RN_m);
   RP_m=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(RP_m);
   z_m=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(z_m);
   du_m=TMPMEMALLOC(n_rows,double);
   Paso_checkPtr(du_m);
   QN_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   QP_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   RN_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   RP_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   uTilde_n_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   u_m_coupler=Paso_Coupler_alloc(Paso_FCTransportProblem_borrowConnector(fctp),blockSize);
   flux_matrix=Paso_SystemMatrix_alloc(fctp->transport_matrix->type,
                                       fctp->transport_matrix->pattern,
                                       fctp->transport_matrix->row_block_size,
                                       fctp->transport_matrix->col_block_size);
   if (Paso_noError()) {
       /*
        *    Preparation:
        *
        */
    
       /* decide on substepping */
       if (fctp->dt_max < LARGE_POSITIVE_FLOAT) {
          n_substeps=ceil(dt/dt_max);
       } else {
          n_substeps=1;
       }
       dt2=dt/n_substeps;
       printf("%d time steps of size is %e (theta = %e, dt_max=%e).\n",n_substeps, dt2,fctp->theta, dt_max);
       /* 
 	* seperate source into positive and negative part:
 	*/
        #pragma omp parallel for private(i,rtmp)
        for (i = 0; i < n_rows; ++i) {
          rtmp=source[i];
          if (rtmp <0) {
             sourceN[i]=-rtmp;
             sourceP[i]=0;
          } else {
             sourceN[i]= 0;
             sourceP[i]= rtmp;
          }
        }
        u_m_coupler->data=u_m;
/*        for (i = 0; i < n_rows; ++i) printf("%d : %e \n",i,source[i],sourceP[i],sourceN[i]) */
        /*
         * let the show begin!!!!
         *
         */
        t=dt2;
        n=0;
        tolerance=options->tolerance;
        while(n<n_substeps && Paso_noError()) {

            printf("substep step %d at t=%e\n",n+1,t);
            /*
             * b^n[i]=m u^n[i] + dt2*(1-theta) sum_{j <> i} l_{ij}*(u^n[j]-u^n[i]) + dt2*sourceP[i]
             *
             * note that iteration_matrix stores the negative values of the 
             * low order transport matrix l therefore a=-dt2*(1-fctp->theta) is used.
             *
             */
             Paso_SolverFCT_setMuPaLuPbQ(b_n,fctp->lumped_mass_matrix,fctp->u_coupler,
                                          -dt2*(1-fctp->theta),fctp->iteration_matrix,dt2,sourceP);
             /*
 	      *   uTilde_n[i]=b[i]/m[i]
 	      *
 	      *   fctp->iteration_matrix[i,i]=m[i]/(dt2 theta) + \frac{1}{\theta} \frac{q^-[i]}-l[i,i]
 	      *
 	      */
             if (fctp->theta > 0) {
                 Paso_solve_free(fctp->iteration_matrix);
                 omega=1./(dt2*fctp->theta);
                 rtmp2=dt2*omega;
                 #pragma omp parallel for private(i,rtmp,rtmp3,rtmp4)
                 for (i = 0; i < n_rows; ++i) {
                      rtmp=fctp->lumped_mass_matrix[i];
                      if (ABS(rtmp)>0.) {
                         rtmp3=b_n[i]/rtmp;
                      } else {
                         rtmp3=fctp->u[i];
                      }
                      rtmp4=rtmp*omega-fctp->main_diagonal_low_order_transport_matrix[i];
                      if (ABS(rtmp3)>0) rtmp4+=sourceN[i]*rtmp2/rtmp3;
                      fctp->iteration_matrix->mainBlock->val[fctp->main_iptr[i]]=rtmp4;
                      uTilde_n[i]=rtmp3;
                 }
             } else {
                 #pragma omp parallel for private(i,rtmp,rtmp3)
                 for (i = 0; i < n_rows; ++i) {
                      rtmp=fctp->lumped_mass_matrix[i];
                      if (ABS(rtmp)>0.) {
                         rtmp3=b_n[i]/rtmp;
                      } else {
                         rtmp3=fctp->u[i];
                      }
                      uTilde_n[i]=rtmp3;
                 }
                 omega=1.;
                 /* no update of iteration_matrix required! */
             } /* end (fctp->theta > 0) */

             /* distribute uTilde_n: */
             Paso_Coupler_startCollect(uTilde_n_coupler,uTilde_n);
             Paso_Coupler_finishCollect(uTilde_n_coupler);
             /* 
              * calculate QP_n[i] max_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
              *           QN_n[i] min_{j} (\tilde{u}[j]- \tilde{u}[i] ) 
              *
              */
              Paso_SolverFCT_setQs(uTilde_n_coupler,QN_n,QP_n,fctp->iteration_matrix);
              Paso_Coupler_startCollect(QN_n_coupler,QN_n);
              Paso_Coupler_startCollect(QP_n_coupler,QP_n);
              Paso_Coupler_finishCollect(QN_n_coupler);
              Paso_Coupler_finishCollect(QP_n_coupler);
              /*
               * now we enter the iteration on a time-step:
               *
               */
               Paso_Coupler_copyAll(fctp->u_coupler, u_m_coupler);
/* REPLACE BY JACOBEAN FREE NEWTON */
               m=0;
               converged=FALSE;
               while ( (!converged) && (m<max_m) && Paso_noError()) {
                    printf("iteration step %d\n",m+1);
                    /*
                     *  set the ant diffusion fluxes:
                     *
                    */
                    Paso_FCTransportProblem_setAntiDiffusionFlux(dt2,fctp,flux_matrix,u_m_coupler);
                    /*
                     *  apply pre flux-correction: f_{ij}:=0 if f_{ij}*(\tilde{u}[i]- \tilde{u}[j])<=0
                     *
                     */
                    Paso_FCTransportProblem_applyPreAntiDiffusionCorrection(flux_matrix,uTilde_n_coupler);
                    /*
                     *  set flux limiters RN_m,RP_m
                     *
                     */
                    Paso_FCTransportProblem_setRs(flux_matrix,fctp->lumped_mass_matrix,QN_n_coupler,QP_n_coupler,RN_m,RP_m);
                    Paso_Coupler_startCollect(RN_m_coupler,RN_m);
                    Paso_Coupler_startCollect(RP_m_coupler,RP_m);
                    /*
                     * z_m[i]=b_n[i] - (m_i*u[i] - dt2*theta*sum_{j<>i} l_{ij} (u[j]-u[i]) + dt2 q^-[i])
                     *
                     * note that iteration_matrix stores the negative values of the 
                     * low order transport matrix l therefore a=dt2*fctp->theta is used.
                     */

                    Paso_SolverFCT_setMuPaLuPbQ(z_m,fctp->lumped_mass_matrix, u_m_coupler,
                                                dt2*fctp->theta,fctp->iteration_matrix,dt2,sourceN);
/* for (i = 0; i < n_rows; ++i) {
* printf("%d: %e %e\n",i,z_m[i], b_n[i]);
* }
*/
                    #pragma omp parallel for private(i)
                    for (i = 0; i < n_rows; ++i) z_m[i]=b_n[i]-z_m[i];

                    Paso_Coupler_finishCollect(RN_m_coupler);
                    Paso_Coupler_finishCollect(RP_m_coupler);
                    /* add corrected fluxes into z_m */
                    Paso_FCTransportProblem_addCorrectedFluxes(z_m,flux_matrix,RN_m_coupler,RP_m_coupler); 
                    /* 
                      * now we solve the linear system to get the correction dt2:
                      *
                     */
                    if (fctp->theta > 0) {
                            /*  set the right hand side of the linear system: */
                            options->tolerance=inner_tol;
                            options->verbose=TRUE;
                            Paso_solve(fctp->iteration_matrix,du_m,z_m,options);
                            /* TODO: check errors ! */
                      } else {
                            #pragma omp parallel for private(i,rtmp,rtmp1)
                            for (i = 0; i < n_rows; ++i) {
                                rtmp=fctp->lumped_mass_matrix[i];
                                if (ABS(rtmp)>0.) {
                                   rtmp1=z_m[i]/rtmp;
                                } else {
                                   rtmp1=0;
                                }
                                du_m[i]=rtmp1;
                            }
                      }
                      /* 
                       * update u and calculate norm of du_m and the new u:
                       * 
                       */
                       norm_u=0.;
                       norm_du=0.;
                       #pragma omp parallel private(local_norm_u,local_norm_du)
                       {
                           local_norm_u=0.;
                           local_norm_du=0.;
                           #pragma omp for schedule(static) private(i)
                           for (i=0;i<n_rows;++i) {
                                u_m[i]+=omega*du_m[i];
                                local_norm_u=MAX(local_norm_u,ABS(u_m[i]));
                                local_norm_du=MAX(local_norm_du,ABS(du_m[i]));
                           }
                           #pragma omp critical 
                           {
                               norm_u=MAX(norm_u,local_norm_u);
                               norm_du=MAX(norm_du,local_norm_du);
                           }
                       }
                       Paso_Coupler_startCollect(u_m_coupler,u_m);
                       #ifdef PASO_MPI
                          local_norm[0]=norm_u;
                          local_norm[1]=norm_du;
	                  MPI_Allreduce(local_norm,norm, 2, MPI_DOUBLE, MPI_MAX, fctp->mpi_info->comm);
                          norm_u=norm[0];
                          norm_du=norm[1];
                       #endif
                       norm_du*=omega;
                       converged=(norm_du <= tolerance * norm_u);
                       m++;
                       printf("iteration step %d: norm u and du_m : %e %e\n",m,norm_u,norm_du);
                       /* TODO: check if du_m has been redu_mced */
                       Paso_Coupler_finishCollect(u_m_coupler);
                    } /* end of inner iteration */
                    SWAP(fctp->u,u_m,double*);
                    SWAP(fctp->u_coupler,u_m_coupler,Paso_Coupler*);
                    n++;
               } /* end of time integration */
               options->tolerance=tolerance;
               #pragma omp parallel for schedule(static) private(i)
               for (i=0;i<n_rows;++i) u[i]=fctp->u[i];
        }
        /* 
         *  clean-up:
         *
         */
        TMPMEMFREE(b_n);
        Paso_SystemMatrix_free(flux_matrix);
        TMPMEMFREE(sourceP);
        TMPMEMFREE(sourceN);
        TMPMEMFREE(uTilde_n);
        TMPMEMFREE(QN_n);
        TMPMEMFREE(QP_n);
        TMPMEMFREE(RN_m);
        TMPMEMFREE(RP_m);
        TMPMEMFREE(z_m);
        TMPMEMFREE(du_m);
        TMPMEMFREE(u_m);
        Paso_Coupler_free(QN_n_coupler);
        Paso_Coupler_free(QP_n_coupler);
        Paso_Coupler_free(RN_m_coupler);
        Paso_Coupler_free(RP_m_coupler);
        Paso_Coupler_free(uTilde_n_coupler);
        Paso_Coupler_free(u_m_coupler);
}
