/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: Flux correction transport solver
 *
 * solves Mu_t=Du+Ku+q
 *
 *  where is D is diffusive (not checked)
 *  		- D is symmetric
 *  		- row sums are equal to zero.
 *  and  K is the advective part.
 *
 *        u(0) >= 0
 *
 * intially fctp->transport_matrix defines the diffusive part 
 * but the matrix is updated by the adevctive part + artificial diffusion 
 *
*/
/**************************************************************/

/* Author: l.gross@uq.edu.au */

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

/***********************************************************************************

        solve (until convergence)

        (*)  [M-theta*dt*L] du = M*u_last + (1-theta)*dt*F(u_last) - M*u + theta*dt F(u) 
             u <- u+du

        with F(u) =  L u + f_a(u)  (f_a anti-diffusion)
        and L = D + K + D_K  stored in transport_matrix
            D = Diffusive part (on input stored in transport_matrix)
            K = flux_matrix 
            D_K = artificial diffusion introduced by K
            f_a = anti-diffusion introduced by K and u
            u_last=u of last time step

        For the case theta==0 (explicit case) no iteration is required. One sets

              M*u= M*u_last + dt*b
 
        with b=F(u) =  L u + f_a(u). 


        For the case theta>0 we solve (*) is the form

            [L-M/theta*dt] du2 =M/(theta*dt)*u_last + ((1-theta)/theta)*F(u_last) - M/(theta*dt)*u + F(u) 
 
        for du2=-du which is solved as

            A du2  = r= b + f(u)
            
        with 

            du=-du2
            A=L-M/(theta*dt) (stored into transport_matrix) 
            f(u) =-M/(theta*dt)*u + F(u) = - M/(theta*dt)*u + Lu + f_a(u) = f_a(u) + Au
            F(u)=f(u)+M/(theta*dt)*u
            b= M/(theta*dt)*u_last+(1-theta)/theta*F(u_last)=
             = M/(theta*dt)*u_last+(1-theta)/theta*(f(u)+M/(theta*dt)*u_last)
             = M/(theta*dt)*u_last+(1-theta)/theta*(f(u_last)+M/(theta*dt)*u_last)
             = M*1./(theta*dt)*(1+(1-theta)/theta)*u_last+(1-theta)/theta*f(u_last)

*/

void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options) {

   index_t i, j;
   int n_substeps,n, iter;
   double fac, fac2, *b=NULL, *f=NULL, *du=NULL;
   double dt2=fctp->dt_max, dt2_loc, rtmp,rtmp2,t;
   double local_norm[2],norm[2],local_norm_u,local_norm_du,norm_u,norm_du, tolerance;
   dim_t n_rows=Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);
   bool_t converged;
   if (dt<=0.) {
       Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
   }
   if (! fctp->valid_matrices) {

        /* extract the row sum of the advective part */
        Paso_SystemMatrix_rowSum(fctp->flux_matrix,fctp->row_sum_flux_matrix);
        /* add the advective part + artificial diffusion to the diffusive part = transport-matrix */
        Paso_FCTransportProblem_addAdvectivePart(fctp,1.);
        /* create a copy of the main diagonal entires of the transport-matrix */
        #pragma omp parallel for schedule(static) private(i) 
        for (i=0;i<n_rows;++i) {
               fctp->transport_matrix_diagonal[i]=
                        fctp->transport_matrix->mainBlock->val[fctp->main_iptr[i]];
        }

 /* Paso_SystemMatrix_saveMM(fctp->flux_matrix,"flux.mm"); 
 Paso_SystemMatrix_saveMM(fctp->transport_matrix,"trans.mm");  */
        /*=================================================================== *
         *
         *  calculate time step size:                                           
         */
        dt2=LARGE_POSITIVE_FLOAT;
        if (fctp->theta < 1.) {
            #pragma omp parallel private(dt2_loc)
            {
               dt2_loc=LARGE_POSITIVE_FLOAT;
               #pragma omp for schedule(static) private(i,rtmp,rtmp2) 
               for (i=0;i<n_rows;++i) {
                    rtmp=fctp->transport_matrix_diagonal[i];
                    rtmp2=fctp->lumped_mass_matrix[i];
                    if ( (rtmp<0 && rtmp2>=0.) || (rtmp>0 && rtmp2<=0.) ) {
                        dt2_loc=MIN(dt2_loc,-rtmp2/rtmp);
                    }
                }
                #pragma omp critical 
                {
                    dt2=MIN(dt2,dt2_loc);
                }
            }
            #ifdef PASO_MPI
               dt2_loc = dt2;
	       MPI_Allreduce(&dt2_loc, &dt2, 1, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
            #endif
            if (dt2<LARGE_POSITIVE_FLOAT) dt2*=1./(1.-fctp->theta);
        }
        if (dt2 <= 0.) {
            Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt must be positive.");
        } else {
           if (dt2<LARGE_POSITIVE_FLOAT) printf("minimum time step size is %e (theta = %e).\n",dt2,fctp->theta);   
        }
        fctp->dt=dt2;
        fctp->dt_max=dt2; /* FIXME: remove*/
        fctp->valid_matrices=Paso_noError();
   }

   /* */
   Paso_SystemMatrix_allocBuffer(fctp->transport_matrix);
   /* 
    * allocate memory:
    *
    */
   b=MEMALLOC(n_rows,double);
   Paso_checkPtr(b);
   if (fctp->theta>0) {
       b=MEMALLOC(n_rows,double);
       du=MEMALLOC(n_rows,double);
       f=MEMALLOC(n_rows,double);
       Paso_checkPtr(du);
       Paso_checkPtr(f);
   }
   if (Paso_noError()) {
       /*
        *    Preparation:
        *
        */
    
       /* decide on substepping */
       if (fctp->dt < LARGE_POSITIVE_FLOAT) {
          n_substeps=ceil(dt/fctp->dt);
       } else {
          n_substeps=1.;
       }
       dt2=dt/n_substeps;
       printf("%d time steps of size is %e (theta = %e, dt_max=%e).\n",n_substeps, dt2,fctp->theta,fctp->dt);
        
       /*
        *  implicit case:
        *
        *   A=L-M/(theta*dt) (stored into transport_matrix) 
        *
        * b= M/(theta*dt)*u_last+(1-theta)/theta)*F(u_last)  
        *
        */
       if (fctp->theta>0) {
          Paso_solve_free(fctp->transport_matrix);
          fac=1./(fctp->theta*dt2);
          fac2=(1.-fctp->theta)/fctp->theta;
          #pragma omp parallel for schedule(static) private(i)
          for (i=0;i<n_rows;++i) {
                fctp->transport_matrix->mainBlock->val[fctp->main_iptr[i]]=
                                            fctp->transport_matrix_diagonal[i]-fctp->lumped_mass_matrix[i]*fac;
          }
       }
       #pragma omp parallel for schedule(static) private(i)
       for (i=0;i<n_rows;++i) {
         u[i]=fctp->u[i];
/* printf("A %d : %e %e\n",i,u[i],fctp->u[i]); */
       }
       /*
        * now the show can begin:
        *
        */
        t=dt2;
        n=0;
        tolerance=options->tolerance;
        while(n<n_substeps && Paso_noError()) {
            printf("substep step %d at t=%e\n",n+1,t);
            if (fctp->theta>0.) {
                /*
                 * implicit scheme:
                 *
                 */
                if (fctp->theta<1.) {
                   Paso_FCTransportProblem_setFlux(fctp,u,b);  /* b=f(u_last) */
                   #pragma omp parallel for schedule(static) private(i)
                   for (i=0;i<n_rows;++i) {
                        b[i]=fctp->lumped_mass_matrix[i]*fac*(1.+fac2)*u[i]+fac2*b[i];
                    }
                } else {
                   #pragma omp parallel for schedule(static) private(i)
                   for (i=0;i<n_rows;++i) {
                       b[i]=fctp->lumped_mass_matrix[i]*fac*u[i];
                   }
                }
                /*
                 * Enter iteration on a time step :
                 *
                 */
                 iter=0;
                 converged=FALSE;
                 while ( (!converged) && (iter<50) && Paso_noError()) {
                    printf("iteration step %d\n",iter+1);
                    Paso_FCTransportProblem_setFlux(fctp,u,f); 
                    #pragma omp parallel for schedule(static) private(i)
                    for (i=0;i<n_rows;++i) {
                       f[i]+=b[i];
                    }
                    options->tolerance=1.e-3;
                    Paso_solve(fctp->transport_matrix,du,f,options);
                    /* update u and calculate norms */
                    norm_u=0.;
                    norm_du=0.;
                    #pragma omp parallel private(local_norm_u,local_norm_du)
                    {
                        local_norm_u=0.;
                        local_norm_du=0.;
                        #pragma omp for schedule(static) private(i)
                        for (i=0;i<n_rows;++i) {
                             u[i]-=du[i];
                             local_norm_u=MAX(local_norm_u,ABS(u[i]));
                             local_norm_du=MAX(local_norm_du,ABS(du[i]));
                        }
                        #pragma omp critical 
                        {
                            norm_u=MAX(norm_u,local_norm_u);
                            norm_du=MAX(norm_du,local_norm_du);
                        }
                    }
                    #ifdef PASO_MPI
                       local_norm[0]=norm_u;
                       local_norm[1]=norm_du;
	               MPI_Allreduce(local_norm,norm, 2, MPI_DOUBLE, MPI_MAX, fctp->mpi_info->comm);
                       norm_u=norm[0];
                       norm_du=norm[1];
                    #endif
                    converged=(norm_du <= tolerance * norm_u);
                    iter++;
                    printf("iteration step %d: norm u and du : %e %e\n",iter,norm_u,norm_du);
                 }
             } else {
                /*
                  * explicit scheme:
                  *
                 */
                #pragma omp parallel for schedule(static) private(i) 
                for (i=0;i<n_rows;++i) {
                    u[i]+=dt2*b[i]/fctp->lumped_mass_matrix[i];
                }
             }
             /* and the next time step */
             t+=dt2;
             n++;
        }
        /*
         * save last u 
         *
         */
        if (Paso_noError()) {
             #pragma omp parallel for schedule(static) private(i)
             for (i=0;i<n_rows;++i) {
                 fctp->u[i]=u[i];
             }
        }
    }
    /*
     *
     * clean-up
     *
     */
    if (fctp->theta>0) {
            MEMFREE(b);
            MEMFREE(du);
            MEMFREE(f);
    }
}
