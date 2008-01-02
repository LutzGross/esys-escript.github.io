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

              M*u= M*u_last + dt*F(u_last)
 
        with F(u) =  L u + f_a(u). 


        For the case theta>0 we set 

            A=L-M/(theta*dt) (stored into transport_matrix) 
            F(u)=L u + f_a(u) = M/(theta*dt) u + A u + f_a(u)
            b(u)= - M/(theta*dt) * u + F(u) =  A u + f_a(u) 
            b_last=M/(theta*dt)*u_last + (1-theta)*dt * F(u_last)
                  =M/(theta*dt)*u_last + (1-theta)/theta * [ M/(theta*dt) u_last + A u_last + f_a(u_last) ]
                  =(1/theta**2*dt)*M u_last + (1-theta)/theta*b(u_last)
        so (*) takes the form

        b_last=(1/theta**2*dt)*M u_last + (1-theta)/theta*b(u_last)
        while (|du| > tol * |u|) :
             A*du=b_last + b(u) 
             u <- u-du      

*/

void Paso_SolverFCT_solve(Paso_FCTransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options) {

   index_t i;
   int n_substeps,n;
   double dt2=fctp->dt_max, dt2_loc, rtmp,rtmp2,t;
   printf("FFFF %e",fctp->dt_max);
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
        for (i=0;i<Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);++i) {
               fctp->transport_matrix_diagonal[i]=
                        fctp->transport_matrix->mainBlock->val[fctp->main_iptr[i]];
        }

Paso_SystemMatrix_saveMM(fctp->flux_matrix,"flux.mm");
Paso_SystemMatrix_saveMM(fctp->transport_matrix,"trans.mm");
        /*=================================================================== *
         *
         *  calculate time step size:                                           
         */
        dt2=fctp->dt_max;
        if (fctp->theta < 1.) {
            dt2=LARGE_POSITIVE_FLOAT;
            #pragma omp parallel private(dt2_loc)
            {
               dt2_loc=LARGE_POSITIVE_FLOAT;
               #pragma omp for schedule(static) private(i,rtmp,rtmp2) 
               for (i=0;i<Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);++i) {
                    rtmp=fctp->transport_matrix_diagonal[i];
                    rtmp2=fctp->lumped_mass_matrix[i];
                    if ( (rtmp<0 && rtmp2>=0.) || (rtmp>0 && rtmp2<=0.) ) {
                        dt2_loc=MIN(dt2_loc,-rtmp2/rtmp);
                    }
                }
                #pragma omp critcal 
                {
                    dt2=MIN(dt2,dt2_loc);
                }
            }
            #ifdef PASO_MPI
               dt2_loc = dt2;
	       MPI_Allreduce(&dt2_loc, &dt2, 1, MPI_DOUBLE, MPI_MIN, fctp->mpi_info->comm);
            #endif
            dt2*=1./(1.-fctp->theta);
            if (fctp->dt_max>0.) dt2=MIN(dt2,fctp->dt_max);
        }
        if (dt2 > 0.) {
         fctp->dt=dt2;
        } else {
         fctp->dt=fctp->dt_max;
        }
        if (fctp->dt < 0.) {
            Paso_setError(TYPE_ERROR,"Paso_SolverFCT_solve: dt_max must be positive.");
        } else {
           printf("minimum time step size is %e (theta = %e).\n",fctp->dt,fctp->theta);   
        }
        /* ===========================================================
         *
         *    
        */
        if (Paso_noError()) {

        }
        fctp->valid_matrices=Paso_noError();
   }

   /* b_last=M*u_last + (1-theta) * F(u_last) */

   /* decide on substepping */
   n_substeps=ceil(dt/fctp->dt);
   dt2=dt/n_substeps;
   printf("%d time steps of size is %e (theta = %e, dt_max=%e).\n",n_substeps, dt2,fctp->theta,fctp->dt);

   /* */
   Paso_SystemMatrix_allocBuffer(fctp->transport_matrix);

   if (fctp->theta>0) {
        #pragma omp parallel for schedule(static) private(i) 
        for (i=0;i<Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);++i) {
               fctp->transport_matrix_diagonal[i]=
                        fctp->transport_matrix->mainBlock->val[fctp->main_iptr[i]];
        }
   } else {
      /* u= u_last + M^-1*dt*F(u_last) */
      t=dt2;
      n=0;
      while(n<n_substeps) {
        printf("substep step %d at t=%e\n",n+1,t);
        Paso_FCTransportProblem_setFlux(fctp,fctp->u,u); /* u stores F(u_last)*/
        #pragma omp parallel for schedule(static) private(i) 
        for (i=0;i<Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);++i) {
            rtmp=fctp->u[i];
            fctp->u[i]=fctp->u[i]+dt2*u[i]/fctp->lumped_mass_matrix[i];
            printf("%d : %e %e %e\n",i,u[i],fctp->u[i],rtmp);
        }
        t+=dt2;
        n++;
      }
  }
  if (Paso_noError()) {
     #pragma omp parallel for schedule(static) private(i) 
     for (i=0;i<Paso_SystemMatrix_getTotalNumRows(fctp->flux_matrix);++i) {
       u[i]=fctp->u[i];
     }
  }
}       
