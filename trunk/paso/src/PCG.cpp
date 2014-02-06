
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/* PCG iterations */

#include "SystemMatrix.h"
#include "Paso.h"
#include "Solver.h"

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
*  PCG solves the linear system A*x = b using the
*  preconditioned conjugate gradient method plus a smoother.
*  A has to be symmetric.
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  r       (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of initial guess x.
*
*  x       (input/output) DOUBLE PRECISION array, dimension N.
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
*          = SOLVER_BREAKDOWN: If parameters RHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters RHO or OMEGA become smaller
*
*  ==============================================================
*/

/* #define PASO_DYNAMIC_SCHEDULING_MVM */

#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined __OPENMP 
#define USE_DYNAMIC_SCHEDULING
#endif

err_t Paso_Solver_PCG(
    Paso_SystemMatrix * A,
    double * r,
    double * x,
    dim_t *iter,
    double * tolerance,
    Paso_Performance* pp) {

  /* Local variables */
  dim_t num_iter=0,maxit,num_iter_global, len,rest, np, ipp;
#ifdef USE_DYNAMIC_SCHEDULING
  dim_t chunk_size=-1;
#endif
  register double ss,ss1;
  dim_t i0, istart, iend;
  bool breakFlag=FALSE, maxIterFlag=FALSE, convergeFlag=FALSE;
  err_t status = SOLVER_NO_ERROR;
  dim_t n = Paso_SystemMatrix_getTotalNumRows(A);
  double *resid = tolerance, *rs=NULL, *p=NULL, *v=NULL, *x2=NULL ;
  double tau_old,tau,beta,delta,gamma_1,gamma_2,alpha,sum_1,sum_2,sum_3,sum_4,sum_5,tol;
#ifdef ESYS_MPI
  double loc_sum[2], sum[2];
#endif
  double norm_of_residual=0,norm_of_residual_global;
  register double d;

  /* There should not be any executable code before this ifdef */

#ifdef USE_DYNAMIC_SCHEDULING

    /* Watch out for these declarations (see above) */
    char* chksz_chr;
    dim_t n_chunks;

    chksz_chr=getenv("PASO_CHUNK_SIZE_PCG");
    if (chksz_chr!=NULL) sscanf(chksz_chr, "%d",&chunk_size);
    np=omp_get_max_threads();
    chunk_size=MIN(MAX(1,chunk_size),n/np);
    n_chunks=n/chunk_size;
    if (n_chunks*chunk_size<n) n_chunks+=1;
#else
    np=omp_get_max_threads();
    len=n/np;
    rest=n-len*np;
#endif
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Start of Calculation :                                        */
/*   ---------------------                                         */
/*                                                                 */
/*                                                                 */
  rs=new double[n];
  p=new double[n];
  v=new double[n];
  x2=new double[n];

  /*     Test the input parameters. */

  if (n < 0) {
    status = SOLVER_INPUT_ERROR;
  } else if (rs==NULL || p==NULL || v==NULL || x2==NULL) {
    status = SOLVER_MEMORY_ERROR;
  } else {
    maxit = *iter;
    tol = *resid;
    Performance_startMonitor(pp,PERFORMANCE_SOLVER);
    /* initialize data */
    #pragma omp parallel private(i0, istart, iend, ipp)
    {
       #ifdef USE_DYNAMIC_SCHEDULING
           #pragma omp for schedule(dynamic, 1)
           for (ipp=0; ipp < n_chunks; ++ipp) {
              istart=chunk_size*ipp;
              iend=MIN(istart+chunk_size,n);
       #else
           #pragma omp for schedule(static)
           for (ipp=0; ipp <np; ++ipp) {
               istart=len*ipp+MIN(ipp,rest);
               iend=len*(ipp+1)+MIN(ipp+1,rest);
       #endif
               #pragma ivdep
               for (i0=istart;i0<iend;i0++) {
                 rs[i0]=r[i0];
                 x2[i0]=x[i0];
                 p[i0]=0;
                 v[i0]=0;
               } 
       #ifdef USE_DYNAMIC_SCHEDULING
         }
       #else
         }
       #endif
    }
    num_iter=0;

    /* PGH */
    /* without this we get a use of an uninitialised var below */
    tau = 0;

    /* start of iteration */
    while (!(convergeFlag || maxIterFlag || breakFlag)) {
           ++(num_iter);

           /* PGH */
           /* The next lines were commented out before I got here */
           /* v=prec(r)  */
           /* tau=v*r; */
           /* leading to the use of an uninitialised var below */

           Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
           Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER);
	   Paso_SystemMatrix_solvePreconditioner(A,v,r);
           Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER);
           Performance_startMonitor(pp,PERFORMANCE_SOLVER);

	   sum_1 = 0;
           #pragma omp parallel private(i0, istart, iend, ipp, ss)
           {
                  ss=0;
                  #ifdef USE_DYNAMIC_SCHEDULING
                      #pragma omp for schedule(dynamic, 1)
                      for (ipp=0; ipp < n_chunks; ++ipp) {
                         istart=chunk_size*ipp;
                         iend=MIN(istart+chunk_size,n);
                  #else
                      #pragma omp for schedule(static) 
                      for (ipp=0; ipp <np; ++ipp) {
                          istart=len*ipp+MIN(ipp,rest);
                          iend=len*(ipp+1)+MIN(ipp+1,rest);
                  #endif
 			  #pragma ivdep
                          for (i0=istart;i0<iend;i0++) ss+=v[i0]*r[i0];
                  #ifdef USE_DYNAMIC_SCHEDULING
                    }
                  #else
                    }
                  #endif
                  #pragma omp critical
                  {
                     sum_1+=ss;
                  }
           }
           #ifdef ESYS_MPI
	        /* In case we have many MPI processes, each of which may have several OMP threads:
	           OMP master participates in an MPI reduction to get global sum_1 */
	        loc_sum[0] = sum_1;
	        MPI_Allreduce(loc_sum, &sum_1, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
           #endif
           tau_old=tau;
           tau=sum_1;
           /* p=v+beta*p */
           #pragma omp parallel private(i0, istart, iend, ipp,beta)
           {
                  #ifdef USE_DYNAMIC_SCHEDULING
                      #pragma omp for schedule(dynamic, 1)
                      for (ipp=0; ipp < n_chunks; ++ipp) {
                         istart=chunk_size*ipp;
                         iend=MIN(istart+chunk_size,n);
                  #else
                      #pragma omp for schedule(static)
                      for (ipp=0; ipp <np; ++ipp) {
                          istart=len*ipp+MIN(ipp,rest);
                          iend=len*(ipp+1)+MIN(ipp+1,rest);
                  #endif
                          if (num_iter==1) {
 			      #pragma ivdep
                              for (i0=istart;i0<iend;i0++) p[i0]=v[i0];
                          } else {
                              beta=tau/tau_old;
 			      #pragma ivdep
                              for (i0=istart;i0<iend;i0++) p[i0]=v[i0]+beta*p[i0];
                          }
                  #ifdef USE_DYNAMIC_SCHEDULING
                    }
                  #else
                    }
                  #endif
           }
           /* v=A*p */
           Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
           Performance_startMonitor(pp,PERFORMANCE_MVM);
	   Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, p,PASO_ZERO,v);
           Performance_stopMonitor(pp,PERFORMANCE_MVM);
           Performance_startMonitor(pp,PERFORMANCE_SOLVER);

           /* delta=p*v */
	   sum_2 = 0;
           #pragma omp parallel private(i0, istart, iend, ipp,ss)
           {
                  ss=0;
                  #ifdef USE_DYNAMIC_SCHEDULING
                      #pragma omp for schedule(dynamic, 1)
                      for (ipp=0; ipp < n_chunks; ++ipp) {
                         istart=chunk_size*ipp;
                         iend=MIN(istart+chunk_size,n);
                  #else
                      #pragma omp for schedule(static)
                      for (ipp=0; ipp <np; ++ipp) {
                          istart=len*ipp+MIN(ipp,rest);
                          iend=len*(ipp+1)+MIN(ipp+1,rest);
                  #endif
 			  #pragma ivdep
                          for (i0=istart;i0<iend;i0++) ss+=v[i0]*p[i0];
                  #ifdef USE_DYNAMIC_SCHEDULING
                    }
                  #else
                    }
                  #endif
                  #pragma omp critical
                  {
                             sum_2+=ss;
                  }
           }
           #ifdef ESYS_MPI
	      loc_sum[0] = sum_2;
	      MPI_Allreduce(loc_sum, &sum_2, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
           #endif
           delta=sum_2;
           alpha=tau/delta;
   
           if (! (breakFlag = (ABS(delta) <= TOLERANCE_FOR_SCALARS))) {
               /* smoother */
	       sum_3 = 0;
	       sum_4 = 0;
               #pragma omp parallel private(i0, istart, iend, ipp,d, ss, ss1)
               {
                  ss=0;
                  ss1=0;
                  #ifdef USE_DYNAMIC_SCHEDULING
                      #pragma omp for schedule(dynamic, 1)
                      for (ipp=0; ipp < n_chunks; ++ipp) {
                         istart=chunk_size*ipp;
                         iend=MIN(istart+chunk_size,n);
                  #else
                      #pragma omp for schedule(static)
                      for (ipp=0; ipp <np; ++ipp) {
                          istart=len*ipp+MIN(ipp,rest);
                          iend=len*(ipp+1)+MIN(ipp+1,rest);
                  #endif
 			  #pragma ivdep
                          for (i0=istart;i0<iend;i0++) {
                                r[i0]-=alpha*v[i0];
                                d=r[i0]-rs[i0];
                                ss+=d*d;
                                ss1+=d*rs[i0];
                          }
                  #ifdef USE_DYNAMIC_SCHEDULING
                    }
                  #else
                    }
                  #endif
                  #pragma omp critical
                  {
                     sum_3+=ss;
                     sum_4+=ss1;
                  }
               }
               #ifdef ESYS_MPI
	           loc_sum[0] = sum_3;
	           loc_sum[1] = sum_4;
	           MPI_Allreduce(loc_sum, sum, 2, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
	           sum_3=sum[0];
	           sum_4=sum[1];
                #endif
	        sum_5 = 0;
                #pragma omp parallel private(i0, istart, iend, ipp, ss, gamma_1,gamma_2)
                {
                  gamma_1= ( (ABS(sum_3)<= PASO_ZERO) ? 0 : -sum_4/sum_3) ;
                  gamma_2= PASO_ONE-gamma_1;
                  ss=0;
                  #ifdef USE_DYNAMIC_SCHEDULING
                      #pragma omp for schedule(dynamic, 1)
                      for (ipp=0; ipp < n_chunks; ++ipp) {
                         istart=chunk_size*ipp;
                         iend=MIN(istart+chunk_size,n);
                  #else
                      #pragma omp for schedule(static)
                      for (ipp=0; ipp <np; ++ipp) {
                          istart=len*ipp+MIN(ipp,rest);
                          iend=len*(ipp+1)+MIN(ipp+1,rest);
                  #endif
 			  #pragma ivdep
                          for (i0=istart;i0<iend;i0++) {
                              rs[i0]=gamma_2*rs[i0]+gamma_1*r[i0];
                              x2[i0]+=alpha*p[i0];
                              x[i0]=gamma_2*x[i0]+gamma_1*x2[i0];
                              ss+=rs[i0]*rs[i0];
                          }
                  #ifdef USE_DYNAMIC_SCHEDULING
                    }
                  #else
                    }
                  #endif
                  #pragma omp critical
                  {
                      sum_5+=ss;
                  }
                }
                #ifdef ESYS_MPI
	           loc_sum[0] = sum_5;
	           MPI_Allreduce(loc_sum, &sum_5, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                #endif
                norm_of_residual=sqrt(sum_5);
                convergeFlag = norm_of_residual <= tol;
                maxIterFlag = num_iter > maxit;
                breakFlag = (ABS(tau) <= TOLERANCE_FOR_SCALARS);
           }
    }
    /* end of iteration */
    num_iter_global=num_iter;
    norm_of_residual_global=norm_of_residual;
    if (maxIterFlag) {
         status = SOLVER_MAXITER_REACHED;
    } else if (breakFlag) {
         status = SOLVER_BREAKDOWN;
    }
    Performance_stopMonitor(pp,PERFORMANCE_SOLVER);
    delete[] rs;
    delete[] x2;
    delete[] v;
    delete[] p;
    *iter=num_iter_global;
    *resid=norm_of_residual_global;
  }
  /*     End of PCG */
  return status;
}

