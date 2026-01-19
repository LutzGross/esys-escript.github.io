
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/*
   Crude modifications and translations for Paso by Matt Davies and Lutz Gross
*/

#include "Solver.h"
#include "SystemMatrix.h"

namespace paso {

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with
*  preconditioning.
*
*  Convergence test: norm( b - A*x )< TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  A       (input)
*
*  R       (input) DOUBLE PRECISION array, dimension N.
*          On entry, residual of initial guess X
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess.
*
*  ITER    (input/output) INT
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x )
*          On output, the final value of this measure.
*
*  return value
*
*          = SOLVER_NO_ERROR: Successful exit. Iterated approximate solution returned.
*          = SOLVER_MAXITER_REACHED
*          = SOLVER_INPUT_ERROR Illegal parameter:
*          = SOLVER_BREAKDOWN: If parameters RHO or OMEGA become smaller
*          = SOLVER_MEMORY_ERROR : If parameters RHO or OMEGA become smaller
*
*  ==============================================================
*/

SolverResult Solver_BiCGStab(SystemMatrix_ptr<double> A, double* r, double* x,
                             dim_t* iter, double* tolerance, Performance* pp)
{
  /* Local variables */
  double *rtld=NULL,*p=NULL,*v=NULL,*t=NULL,*phat=NULL,*shat=NULL,*s=NULL;/*, *buf1=NULL, *buf0=NULL;*/
  double beta,norm_of_residual=0,sum_1,sum_2,sum_3,sum_4,norm_of_residual_global=0;
  double alpha=0, omega=0, omegaNumtr, omegaDenumtr, rho, tol, rho1=0;
#ifdef ESYS_MPI
  double loc_sum[2], sum[2];
#endif
  dim_t num_iter=0,maxit,num_iter_global=0;
  dim_t i0;
  bool breakFlag=false, maxIterFlag=false, convergeFlag=false;
  SolverResult status = NoError;
  double *resid = tolerance;
  dim_t n = A->getTotalNumRows();

  /* Test the input parameters. */

  if (n < 0) {
    status = InputError;
  } else {
    /* allocate memory: */
    rtld=new double[n];
    p=new double[n];
    v=new double[n];
    t=new double[n];
    phat=new double[n];
    shat=new double[n];
    s=new double[n];

    /* now bicgstab starts : */
    maxit = *iter;
    tol = *resid;

    num_iter =0;
    convergeFlag=false;
    maxIterFlag=false;
    breakFlag=false;

    /* initialise arrays */

    #pragma omp parallel for private(i0) schedule(static)
    for (i0 = 0; i0 < n; i0++) {
        rtld[i0]=0;
        p[i0]=0;
        v[i0]=0;
        t[i0]=0;
        phat[i0]=0;
        shat[i0]=0;
        rtld[i0] = r[i0];
    }

    /*     Perform BiConjugate Gradient Stabilized iteration. */

    L10:
      ++(num_iter);
        sum_1 = 0;
        sum_2 = 0;
        sum_3 = 0;
        sum_4 = 0;
        omegaNumtr = 0.0;
      omegaDenumtr = 0.0;
      #pragma omp parallel for private(i0) reduction(+:sum_1) schedule(static)
      for (i0 = 0; i0 < n; i0++) sum_1 += rtld[i0] * r[i0];
      #ifdef ESYS_MPI
          loc_sum[0] = sum_1;
          MPI_Allreduce(loc_sum, &sum_1, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
      #endif
      rho = sum_1;

      if (! (breakFlag = (std::abs(rho) <= TOLERANCE_FOR_SCALARS))) {
        /*        Compute vector P. */

        if (num_iter > 1) {
          beta = rho / rho1 * (alpha / omega);
          #pragma omp parallel for private(i0) schedule(static)
          for (i0 = 0; i0 < n; i0++) p[i0] = r[i0] + beta * (p[i0] - omega * v[i0]);
        } else {
          #pragma omp parallel for private(i0) schedule(static)
          for (i0 = 0; i0 < n; i0++) p[i0] = r[i0];
        }

        /*        Compute direction adjusting vector PHAT and scalar ALPHA. */

        A->solvePreconditioner(&phat[0], &p[0]);
        A->MatrixVector_CSR_OFFSET0(PASO_ONE, &phat[0], PASO_ZERO, &v[0]);

        #pragma omp parallel for private(i0) reduction(+:sum_2) schedule(static)
        for (i0 = 0; i0 < n; i0++) sum_2 += rtld[i0] * v[i0];
        #ifdef ESYS_MPI
           loc_sum[0] = sum_2;
            MPI_Allreduce(loc_sum, &sum_2, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
        #endif
        if (! (breakFlag = (std::abs(sum_2) <= TOLERANCE_FOR_SCALARS))) {
           alpha = rho / sum_2;

           #pragma omp parallel for private(i0) reduction(+:sum_3) schedule(static)
           for (i0 = 0; i0 < n; i0++) {
             r[i0] -= alpha * v[i0];
             s[i0] = r[i0];
             sum_3 += s[i0] * s[i0];
           }
           #ifdef ESYS_MPI
               loc_sum[0] = sum_3;
               MPI_Allreduce(loc_sum, &sum_3, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
           #endif
           norm_of_residual = sqrt(sum_3);

           /*        Early check for tolerance. */
           if ( (convergeFlag = (norm_of_residual <= tol)) ) {
             #pragma omp parallel for  private(i0) schedule(static)
             for (i0 = 0; i0 < n; i0++) x[i0] += alpha * phat[i0];
             maxIterFlag = false;
             breakFlag = false;
           } else {
             /*           Compute stabilizer vector SHAT and scalar OMEGA. */
             A->solvePreconditioner(&shat[0], &s[0]);
             A->MatrixVector_CSR_OFFSET0(PASO_ONE, &shat[0],PASO_ZERO,&t[0]);

             #pragma omp parallel for private(i0) reduction(+:omegaNumtr,omegaDenumtr) schedule(static)
             for (i0 = 0; i0 < n; i0++) {
               omegaNumtr +=t[i0] * s[i0];
               omegaDenumtr += t[i0] * t[i0];
             }
             #ifdef ESYS_MPI
                loc_sum[0] = omegaNumtr;
                loc_sum[1] = omegaDenumtr;
                MPI_Allreduce(loc_sum, sum, 2, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                omegaNumtr=sum[0];
                omegaDenumtr=sum[1];
             #endif
             if (! (breakFlag = (std::abs(omegaDenumtr) <= TOLERANCE_FOR_SCALARS))) {
                omega = omegaNumtr / omegaDenumtr;

                #pragma omp parallel for private(i0) reduction(+:sum_4) schedule(static)
                for (i0 = 0; i0 < n; i0++) {
                  x[i0] += alpha * phat[i0] + omega * shat[i0];
                  r[i0] = s[i0]-omega * t[i0];
                  sum_4 += r[i0] * r[i0];
                }
                #ifdef ESYS_MPI
                   loc_sum[0] = sum_4;
                    MPI_Allreduce(loc_sum, &sum_4, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                #endif
                norm_of_residual = sqrt(sum_4);
                convergeFlag = norm_of_residual <= tol;
                maxIterFlag = num_iter > maxit;
                breakFlag = (std::abs(omega) <= TOLERANCE_FOR_SCALARS);
              }
           }
        }
        if (!(convergeFlag || maxIterFlag || breakFlag)) {
          rho1 = rho;
          goto L10;
        }
      }
      /* end of iteration */
      num_iter_global=num_iter;
      norm_of_residual_global=norm_of_residual;
      if (maxIterFlag) {
            status = MaxIterReached;
      } else if (breakFlag) {
            status = Breakdown;
      }
    }
    delete[] rtld;
    delete[] p;
    delete[] v;
    delete[] t;
    delete[] phat;
    delete[] shat;
    delete[] s;
    *iter=num_iter_global;
    *resid=norm_of_residual_global;

    return status;
}

} // namespace paso

