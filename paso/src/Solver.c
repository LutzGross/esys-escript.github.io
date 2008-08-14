
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: SystemMatrix: controls iterative linear system solvers  */ 

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"
#include "escript/blocktimer.h"

/***********************************************************************************/

/*  free space */

void Paso_Solver_free(Paso_SystemMatrix* A) {
   Paso_Preconditioner_free(A->solver);
   A->solver=NULL;
}
/*  call the iterative solver: */

void Paso_Solver(Paso_SystemMatrix* A,double* x,double* b,
                 Paso_Options* options,Paso_Performance* pp) {

   double norm2_of_b,tol,tolerance,time_iter;
   double *r=NULL,norm2_of_residual,last_norm2_of_residual,norm_max_of_b;
   double norm2_of_b_local,norm_max_of_b_local,norm2_of_residual_local;
   double norm_max_of_residual_local,norm_max_of_residual;
   double last_norm_max_of_residual,*scaling;
#ifdef PASO_MPI
   double loc_norm;
#endif
   dim_t i,totIter,cntIter,method;
   bool_t finalizeIteration;
   err_t errorCode;
   dim_t numSol = Paso_SystemMatrix_getTotalNumCols(A);
   dim_t numEqua = Paso_SystemMatrix_getTotalNumRows(A);
   double blocktimer_start = blocktimer_time();
 
     tolerance=MAX(options->tolerance,EPSILON);
     Paso_resetError();
     method=Paso_Options_getSolver(options->method,PASO_PASO,options->symmetric);
     /* check matrix type */
     if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1) || (A->type & MATRIX_FORMAT_SYM) ) {
       Paso_setError(TYPE_ERROR,"Iterative solver requires CSR format with unsymmetric storage scheme and index offset 0.");
     }
     if (A->col_block_size != A->row_block_size) {
        Paso_setError(TYPE_ERROR,"Iterative solver requires row and column block sizes to be equal.");
     }
     if (Paso_SystemMatrix_getGlobalNumCols(A) != Paso_SystemMatrix_getGlobalNumRows(A)) {
        Paso_setError(TYPE_ERROR,"Iterative solver requires requires a square matrix.");
        return;
     }
     /* this for testing only */
     if (method==PASO_NONLINEAR_GMRES) {
        Paso_Function* F=NULL;
        F=Paso_Function_LinearSystem_alloc(A,b,options);
        errorCode=Paso_Solver_NewtonGMRES(F,x,options,pp);
        if (errorCode==NO_ERROR) {
           Paso_setError(SYSTEM_ERROR,"Paso_Solver_NewtonGMRES: an error has occured.");
        }
        Paso_Function_LinearSystem_free(F);
     }
     /* ========================= */
     Performance_startMonitor(pp,PERFORMANCE_ALL);
     if (Paso_noError()) {
        /* get normalization */
        scaling=Paso_SystemMatrix_borrowNormalization(A);
        if (Paso_noError()) {
           /* get the norm of the right hand side */
           norm2_of_b=0.;
           norm_max_of_b=0.;
           #pragma omp parallel private(norm2_of_b_local,norm_max_of_b_local) 
           { 
                norm2_of_b_local=0.;
                norm_max_of_b_local=0.;
                #pragma omp for private(i) schedule(static)
                for (i = 0; i < numEqua ; ++i) {
                      norm2_of_b_local += b[i] * b[i];
                      norm_max_of_b_local = MAX(ABS(scaling[i]*b[i]),norm_max_of_b_local);
                }
                #pragma omp critical
                {
                   norm2_of_b += norm2_of_b_local;
                   norm_max_of_b = MAX(norm_max_of_b_local,norm_max_of_b);
                }
           }
           /* TODO: use one call */
           #ifdef PASO_MPI
           {
               loc_norm = norm2_of_b;
               MPI_Allreduce(&loc_norm,&norm2_of_b, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
               loc_norm = norm_max_of_b;
               MPI_Allreduce(&loc_norm,&norm_max_of_b, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
           }
           #endif
           norm2_of_b=sqrt(norm2_of_b);
    
         /* if norm2_of_b==0 we are ready: x=0 */
         if (norm2_of_b <=0.) {
#pragma omp parallel for private(i) schedule(static)
            for (i = 0; i < numSol; i++) x[i]=0.;
            if (options->verbose) printf("right hand side is identical zero.\n");
         } else {
            if (options->verbose) {
               printf("l2/lmax-norm of right hand side is  %e/%e.\n",norm2_of_b,norm_max_of_b);
               printf("l2/lmax-stopping criterion is  %e/%e.\n",norm2_of_b*tolerance,norm_max_of_b*tolerance);
               switch (method) {
               case PASO_BICGSTAB:
                  printf("Iterative method is BiCGStab.\n");
                  break;
               case PASO_PCG:
                  printf("Iterative method is PCG.\n");
                  break;
               case PASO_TFQMR:
                  printf("Iterative method is TFQMR.\n");
                  break;
               case PASO_PRES20:
                  printf("Iterative method is PRES20.\n");
                  break;
               case PASO_GMRES:
                  if (options->restart>0) {
                     printf("Iterative method is GMRES(%d,%d)\n",options->truncation,options->restart);
                  } else {
                     printf("Iterative method is GMRES(%d)\n",options->truncation);
                  }
                  break;
               }
            }
            /* construct the preconditioner */
              
            Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
            Paso_Solver_setPreconditioner(A,options);
            Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
            if (! Paso_noError()) return;
              time_iter=Paso_timer();
              /* get an initial guess by evaluating the preconditioner */
              Paso_Solver_solvePreconditioner(A,x,b);
              /* start the iteration process :*/
              r=TMPMEMALLOC(numEqua,double);
              Paso_checkPtr(r);
              if (Paso_noError()) {
                 totIter = 0;
                 finalizeIteration = FALSE;
                 last_norm2_of_residual=norm2_of_b;
                 last_norm_max_of_residual=norm_max_of_b;
                 while (! finalizeIteration) {
                    cntIter = options->iter_max - totIter;
                    finalizeIteration = TRUE;
                    /*     Set initial residual. */
                    norm2_of_residual = 0;
                    norm_max_of_residual = 0;
                    #pragma omp parallel for private(i) schedule(static)
                    for (i = 0; i < numEqua; i++) r[i]=b[i];
                    Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1), A, x, DBLE(1), r);
             
                    #pragma omp parallel private(norm2_of_residual_local,norm_max_of_residual_local)
                    {
                       norm2_of_residual_local = 0;
                       norm_max_of_residual_local = 0;
                       #pragma omp for private(i) schedule(static)
                       for (i = 0; i < numEqua; i++) {
                              norm2_of_residual_local+= r[i] * r[i];
                              norm_max_of_residual_local=MAX(ABS(scaling[i]*r[i]),norm_max_of_residual_local);
                       }
                       #pragma omp critical
                       {
                          norm2_of_residual += norm2_of_residual_local;
                          norm_max_of_residual = MAX(norm_max_of_residual_local,norm_max_of_residual);
                       }
                    }
                    /* TODO: use one call */
                    #ifdef PASO_MPI
                        loc_norm = norm2_of_residual;
                        MPI_Allreduce(&loc_norm,&norm2_of_residual, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                        loc_norm = norm_max_of_residual;
                        MPI_Allreduce(&loc_norm,&norm_max_of_residual, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
                    #endif
                    norm2_of_residual =sqrt(norm2_of_residual);
             
                  if (options->verbose) printf("Step %5d: l2/lmax-norm of residual is  %e/%e",totIter,norm2_of_residual,norm_max_of_residual);
                  if (totIter>0 && norm2_of_residual>=last_norm2_of_residual &&  norm_max_of_residual>=last_norm_max_of_residual) {
                     if (options->verbose) printf(" divergence!\n");
                     Paso_setError(WARNING, "No improvement during iteration. Iterative solver gives up.");
                  } else {
                     /* */
                     if (norm2_of_residual>tolerance*norm2_of_b || norm_max_of_residual>tolerance*norm_max_of_b ) {
                        tol=tolerance*MIN(norm2_of_b,0.1*norm2_of_residual/norm_max_of_residual*norm_max_of_b);
                        if (options->verbose) printf(" (new tolerance = %e).\n",tol);
                        last_norm2_of_residual=norm2_of_residual;
                        last_norm_max_of_residual=norm_max_of_residual;
                        /* call the solver */
                        switch (method) {
                        case PASO_BICGSTAB:
                           errorCode = Paso_Solver_BiCGStab(A, r, x, &cntIter, &tol, pp); 
                           break;
                        case PASO_PCG:
                           errorCode = Paso_Solver_PCG(A, r, x, &cntIter, &tol, pp); 
                           break;
                        case PASO_TFQMR:
                           errorCode = Paso_Solver_TFQMR(A, r, x, &cntIter, &tol, pp); 
                           break;
                        case PASO_PRES20:
                           errorCode = Paso_Solver_GMRES(A, r, x, &cntIter, &tol,5,20, pp); 
                           break;
                        case PASO_GMRES:
                           errorCode = Paso_Solver_GMRES(A, r, x, &cntIter, &tol,options->truncation,options->restart, pp); 
                           break;
                        }
                        totIter += cntIter;
                      
                        /* error handling  */
                        if (errorCode==NO_ERROR) {
                           finalizeIteration = FALSE;
                        } else if (errorCode==SOLVER_MAXITER_REACHED) {
                           Paso_setError(WARNING,"maximum number of iteration step reached.\nReturned solution does not fulfil stopping criterion.");
                           if (options->verbose) printf("Maximum number of iterations reached.!\n");
                        } else if (errorCode == SOLVER_INPUT_ERROR ) {
                           Paso_setError(SYSTEM_ERROR,"illegal dimension in iterative solver.");
                           if (options->verbose) printf("Internal error!\n");
                        } else if ( errorCode == SOLVER_BREAKDOWN ) {
                           if (cntIter <= 1) {
                              Paso_setError(ZERO_DIVISION_ERROR, "fatal break down in iterative solver.");
                              if (options->verbose) printf("Uncurable break down!\n");
                           } else {
                              if (options->verbose) printf("Breakdown at iter %d (residual = %e). Restarting ...\n", cntIter+totIter, tol);
                              finalizeIteration = FALSE;
                           }
                      } else {
                         if (options->verbose) printf(". convergence! \n");
                      }
                   }
                 } /* while */
              }
              MEMFREE(r);
              time_iter=Paso_timer()-time_iter;
              if (options->verbose)  {
                 printf("\ntiming: solver: %.4e sec\n",time_iter);
                 if (totIter>0) printf("timing: per iteration step: %.4e sec\n",time_iter/totIter);
              }
           }
        }
      }
   }
   Performance_stopMonitor(pp,PERFORMANCE_ALL);
   blocktimer_increment("Paso_Solver()", blocktimer_start);
}
