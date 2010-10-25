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

/* Paso: SystemMatrix: controls iterative linear system solvers  */ 

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"
#include "esysUtils/blocktimer.h"

/***********************************************************************************/

/*  free space */

void Paso_Solver_free(Paso_SystemMatrix* A) {
   Paso_SystemMatrix_freePreconditioner(A);
}
/*  call the iterative solver: */

void Paso_Solver(Paso_SystemMatrix* A,double* x,double* b,
                 Paso_Options* options,Paso_Performance* pp) {

   double norm2_of_b,tol,tolerance,time_iter,net_time_start;
   double *r=NULL,norm2_of_residual,last_norm2_of_residual,norm_max_of_b;
   double norm2_of_b_local,norm_max_of_b_local,norm2_of_residual_local;
   double norm_max_of_residual_local,norm_max_of_residual;
   double last_norm_max_of_residual,*scaling;
#ifdef ESYS_MPI
   double loc_norm;
#endif
   dim_t i,totIter=0,cntIter,method;
   bool_t finalizeIteration;
   err_t errorCode=SOLVER_NO_ERROR;
   dim_t numSol = Paso_SystemMatrix_getTotalNumCols(A);
   dim_t numEqua = Paso_SystemMatrix_getTotalNumRows(A);
   double blocktimer_precond, blocktimer_start = blocktimer_time();
   double *x0=NULL;

 
     Esys_resetError();
     tolerance=options->tolerance;
     if (tolerance < 100.* EPSILON) {
       Esys_setError(VALUE_ERROR,"Paso_Solver: Tolerance is too small.");
     }
     if (tolerance >1.) {
       Esys_setError(VALUE_ERROR,"Paso_Solver: Tolerance mut be less than one.");
     }
     method=Paso_Options_getSolver(options->method,PASO_PASO,options->symmetric,A->mpi_info);
     /* check matrix type */
     if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1) ) {
       Esys_setError(TYPE_ERROR,"Paso_Solver: Iterative solver requires CSR format with unsymmetric storage scheme and index offset 0.");
     }
     if (A->col_block_size != A->row_block_size) {
        Esys_setError(TYPE_ERROR,"Paso_Solver: Iterative solver requires row and column block sizes to be equal.");
     }
     if (Paso_SystemMatrix_getGlobalNumCols(A) != Paso_SystemMatrix_getGlobalNumRows(A)) {
        Esys_setError(TYPE_ERROR,"Paso_Solver: Iterative solver requires a square matrix.");
        return;
     }
     /*if (A->block_size != 1 && options->preconditioner==PASO_AMG) {
        Esys_setError(TYPE_ERROR,"Paso_Solver: AMG on systems not supported yet.");
     }
     */
     time_iter=Esys_timer();
     /* this for testing only */
     if (method==PASO_NONLINEAR_GMRES) {
        Paso_Function* F=NULL;
        F=Paso_Function_LinearSystem_alloc(A,b,options);
        Paso_SystemMatrix_solvePreconditioner(A,x,b);
        errorCode=Paso_Solver_NewtonGMRES(F,x,options,pp);
        if (errorCode!=NO_ERROR) {
           Esys_setError(SYSTEM_ERROR,"Paso_Solver_NewtonGMRES: an error has occured.");
        }
        Paso_Function_LinearSystem_free(F);
        return;
     }
     
     /* ========================= */
     Performance_startMonitor(pp,PERFORMANCE_ALL);
     if (Esys_noError()) {
        /* get normalization */
        scaling=Paso_SystemMatrix_borrowNormalization(A);
        if (Esys_noError()) {
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
           #ifdef ESYS_MPI
           {
               loc_norm = norm2_of_b;
               MPI_Allreduce(&loc_norm,&norm2_of_b, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
               loc_norm = norm_max_of_b;
               MPI_Allreduce(&loc_norm,&norm_max_of_b, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
           }
           #endif
           norm2_of_b=sqrt(norm2_of_b);
         /* if norm2_of_b==0 we are ready: x=0 */
         if ( IS_NAN(norm2_of_b) || IS_NAN(norm_max_of_b) ) {
            Esys_setError(VALUE_ERROR,"Paso_Solver: Matrix or right hand side contains undefined values.");
         } else if (norm2_of_b <=0.) {

            #pragma omp parallel for private(i) schedule(static)
            for (i = 0; i < numSol; i++) x[i]=0.;
            if (options->verbose) printf("right hand side is identical zero.\n");

         } else {

            if (options->verbose) {
               printf("Paso_Solver: l2/lmax-norm of right hand side is  %e/%e.\n",norm2_of_b,norm_max_of_b);
               printf("Paso_Solver: l2/lmax-stopping criterion is  %e/%e.\n",norm2_of_b*tolerance,norm_max_of_b*tolerance);
               switch (method) {
               case PASO_BICGSTAB:
                  printf("Paso_Solver: Iterative method is BiCGStab.\n");
                  break;
               case PASO_PCG:
                  printf("Paso_Solver: Iterative method is PCG.\n");
                  break;
               case PASO_TFQMR:
                  printf("Paso_Solver: Iterative method is TFQMR.\n");
                  break;
               case PASO_MINRES:
                  printf("Paso_Solver: Iterative method is MINRES.\n");
                  break;
               case PASO_PRES20:
                  printf("Paso_Solver: Iterative method is PRES20.\n");
                  break;
               case PASO_GMRES:
                  if (options->restart>0) {
                     printf("Paso_Solver: Iterative method is GMRES(%d,%d)\n",options->truncation,options->restart);
                  } else {
                     printf("Paso_Solver: Iterative method is GMRES(%d)\n",options->truncation);
                  }
                  break;
               }

            }

            /* construct the preconditioner */
            blocktimer_precond = blocktimer_time();
            Performance_startMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
	    Paso_SystemMatrix_setPreconditioner(A,options);
            Performance_stopMonitor(pp,PERFORMANCE_PRECONDITIONER_INIT);
            blocktimer_increment("Paso_Solver_setPreconditioner()", blocktimer_precond);
            options->set_up_time=Esys_timer()-time_iter;
            if (! Esys_noError()) return;


              /* get an initial guess by evaluating the preconditioner */
	      Paso_SystemMatrix_solvePreconditioner(A,x,b);

              /* start the iteration process :*/
              r=TMPMEMALLOC(numEqua,double);
              x0=TMPMEMALLOC(numEqua,double);
              Esys_checkPtr(r);
	      Esys_checkPtr(x0);
              if (Esys_noError()) {

                 totIter = 1;
                 finalizeIteration = FALSE;
                 last_norm2_of_residual=norm2_of_b;
                 last_norm_max_of_residual=norm_max_of_b;
		 net_time_start=Esys_timer();

                 /* Loop */
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
                    #ifdef ESYS_MPI
                        loc_norm = norm2_of_residual;
                        MPI_Allreduce(&loc_norm,&norm2_of_residual, 1, MPI_DOUBLE, MPI_SUM, A->mpi_info->comm);
                        loc_norm = norm_max_of_residual;
                        MPI_Allreduce(&loc_norm,&norm_max_of_residual, 1, MPI_DOUBLE, MPI_MAX, A->mpi_info->comm);
                    #endif
                    norm2_of_residual =sqrt(norm2_of_residual);
                    options->residual_norm=norm2_of_residual;

                  if (options->verbose) printf("Paso_Solver: Step %5d: l2/lmax-norm of residual is  %e/%e",totIter,norm2_of_residual,norm_max_of_residual);

                  if (totIter>1 && norm2_of_residual>=last_norm2_of_residual &&  norm_max_of_residual>=last_norm_max_of_residual) {

                     if (options->verbose) printf(" divergence!\n");
                     Esys_setError(DIVERGED, "Paso_Solver: No improvement during iteration. Iterative solver gives up.");

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
                           tol=tolerance*norm2_of_residual/norm2_of_b; 
                           errorCode = Paso_Solver_TFQMR(A, r, x0, &cntIter, &tol, pp);
                           #pragma omp for private(i) schedule(static)
                           for (i = 0; i < numEqua; i++) {
                            x[i]+= x0[i];
                           }
                           break;
                        case PASO_MINRES:
                           /* tol=tolerance*norm2_of_residual/norm2_of_b; */
                           errorCode = Paso_Solver_MINRES(A, r, x, &cntIter, &tol, pp);
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
                        if (errorCode==SOLVER_NO_ERROR) {
                           finalizeIteration = FALSE;
                        } else if (errorCode==SOLVER_MAXITER_REACHED) {
                           Esys_setError(DIVERGED,"Paso_Solver: maximum number of iteration step reached.\nReturned solution does not fulfil stopping criterion.");
                           if (options->verbose) printf("Paso_Solver: Maximum number of iterations reached.!\n");
                        } else if (errorCode == SOLVER_INPUT_ERROR ) {
                           Esys_setError(SYSTEM_ERROR,"Paso_Solver: illegal dimension in iterative solver.");
                           if (options->verbose) printf("Paso_Solver: Internal error!\n");
			} else if (errorCode == SOLVER_NEGATIVE_NORM_ERROR) {
			   Esys_setError(VALUE_ERROR,"Paso_Solver: negative energy norm (try other solver or preconditioner).");
			   if (options->verbose) printf("Paso_Solver: negative energy norm (try other solver or preconditioner)!\n");
                        } else if ( errorCode == SOLVER_BREAKDOWN ) {
                           if (cntIter <= 1) {
                              Esys_setError(ZERO_DIVISION_ERROR, "Paso_Solver: fatal break down in iterative solver.");
                              if (options->verbose) printf("Paso_Solver: Uncurable break down!\n");
                           } else {
                              if (options->verbose) printf("Paso_Solver: Breakdown at iter %d (residual = %e). Restarting ...\n", totIter, tol);
                              finalizeIteration = FALSE;
                           }
                        } else {
			   Esys_setError(SYSTEM_ERROR,"Paso_Solver:generic error in solver.");
			   if (options->verbose) printf("Paso_Solver: generic error in solver!\n");
                        }
                      } else {
                         if (options->verbose) printf(". convergence! \n");
                         options->converged=TRUE;
                      }
                   }
                 } /* while */
                 options->net_time=Esys_timer()-net_time_start;
              }
              MEMFREE(r);
              MEMFREE(x0);
              options->num_level=0;
              options->num_inner_iter=0;
              options->num_iter=totIter;
           }
        }
      }
   options->time=Esys_timer()-time_iter;
   Performance_stopMonitor(pp,PERFORMANCE_ALL);
   blocktimer_increment("Paso_Solver()", blocktimer_start);

}
