/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SystemMatrix: controls iterative linear system solvers  */ 

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"

/***********************************************************************************/

/*  free space */

void Paso_Solver_free(Paso_SystemMatrix* A) {
    Paso_Preconditioner_free(A->solver);
    A->solver=NULL;
}
/*  call the iterative solver: */

void Paso_Solver(Paso_SystemMatrix* A,double* x,double* b,Paso_Options* options,Paso_Performance* pp) {
    double norm2_of_b,tol,tolerance,time_iter,*r=NULL,norm2_of_residual,last_norm2_of_residual,norm_max_of_b,
           norm2_of_b_local,norm_max_of_b_local,norm2_of_residual_local,norm_max_of_residual_local,norm_max_of_residual,
           last_norm_max_of_residual,*scaling;
     dim_t i,totIter,cntIter,method;
     bool_t finalizeIteration;
     err_t errorCode;
     dim_t n_col = A->num_cols * A-> col_block_size;
     dim_t n_row = A->num_rows * A-> row_block_size;
 
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
     if (A->num_cols != A->num_rows) {
        Paso_setError(TYPE_ERROR,"Iterative solver requires requires a square matrix.");
        return;
     }
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
                for (i = 0; i < n_row ; ++i) {
                      norm2_of_b_local += b[i] * b[i];
                      norm_max_of_b_local = MAX(ABS(scaling[i]*b[i]),norm_max_of_b_local);
                }
                #pragma omp critical
                {
                   norm2_of_b += norm2_of_b_local;
                   norm_max_of_b = MAX(norm_max_of_b_local,norm_max_of_b);
                }
           }
           norm2_of_b=sqrt(norm2_of_b);
    
           /* if norm2_of_b==0 we are ready: x=0 */
           if (norm2_of_b <=0.) {
               #pragma omp parallel for private(i) schedule(static)
               for (i = 0; i < n_col; i++) x[i]=0.;
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
          
              /* get an initial guess by evaluating the preconditioner */
              time_iter=Paso_timer();
              #pragma omp parallel
              Paso_Solver_solvePreconditioner(A,x,b);
              if (! Paso_noError()) return;
              /* start the iteration process :*/
              r=TMPMEMALLOC(n_row,double);
              if (Paso_checkPtr(r)) return;
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
                 #pragma omp parallel private(norm2_of_residual_local,norm_max_of_residual_local)
                 {
                    #pragma omp for private(i) schedule(static)
                    for (i = 0; i < n_row; i++) r[i]=b[i];
          
                    Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1), A, x, DBLE(1), r);
          
                    norm2_of_residual_local = 0;
                    norm_max_of_residual_local = 0;
                    #pragma omp for private(i) schedule(static)
                    for (i = 0; i < n_row; i++) {
                           norm2_of_residual_local+= r[i] * r[i];
                           norm_max_of_residual_local=MAX(ABS(scaling[i]*r[i]),norm_max_of_residual_local);
                    }
                    #pragma omp critical
                    {
                       norm2_of_residual += norm2_of_residual_local;
                       norm_max_of_residual = MAX(norm_max_of_residual_local,norm_max_of_residual);
                    }
                 }
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
                                  errorCode = Paso_Solver_BiCGStab(A, r, x, &cntIter, &tol,pp); 
                                  break;
                              case PASO_PCG:
                                  errorCode = Paso_Solver_PCG(A, r, x, &cntIter, &tol,pp); 
                                  break;
                              case PASO_PRES20:
                                  errorCode = Paso_Solver_GMRES(A, r, x, &cntIter, &tol,5,20,pp); 
                                  break;
                              case PASO_GMRES:
                                  errorCode = Paso_Solver_GMRES(A, r, x, &cntIter, &tol,options->truncation,options->restart,pp); 
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
                        } else if (errorCode == SOLVER_MEMORY_ERROR) {
            	          Paso_setError(MEMORY_ERROR,"memory allocation failed.");
                             if (options->verbose) printf("Memory allocation failed!\n");
                        } else if (errorCode !=SOLVER_NO_ERROR ) {
            	          Paso_setError(SYSTEM_ERROR,"unidentified error in iterative solver.");
                             if (options->verbose) printf("Unidentified error!\n");
                        }
                   } else {
                      if (options->verbose) printf(". convergence! \n");
                   }
                }
              } /* while */
              MEMFREE(r);
              time_iter=Paso_timer()-time_iter;
              if (options->verbose)  {
                 printf("timing: solver: %.4e sec\n",time_iter);
                 if (totIter>0) printf("timing: per iteration step: %.4e sec\n",time_iter/totIter);
              }
           }
        }
      }
      Performance_stopMonitor(pp,PERFORMANCE_ALL);
}
