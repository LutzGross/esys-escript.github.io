/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix: controls iterative linear system solvers  */ 

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Solver.h"
#include "Common.h"
#if ITERATIVE_SOLVER == NO_LIB
#if PTR_OFFSET !=0 || INDEX_OFFSET!=0
#error Finley library usage requires PTR_OFFSET=0 and INDEX_OFFSET=0
#endif
#endif


/***********************************************************************************/

/*  free space */

void Finley_Solver_free(Finley_SystemMatrix* A) {
#if ITERATIVE_SOLVER == NO_LIB
    Finley_Preconditioner_free(A->iterative);
    A->iterative=NULL;
#endif
}
/*  call the iterative solver: */

void Finley_Solver(Finley_SystemMatrix* A,double* x,double* b,Finley_SolverOptions* options) {
#if ITERATIVE_SOLVER == NO_LIB
    double norm2OfB,tol,tolerance,time0,time1,*r=NULL,norm_of_residual,last_norm_of_residual;
    int i,totIter,cntIter,finalizeIteration,errorCode;
    int n_col = A->num_cols * A-> col_block_size;
    int n_row = A->num_rows * A-> row_block_size;

    tolerance=MAX(options->tolerance,EPSILON);
    Finley_ErrorCode=NO_ERROR;
    /* check matrix type */
    if (A->type!=CSR) {
      Finley_ErrorCode = TYPE_ERROR;
      sprintf(Finley_ErrorMsg, "Iterative solver requires CSR format.");
      return;
    }
    if (A->symmetric) {
      Finley_ErrorCode = TYPE_ERROR;
      sprintf(Finley_ErrorMsg, "Iterative solver does not accept matrices stored in symmetric format.");
      return;
    }

    if (A->col_block_size != A->row_block_size) {
       Finley_ErrorCode = TYPE_ERROR;
       sprintf(Finley_ErrorMsg, "preconditioner requires row and column block sizes to be equal.");
       return;
     }
    /* check matrix is square */
    if (A->col_block_size != A->row_block_size) {
       Finley_ErrorCode = TYPE_ERROR;
       sprintf(Finley_ErrorMsg, "preconditioner requires row and column block sizes to be equal.");
       return;
     }
     if (A->num_cols != A->num_rows) {
       Finley_ErrorCode = TYPE_ERROR;
       sprintf(Finley_ErrorMsg, "This method requires a square matrix.");
       return;
     }
    /* get the norm of the right hand side */
    norm2OfB=0;
    #pragma omp parallel for private(i) reduction(+:norm2OfB) schedule(static)
    for (i = 0; i < n_row ; i++) norm2OfB += b[i] * b[i];
    norm2OfB=sqrt(norm2OfB);

    /* if norm2OfB==0 we are ready: x=0 */

    if (norm2OfB <=0.) {
        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n_col; i++) x[i]=0.;
        return;
    }
  

    /* construct the preconditioner */
    
    time1=Finley_timer();
    Finley_Solver_setPreconditioner(A,options);
    time1=Finley_timer()-time1;
    if (Finley_ErrorCode!=NO_ERROR) return;

    /* get an initial guess by evaluating the preconditioner */
    time0=Finley_timer();
    #pragma omp parallel
    Finley_Solver_solvePreconditioner(A,x,b);
    if (Finley_ErrorCode!=NO_ERROR) return;
    /* start the iteration process :*/
    r=(double*)TMPMEMALLOC(n_row*sizeof(double));
    if (Finley_checkPtr(r)) return;
    totIter = 0;
    finalizeIteration = FALSE;
    norm_of_residual=norm2OfB;
    while (! finalizeIteration) {
       cntIter = options->iter_max - totIter;
       finalizeIteration = TRUE;
       /*     Set initial residual. */
       norm_of_residual = 0;
       #pragma omp parallel
       {
          #pragma omp for private(i) schedule(static)
          for (i = 0; i < n_row; i++) r[i]=b[i];
          Finley_RawScaledSystemMatrixVector(DBLE(-1), A, x, DBLE(1), r);
          #pragma omp for private(i) reduction(+:norm_of_residual) schedule(static)
          for (i = 0; i < n_row; i++) norm_of_residual+= r[i] * r[i];
       }
       norm_of_residual =sqrt(norm_of_residual);
       if (totIter>0 && norm_of_residual>=last_norm_of_residual) {
            Finley_ErrorCode=WARNING;
	    sprintf(Finley_ErrorMsg, "No improvement during iteration.\nIterative solver gives up.");
       } else {
          /* */
          if (norm_of_residual>tolerance*norm2OfB) {
      	     if (totIter>0) {
                  tol=MAX(tolerance*tol/norm_of_residual,EPSILON*10.)*norm2OfB;
                  printf("new stopping criterion is %e.\n",tol);
      	     } else {
                  tol=tolerance*norm2OfB;
                  printf("Stopping criterion is %e.\n",tol);
      	     }
             last_norm_of_residual=norm_of_residual;
             /* call the solver */
             switch (options->iterative_method) {
                    default:
                        printf("Information: selceted iterative method not yet implemented. BiCGStab used instead.\n");
                    case BICGSTAB:
                        errorCode = Finley_Solver_BiCGStab(A, r, x, &cntIter, &tol); 
                        break;
                    case PCG:
                        errorCode = Finley_Solver_PCG(A, r, x, &cntIter, &tol); 
                        break;
              }
              totIter += cntIter;
         
              /* error handling  */
              if (errorCode==NO_ERROR) {
                    finalizeIteration = FALSE;
      	      } else if (errorCode==SOLVER_MAXITER_REACHED) {
                    Finley_ErrorCode=WARNING;
                    sprintf(Finley_ErrorMsg,"maximum number of iteration step reached.\nReturned solution does not fulfil stopping criterion.");
              } else if (errorCode == SOLVER_INPUT_ERROR ) {
                    Finley_ErrorCode=SYSTEM_ERROR;
         	   sprintf(Finley_ErrorMsg,"illegal dimension in iterative solver.");
              } else if ( errorCode == SOLVER_BREAKDOWN ) {
         	  if (cntIter <= 1) {
                       Finley_ErrorCode =ZERO_DIVISION_ERROR;
         	      sprintf(Finley_ErrorMsg, "fatal break down in iterative solver.");
                  } else {
         	      printf("Breakdown at iter %d (residual = %e). Restarting ...\n", cntIter+totIter, tol);
                      finalizeIteration = FALSE;
                  }
              } else if (errorCode == SOLVER_MEMORY_ERROR) {
         	      Finley_ErrorCode=MEMORY_ERROR;
         	      sprintf(Finley_ErrorMsg,"memory allocation failed.");
              } else if (errorCode !=SOLVER_NO_ERROR ) {
         	      Finley_ErrorCode=SYSTEM_ERROR;
         	      sprintf(Finley_ErrorMsg,"unidentified error in iterative solver.");
              }
         }
      }
    } /* while */
    MEMFREE(r);
    time0=Finley_timer()-time0;
    printf("Iterative solver reached tolerance %.5e after %d steps.\n",norm_of_residual,totIter);
    printf("timing: solver: %.4e sec\n",time1+time0);
    printf("timing: preconditioner: %.4e sec\n",time1);
    if (totIter>0) printf("timing: per iteration step: %.4e sec\n",time0/totIter);

#endif
}

/*
* $Log$
* Revision 1.3  2004/12/15 03:48:47  jgs
* *** empty log message ***
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1  2004/07/02 04:21:14  gross
* Finley C code has been included
*
*
*/
