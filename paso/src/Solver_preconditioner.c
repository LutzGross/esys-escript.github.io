
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix: sets-up the preconditioner           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"

/***********************************************************************************/

/*  free space */

void Paso_Preconditioner_free(Paso_Solver_Preconditioner* in) {
    if (in!=NULL) {
      Paso_Solver_ILU_free(in->ilu);
      Paso_Solver_RILU_free(in->rilu);
      Paso_Solver_Jacobi_free(in->jacobi);
      Paso_Solver_GS_free(in->gs);
      Paso_Solver_AMG_free(in->amg);
      Paso_Solver_AMG_System_free(in->amgSystem);
      MEMFREE(in);
    }
}
/*  call the iterative solver: */

void Paso_Solver_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options) {

    Paso_Solver_Preconditioner* prec=NULL;
    dim_t i;
    if (A->solver==NULL) {
        /* allocate structure to hold preconditioner */
        prec=MEMALLOC(1,Paso_Solver_Preconditioner);
        prec->amgSystem=MEMALLOC(1,Paso_Solver_AMG_System);
        if (Paso_checkPtr(prec)) return;
        prec->type=UNKNOWN;
        prec->rilu=NULL;
        prec->ilu=NULL;
        prec->jacobi=NULL;
        prec->gs=NULL;
        prec->amg=NULL;
        
        prec->amgSystem->block_size=A->row_block_size;
        for (i=0;i<A->row_block_size;++i) {
          prec->amgSystem->amgblock[i]=NULL;
          prec->amgSystem->block[i]=NULL;
        }
        

        A->solver=prec;
        switch (options->preconditioner) {
           default:
           case PASO_JACOBI:
              if (options->verbose) printf("Jacobi preconditioner is used.\n");
              prec->jacobi=Paso_Solver_getJacobi(A->mainBlock);
              prec->type=PASO_JACOBI;
              break;
           case PASO_ILU0:
              if (options->verbose) printf("ILU preconditioner is used.\n");
              prec->ilu=Paso_Solver_getILU(A->mainBlock,options->verbose);
              prec->type=PASO_ILU0;
              break;
           case PASO_RILU:
              if (options->verbose) printf("RILU preconditioner is used.\n");
              prec->rilu=Paso_Solver_getRILU(A->mainBlock,options->verbose);
              prec->type=PASO_RILU;
              break;
            case PASO_GS:
              if (options->verbose) printf("Gauss-Seidel preconditioner is used.\n");
              prec->gs=Paso_Solver_getGS(A->mainBlock,options->verbose);
              prec->gs->sweeps=options->sweeps;
              prec->type=PASO_GS;
              break;
            case PASO_AMG:
              if (options->verbose) printf("AMG preconditioner is used.\n");
              
              /*For performace reasons we check if block_size is one. If yes, then we do not need to separate blocks.*/
              if (A->row_block_size==1) {
                prec->amg=Paso_Solver_getAMG(A->mainBlock,options->level_max,options);  
              }
              else {
                #pragma omp parallel for private(i) schedule(static)
                for (i=0;i<A->row_block_size;++i) {
                prec->amgSystem->block[i]=Paso_SparseMatrix_getBlock(A->mainBlock,i+1);
                prec->amgSystem->amgblock[i]=Paso_Solver_getAMG(prec->amgSystem->block[i],options->level_max,options);
                }
              }

              prec->type=PASO_AMG;
              break;
 
        }
        if (! Paso_MPIInfo_noError(A->mpi_info ) ){
            Paso_Preconditioner_free(prec);
            A->solver=NULL;
        }
    }
}

/* applies the preconditioner */
/* has to be called within a parallel reqion */
/* barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Paso_Solver_solvePreconditioner(Paso_SystemMatrix* A,double* x,double* b){
    Paso_Solver_Preconditioner* prec=(Paso_Solver_Preconditioner*) A->solver;
    dim_t i,j;
    dim_t n=Paso_SystemMatrix_getGlobalNumRows(A);
    double **xx;
    double **bb;

    switch (prec->type) {
        default:
        case PASO_JACOBI:
           Paso_Solver_solveJacobi(prec->jacobi,x,b);
           break;
        case PASO_ILU0:
           Paso_Solver_solveILU(prec->ilu,x,b);
           break;
        case PASO_RILU:
           Paso_Solver_solveRILU(prec->rilu,x,b);
           break;
        case PASO_GS:
            /* Gauss-Seidel preconditioner P=U^{-1}DL^{-1} is used here with sweeps paramenter.
              We want to solve x_new=x_old+P^{-1}(b-Ax_old). So for fisrt 3 we will have the following:
              x_0=0
              x_1=P^{-1}(b)
              
              b_old=b
              
              b_new=b_old+b
              x_2=x_1+P^{-1}(b-Ax_1)=P^{-1}(b)+P^{-1}(b)-P^{-1}AP^{-1}(b))
                 =P^{-1}(2b-AP^{-1}b)=P^{-1}(b_new-AP^{-1}b_old)
              b_old=b_new
              
              b_new=b_old+b
              x_3=....=P^{-1}(b_new-AP^{-1}b_old)
              b_old=b_new
              
              So for n- steps we will use loop, but we always make sure that every value calculated only once!
            */

           Paso_Solver_solveGS(prec->gs,x,b);
           if (prec->gs->sweeps>1) {
                double *bold=MEMALLOC(prec->gs->n*prec->gs->n_block,double);
                double *bnew=MEMALLOC(prec->gs->n*prec->gs->n_block,double);
           
                #pragma omp parallel for private(i) schedule(static)
                for (i=0;i<prec->gs->n*prec->gs->n_block;++i) bold[i]=b[i];
           
                while(prec->gs->sweeps>1) {
                   #pragma omp parallel for private(i) schedule(static)
                   for (i=0;i<prec->gs->n*prec->gs->n_block;++i) bnew[i]=bold[i]+b[i];
                    /* Compute the residual b=b-Ax*/
                   Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1), A, x, DBLE(1), bnew);
                   /* Go round again*/
                   Paso_Solver_solveGS(prec->gs,x,bnew);
                   #pragma omp parallel for private(i) schedule(static)
                   for (i=0;i<prec->gs->n*prec->gs->n_block;++i) bold[i]=bnew[i];
                   prec->gs->sweeps=prec->gs->sweeps-1;
                }
                MEMFREE(bold);
                MEMFREE(bnew); 
           }
           break;
        case PASO_AMG:
            
            /*For performace reasons we check if block_size is one. If yes, then we do not need to do unnecessary copying.*/
            if (A->row_block_size==1) {
                Paso_Solver_solveAMG(prec->amg,x,b);
            }
            else {
           
                xx=MEMALLOC(A->row_block_size,double*);
                bb=MEMALLOC(A->row_block_size,double*);
                
                for (i=0;i<A->row_block_size;i++) {
                    xx[i]=MEMALLOC(n,double);
                    bb[i]=MEMALLOC(n,double);
                }

                #pragma omp parallel for private(i,j) schedule(static)
                for (i=0;i<n;i++) {
                    for (j=0;j<A->row_block_size;j++) {
                     bb[j][i]=b[A->row_block_size*i+j];
                    }
                 }
                                
                #pragma omp parallel for private(i) schedule(static)
                for (i=0;i<A->row_block_size;i++) {
                Paso_Solver_solveAMG(prec->amgSystem->amgblock[i],xx[i],bb[i]);
                }
                   
                #pragma omp parallel for private(i,j) schedule(static)
                for (i=0;i<n;i++) {
                    for (j=0;j<A->row_block_size;j++) {
                    x[A->row_block_size*i+j]=xx[j][i];
                    }
                 }
                               
                for (i=0;i<A->row_block_size;i++) {
                MEMFREE(xx[i]);
                MEMFREE(bb[i]);
                }
                
                MEMFREE(xx);
                MEMFREE(bb);
            }
            
        break;
    }

}
