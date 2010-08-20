
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

/* Paso: SystemMatrix: sets-up the preconditioner           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"
#include "PasoUtil.h"

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
      Paso_Solver_AMLI_free(in->amli);
      Paso_Solver_AMLI_System_free(in->amliSystem);
      MEMFREE(in);
    }
}


/*  call the iterative solver: */

void Paso_Solver_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options) {

    Paso_Solver_Preconditioner* prec=NULL;
    dim_t i;
    /*char filename[7];*/
    if (A->solver==NULL) {
        /* allocate structure to hold preconditioner */
        prec=MEMALLOC(1,Paso_Solver_Preconditioner);

        if (Paso_checkPtr(prec)) return;
        
        prec->type=UNKNOWN;
        prec->rilu=NULL;
        prec->ilu=NULL;
        prec->jacobi=NULL;
        prec->gs=NULL;
        prec->amg=NULL;
        prec->amli=NULL;
        prec->amgSystem=NULL;
        prec->amliSystem=NULL;
	if (options->verbose && options->use_local_preconditioner) printf("Apply preconditioner locally only.\n");

        A->solver=prec;
        switch (options->preconditioner) {
           default:
           case PASO_JACOBI:
              if (options->verbose) printf("Jacobi preconditioner is used.\n");
              prec->jacobi=Paso_Solver_getJacobi(A);
              prec->type=PASO_JACOBI;
              break;
	   case PASO_GS:
	      if (options->verbose) printf("Gauss-Seidel(%d) preconditioner is used.\n",options->sweeps);
	      prec->gs=Paso_Solver_getGS(A, options->sweeps, options->use_local_preconditioner, options->verbose);
	      prec->type=PASO_GS;
	      break;
	      
	   /***************************************************************************************/   
           case PASO_ILU0:
              if (options->verbose) printf("ILU preconditioner is used.\n");
              prec->ilu=Paso_Solver_getILU(A->mainBlock,options->verbose);
              prec->type=PASO_ILU0;
	      Paso_MPIInfo_noError(A->mpi_info);
              break;
           case PASO_RILU:
              if (options->verbose) printf("RILU preconditioner is used.\n");
              prec->rilu=Paso_Solver_getRILU(A->mainBlock,options->verbose);
	      Paso_MPIInfo_noError(A->mpi_info);
              prec->type=PASO_RILU;
              break;

            case PASO_AMG:
              if (options->verbose) printf("AMG preconditioner is used.\n");
              prec->amg=Paso_Solver_getAMG(A->mainBlock,options->level_max,options); 
	      Paso_MPIInfo_noError(A->mpi_info);
	      /*
              prec->amgSystem=MEMALLOC(1,Paso_Solver_AMG_System);
	      if (Paso_checkPtr(prec->amgSystem)) return;
	      prec->amgSystem->block_size=A->row_block_size;
	      for (i=0;i<A->row_block_size;++i) {
		prec->amgSystem->amgblock[i]=NULL;
		prec->amgSystem->block[i]=NULL;
	      }
              */
                       
              /*For performace reasons we check if block_size is one. If yes, then we do not need to separate blocks.*/
              /*if (A->row_block_size==1) {
                prec->amg=Paso_Solver_getAMG(A->mainBlock,options->level_max,options);  
              }
              else {
                for (i=0;i<A->row_block_size;++i) {
                prec->amgSystem->block[i]=Paso_SparseMatrix_getBlock(A->mainBlock,i+1);
                prec->amgSystem->amgblock[i]=Paso_Solver_getAMG(prec->amgSystem->block[i],options->level_max,options);
                }
              }
              */
                          
              prec->type=PASO_AMG;
              break;

            case PASO_AMLI:
	      
	      prec->amliSystem=MEMALLOC(1,Paso_Solver_AMLI_System);
              if (Paso_checkPtr(prec->amliSystem)) return;
                
              prec->amliSystem->block_size=A->row_block_size;
              
	      for (i=0;i<A->row_block_size;++i) {
		  prec->amliSystem->amliblock[i]=NULL;
		  prec->amliSystem->block[i]=NULL;
	      }

	      if (options->verbose) printf("AMLI preconditioner is used.\n");
              
              /*For performace reasons we check if block_size is one. If yes, then we do not need to separate blocks.*/
              if (A->row_block_size==1) {
                prec->amli=Paso_Solver_getAMLI(A->mainBlock,options->level_max,options);  
              }
              else {
                for (i=0;i<A->row_block_size;++i) {
                prec->amliSystem->block[i]=Paso_SparseMatrix_getBlock(A->mainBlock,i+1);
                prec->amliSystem->amliblock[i]=Paso_Solver_getAMLI(prec->amliSystem->block[i],options->level_max,options);
                }
              }
              prec->type=PASO_AMLI;
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
    dim_t n=A->mainBlock->numRows;

    


    switch (prec->type) {
        default:
        case PASO_JACOBI:
           Paso_Solver_solveJacobi(A, prec->jacobi,x,b);
           break;
	   
	case PASO_GS:
	   Paso_Solver_solveGS(A, prec->gs,x,b);
	   break;	   
        case PASO_ILU0:
           Paso_Solver_solveILU(A->mainBlock, prec->ilu,x,b);
           break;
        case PASO_RILU:
	   Paso_Solver_solveRILU(prec->rilu,x,b);
           break;

        case PASO_AMG:
            Paso_Solver_solveAMG(prec->amg,x,b);

            /*For performace reasons we check if block_size is one. If yes, then we do not need to do unnecessary copying.*/
            /*
            if (A->row_block_size<4) {
                Paso_Solver_solveAMG(prec->amg,x,b);
            }
            else {
           
                
                 for (i=0;i<A->row_block_size;i++) {
                    xx[i]=MEMALLOC(n,double);
                    bb[i]=MEMALLOC(n,double);
                    if (Paso_checkPtr(xx[i]) && Paso_checkPtr(bb[i])) return;
                }
                
                for (i=0;i<n;i++) {
                    for (j=0;j<A->row_block_size;j++) {
                     bb[j][i]=b[A->row_block_size*i+j];
                    }
                 }
                
                
                for (i=0;i<A->row_block_size;i++) {
                Paso_Solver_solveAMG(prec->amgSystem->amgblock[i],xx[i],bb[i]);
                }
                               
                for (i=0;i<n;i++) {
                    for (j=0;j<A->row_block_size;j++) {
                    x[A->row_block_size*i+j]=xx[j][i];
                    }
                 }
                
                for (i=0;i<A->row_block_size;i++) {
                MEMFREE(xx[i]);
                MEMFREE(bb[i]);
                }
            }
            */
        break;
        case PASO_AMLI:
            
            /*For performace reasons we check if block_size is one. If yes, then we do not need to do unnecessary copying.*/
            if (A->row_block_size==1) {
                Paso_Solver_solveAMLI(prec->amli,x,b);
            }
            else {
	       double **xx;
	       double **bb;
	       xx=MEMALLOC(A->row_block_size,double*);
	       bb=MEMALLOC(A->row_block_size,double*);
	       if (Paso_checkPtr(xx) || Paso_checkPtr(bb)) return;
                 for (i=0;i<A->row_block_size;i++) {
                    xx[i]=MEMALLOC(n,double);
                    bb[i]=MEMALLOC(n,double);
                    if (Paso_checkPtr(xx[i]) && Paso_checkPtr(bb[i])) return;
                }
                
                /*#pragma omp parallel for private(i,j) schedule(static)*/
                for (i=0;i<n;i++) {
                    for (j=0;j<A->row_block_size;j++) {
                     bb[j][i]=b[A->row_block_size*i+j];
                    }
                 }
                
                
                for (i=0;i<A->row_block_size;i++) {
                Paso_Solver_solveAMLI(prec->amliSystem->amliblock[i],xx[i],bb[i]);
                }
                               
                /*#pragma omp parallel for private(i,j) schedule(static)*/
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
