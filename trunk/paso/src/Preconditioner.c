
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
#include "PasoUtil.h"
#include "Preconditioner.h"

/***********************************************************************************/

/*  free space */

void Paso_Preconditioner_free(Paso_Preconditioner* in) {
    if (in!=NULL) {
      Paso_Preconditioner_Smoother_free(in->jacobi);
      Paso_Preconditioner_Smoother_free(in->gs);
      Paso_Preconditioner_LocalAMG_free(in->localamg);
      
      /*********************************/
      Paso_Solver_ILU_free(in->ilu);
      Paso_Solver_RILU_free(in->rilu);
      Paso_Solver_AMLI_free(in->amli);
      Paso_Solver_AMLI_System_free(in->amliSystem);
      /*********************************/
      
      MEMFREE(in);
    }
}

Paso_Preconditioner* Paso_Preconditioner_alloc(Paso_SystemMatrix* A,Paso_Options* options) {

    Paso_Preconditioner* prec=NULL;
    dim_t i;

    prec=MEMALLOC(1,Paso_Preconditioner);

    if (! Esys_checkPtr(prec)) {
        
        prec->type=UNKNOWN;
	
	
	prec->jacobi=NULL;
	prec->gs=NULL;
	prec->localamg=NULL;
	/*********************************/	
        prec->rilu=NULL;
        prec->ilu=NULL;

        prec->amli=NULL;
        prec->amliSystem=NULL;
	/*********************************/
	
	if (options->verbose && options->use_local_preconditioner) printf("Apply preconditioner locally only.\n");

        switch (options->preconditioner) {
           default:
           case PASO_JACOBI:
	      if (options->verbose) printf("Jacobi(%d) preconditioner is used.\n",options->sweeps);
	      prec->jacobi=Paso_Preconditioner_Smoother_alloc(A, TRUE, options->use_local_preconditioner, options->verbose);
              prec->type=PASO_JACOBI;
	      prec->sweeps=options->sweeps;
              break;
	   case PASO_GS:
	      if (options->verbose) printf("Gauss-Seidel(%d) preconditioner is used.\n",options->sweeps);
	      prec->gs=Paso_Preconditioner_Smoother_alloc(A, FALSE, options->use_local_preconditioner, options->verbose);
	      prec->type=PASO_GS;
	      prec->sweeps=options->sweeps;
	      break;
	   case PASO_AMG:
	      if (options->verbose) printf("AMG preconditioner is used.\n");
	      prec->localamg=Paso_Preconditioner_LocalAMG_alloc(A->mainBlock,options->level_max,options); 
	      Esys_MPIInfo_noError(A->mpi_info);
	      prec->type=PASO_AMG;
	      break;
	      
	   /***************************************************************************************/   
           case PASO_ILU0:
              if (options->verbose) printf("ILU preconditioner is used.\n");
              prec->ilu=Paso_Solver_getILU(A->mainBlock,options->verbose);
              prec->type=PASO_ILU0;
	      Esys_MPIInfo_noError(A->mpi_info);
              break;
           case PASO_RILU:
              if (options->verbose) printf("RILU preconditioner is used.\n");
              prec->rilu=Paso_Solver_getRILU(A->mainBlock,options->verbose);
	      Esys_MPIInfo_noError(A->mpi_info);
              prec->type=PASO_RILU;
              break;



            case PASO_AMLI:
	      
	      prec->amliSystem=MEMALLOC(1,Paso_Solver_AMLI_System);
              if (! Esys_checkPtr(prec->amliSystem)) {
                
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
	      }
              break;
 
        }
    }
    if (! Esys_MPIInfo_noError(A->mpi_info ) ){
         Paso_Preconditioner_free(prec);
	return NULL;
    } else {
	return prec;
    }
}

/* applies the preconditioner */
/* has to be called within a parallel reqion */
/* barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Paso_Preconditioner_solve(Paso_Preconditioner* prec, Paso_SystemMatrix* A,double* x,double* b){
    dim_t i,j;
    dim_t n=A->mainBlock->numRows;
    switch (prec->type) {
        default:
        case PASO_JACOBI:
	   Paso_Preconditioner_Smoother_solve(A, prec->jacobi,x,b,prec->sweeps, FALSE);
           break;   
	case PASO_GS:
	   Paso_Preconditioner_Smoother_solve(A, prec->gs,x,b,prec->sweeps, FALSE);
	   break;
	case PASO_AMG:
	   Paso_Preconditioner_LocalAMG_solve(prec->localamg,x,b);
	   break;
	   
	/*=========================================================*/   
        case PASO_ILU0:
           Paso_Solver_solveILU(A->mainBlock, prec->ilu,x,b);
           break;
        case PASO_RILU:
	   Paso_Solver_solveRILU(prec->rilu,x,b);
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
	       if (Esys_checkPtr(xx) || Esys_checkPtr(bb)) return;
                 for (i=0;i<A->row_block_size;i++) {
                    xx[i]=MEMALLOC(n,double);
                    bb[i]=MEMALLOC(n,double);
                    if (Esys_checkPtr(xx[i]) && Esys_checkPtr(bb[i])) return;
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
	/*=========================================================*/
    }

}
