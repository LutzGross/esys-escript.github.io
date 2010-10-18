
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
      Paso_Preconditioner_LocalSmoother_free(in->localamgsubstitute);
      /*********************************/
      Paso_Solver_ILU_free(in->ilu);
      Paso_Solver_RILU_free(in->rilu);
      /*********************************/
      
      MEMFREE(in);
    }
}

Paso_Preconditioner* Paso_Preconditioner_alloc(Paso_SystemMatrix* A,Paso_Options* options) {

    Paso_Preconditioner* prec=NULL;

    prec=MEMALLOC(1,Paso_Preconditioner);

    if (! Esys_checkPtr(prec)) {
        
        prec->type=UNKNOWN;
	
	
	prec->jacobi=NULL;
	prec->gs=NULL;
	prec->localamg=NULL;
	prec->localamgsubstitute=NULL;
	
	/*********************************/	
        prec->rilu=NULL;
        prec->ilu=NULL;
	/*********************************/
	
	if (options->verbose && options->use_local_preconditioner) printf("Paso: Apply preconditioner locally only.\n");

        switch (options->preconditioner) {
           default:
           case PASO_JACOBI:
	      if (options->verbose) printf("Paso: Jacobi(%d) preconditioner is used.\n",options->sweeps);
	      prec->jacobi=Paso_Preconditioner_Smoother_alloc(A, TRUE, options->use_local_preconditioner, options->verbose);
              prec->type=PASO_JACOBI;
	      prec->sweeps=options->sweeps;
              break;
	   case PASO_GS:
	      if (options->verbose) printf("Paso: Gauss-Seidel(%d) preconditioner is used.\n",options->sweeps);
	      prec->gs=Paso_Preconditioner_Smoother_alloc(A, FALSE, options->use_local_preconditioner, options->verbose);
	      prec->type=PASO_GS;
	      prec->sweeps=options->sweeps;
	      break;
	   case PASO_AMLI:
	   case PASO_AMG:
	      if (options->verbose) printf("Paso: AMG preconditioner is used.\n");
	      prec->localamg=Paso_Preconditioner_LocalAMG_alloc(A->mainBlock,1,options);
	      prec->sweeps=options->sweeps;
	      /* if NULL is returned (and no error) no AMG has been constructed because the system is too small or not big enough 
	         we now use the Smoother as a preconditioner                                                                      */
	      if ( (Esys_noError()) && (prec->localamg == NULL) ) {
		 if (options->verbose) { 
		    if (options->smoother == PASO_JACOBI) {
			printf("Paso: Jacobi(%d) preconditioner is used.\n",prec->sweeps);
		    } else {
		       printf("Paso: Gauss-Seidel(%d) preconditioner is used.\n",prec->sweeps);
		    }
		 }
		 prec->localamgsubstitute=Paso_Preconditioner_LocalSmoother_alloc(A->mainBlock, (options->smoother == PASO_JACOBI), options->verbose);
	      } else {
		 if (options->verbose) printf("Paso: AMG preconditioner is used.\n");
	      }
	      prec->type=PASO_AMG;
	      Esys_MPIInfo_noError(A->mpi_info);
	      break;
	      
	   /***************************************************************************************/   
           case PASO_ILU0:
	      if (options->verbose) printf("Paso: ILU preconditioner is used.\n");
              prec->ilu=Paso_Solver_getILU(A->mainBlock,options->verbose);
              prec->type=PASO_ILU0;
	      Esys_MPIInfo_noError(A->mpi_info);
              break;
           case PASO_RILU:
	      if (options->verbose) printf("Paso: RILU preconditioner is used.\n");
              prec->rilu=Paso_Solver_getRILU(A->mainBlock,options->verbose);
	      Esys_MPIInfo_noError(A->mpi_info);
              prec->type=PASO_RILU;
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

    switch (prec->type) {
        default:
        case PASO_JACOBI:
	   Paso_Preconditioner_Smoother_solve(A, prec->jacobi,x,b,prec->sweeps, FALSE);
           break;   
	case PASO_GS:
	   Paso_Preconditioner_Smoother_solve(A, prec->gs,x,b,prec->sweeps, FALSE);
	   break;
	case PASO_AMG:
	   if (prec->localamg == NULL) {
	      Paso_Preconditioner_LocalSmoother_solve(A->mainBlock, prec->localamgsubstitute,x,b,prec->sweeps, FALSE);
	   } else {
	      Paso_Preconditioner_LocalAMG_solve(A->mainBlock, prec->localamg,x,b);
	   }
	   break;
	   
	/*=========================================================*/   
        case PASO_ILU0:
           Paso_Solver_solveILU(A->mainBlock, prec->ilu,x,b);
           break;
        case PASO_RILU:
	   Paso_Solver_solveRILU(prec->rilu,x,b);
           break;

    }

}
