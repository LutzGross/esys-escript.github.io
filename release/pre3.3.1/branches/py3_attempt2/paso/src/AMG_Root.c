
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

/* Paso: AMG set-ups                                          */

/**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Preconditioner.h"
#include "BOOMERAMG.h"

/***********************************************************************************/

/*  free space */

void Paso_Preconditioner_AMG_Root_free(Paso_Preconditioner_AMG_Root* in) {
    if (in!=NULL) {
      Paso_Preconditioner_AMG_free(in->amg);
      Paso_Preconditioner_LocalAMG_free(in->localamg);
      Paso_Preconditioner_BoomerAMG_free(in->boomeramg);
      Paso_Preconditioner_Smoother_free(in->amgsubstitute);
      MEMFREE(in);
    }
}

Paso_Preconditioner_AMG_Root* Paso_Preconditioner_AMG_Root_alloc(Paso_SystemMatrix* A,Paso_Options* options) {
    Paso_Preconditioner_AMG_Root* prec=MEMALLOC(1,Paso_Preconditioner_AMG_Root);

    if (! Esys_checkPtr(prec)) {
        
	prec->amg=NULL;
	prec->localamg=NULL;
	prec->amgsubstitute=NULL;
        prec->boomeramg=NULL;
	if (options->preconditioner == PASO_BOOMERAMG){
	  prec->boomeramg = Paso_Preconditioner_BoomerAMG_alloc(A,options);
	} else {
          prec->is_local=( A->mpi_info->size == 1 ) | options->use_local_preconditioner;
        
          if (prec->is_local) {
	      prec->localamg=Paso_Preconditioner_LocalAMG_alloc(A->mainBlock,1,options);
	      Esys_MPIInfo_noError(A->mpi_info);
          } else {
              prec->amg=Paso_Preconditioner_AMG_alloc(A,1,options);
          }
	}
	if ( Esys_noError() ) {
	      if (options->verbose) {
                if ( (prec->localamg != NULL) || ( prec->amg != NULL) ||
		     (prec->boomeramg != NULL) )  {
	           printf("Paso_Preconditioner_AMG_Root:  Smoother is ");
	           if (options->smoother == PASO_JACOBI) {
		           printf("Jacobi");
		   } else {
		           printf("Gauss-Seidel");
		   }
		   printf(" with %d/%d pre/post sweeps",options->pre_sweeps, options->post_sweeps);
                   if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION) {
                               printf( " and classical interpolation.\n");
                   } else if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) {
                               printf( " and classical interpolation with enforced FF coupling.\n");
                   } else {
                               printf( " and direct interpolation.\n");
                   }
                } else {
	           printf("Paso_Preconditioner_AMG_Root:  no coarsening constructed.\n");
                }
             }

             if (prec->localamg != NULL) {
		     options->num_level=Paso_Preconditioner_LocalAMG_getMaxLevel(prec->localamg);
		     options->coarse_level_sparsity=Paso_Preconditioner_LocalAMG_getCoarseLevelSparsity(prec->localamg);
		     options->num_coarse_unknowns=Paso_Preconditioner_LocalAMG_getNumCoarseUnknwons(prec->localamg);
             } else if ( prec->amg != NULL) {
		     options->num_level=Paso_Preconditioner_AMG_getMaxLevel(prec->amg);
		     options->coarse_level_sparsity=Paso_Preconditioner_AMG_getCoarseLevelSparsity(prec->amg);
		     options->num_coarse_unknowns=Paso_Preconditioner_AMG_getNumCoarseUnknwons(prec->amg);
             } else if (prec->boomeramg == NULL) {
                     prec->sweeps=options->sweeps;
		     prec->amgsubstitute=Paso_Preconditioner_Smoother_alloc(A, (options->smoother == PASO_JACOBI), prec->is_local, options->verbose);
		     options->num_level=0;
		     if (options->verbose) { 
			if (options->smoother == PASO_JACOBI) {
			   printf("Paso_Preconditioner: Jacobi(%d) preconditioner is used.\n",prec->sweeps);
			} else {
			   printf("Paso_Preconditioner: Gauss-Seidel(%d) preconditioner is used.\n",prec->sweeps);
			}
		     }
	      }
         }
      }
      if (! Esys_noError() ){
             Paso_Preconditioner_AMG_Root_free(prec);
	     return NULL;
      } else {
	    return prec;
      }
}

/* Applies the preconditioner. */
/* Has to be called within a parallel region. */
/* Barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Paso_Preconditioner_AMG_Root_solve(Paso_SystemMatrix* A, Paso_Preconditioner_AMG_Root* prec, double* x,double* b)
{
             if (prec->localamg != NULL) {
	             Paso_Preconditioner_LocalAMG_solve(A->mainBlock, prec->localamg,x,b);
             } else if ( prec->amg != NULL) {
	             Paso_Preconditioner_AMG_solve(A, prec->amg,x,b);
	     } else if ( prec->boomeramg != NULL) {
		     Paso_Preconditioner_BoomerAMG_solve(A, prec->boomeramg,x,b);
             } else {
	             Paso_Preconditioner_Smoother_solve(A, prec->amgsubstitute,x,b,prec->sweeps,FALSE);
	    }
}

