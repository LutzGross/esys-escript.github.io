
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: interface to the direct solvers                      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003                        */
/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/**************************************************************/

#include "Paso.h"
#include "performance.h"
#include "Preconditioner.h"
#include "Solver.h"
#include "MKL.h"
#include "UMFPACK.h"

/**************************************************************/

void Paso_solve(Paso_SystemMatrix* A,
                double* out,
                double* in,
                Paso_Options* options) {

  Paso_Performance pp;
  index_t package;
  Esys_resetError();
  if (Paso_SystemMatrix_getGlobalNumCols(A) != Paso_SystemMatrix_getGlobalNumRows(A)
                || A->col_block_size!=A->row_block_size) {
       Esys_setError(VALUE_ERROR,"Paso_solve: matrix has to be a square matrix.");
       return;
  }
  /* Paso_Options_show(options); */
  Performance_open(&pp,options->verbose);
  package=Paso_Options_getPackage(options->method,options->package,options->symmetric, A->mpi_info);
  if (Esys_noError()) {
     switch(package) {

        case PASO_PASO:
          Paso_Solver(A,out,in,options,&pp);
          A->solver_package=PASO_PASO;
          break;

        case PASO_MKL:
          if (A->mpi_info->size>1) {
              Esys_setError(VALUE_ERROR,"Paso_solve: MKL package does not support MPI.");
              return;
          }
          options->converged=FALSE;
	  options->time=Esys_timer();
	  Performance_startMonitor(&pp,PERFORMANCE_ALL);
	  Paso_MKL(A->mainBlock, out, in, options->reordering, options->refinements, options->verbose);
          A->solver_package=PASO_MKL;
	  Performance_stopMonitor(&pp,PERFORMANCE_ALL);
	  options->time=Esys_timer()-options->time;
	  options->set_up_time=0;
	  options->residual_norm=0.;
	  options->num_iter=0;
	  if (Esys_MPIInfo_noError(A->mpi_info)) options->converged=TRUE;
          break;

	 case PASO_UMFPACK:
	    if (A->mpi_info->size>1) {
	       Esys_setError(VALUE_ERROR,"Paso_solve: UMFPACK package does not support MPI.");
	       return;
	    }
	    options->converged=FALSE;
	    options->time=Esys_timer();
	    Performance_startMonitor(&pp,PERFORMANCE_ALL);
	    Paso_UMFPACK(A->mainBlock, out, in, options->refinements, options->verbose);
	    A->solver_package=PASO_UMFPACK;
	    Performance_stopMonitor(&pp,PERFORMANCE_ALL);
	    options->time=Esys_timer()-options->time;
	    options->set_up_time=0;
	    options->residual_norm=0.;
	    options->num_iter=0;
	    if (Esys_MPIInfo_noError(A->mpi_info)) options->converged=TRUE;
            break;

        default:
           Esys_setError(VALUE_ERROR,"Paso_solve: unknown package code");
           break;
     }
  }
  /*
       cancel divergence errors
  */
  if (options->accept_failed_convergence) {
         if (Esys_getErrorType() == DIVERGED) {
             Esys_resetError();
             if (options->verbose) printf("PASO: failed convergence error has been canceled as requested.\n");
         } 
  }
  Performance_close(&pp,options->verbose);
  /* Paso_Options_showDiagnostics(options); */
  return;
}

/*  free memory possibly reserved for a recall */

void Paso_solve_free(Paso_SystemMatrix* in) { 

     if (in==NULL) return;

     switch(in->solver_package) {

        case PASO_PASO:
          Paso_Solver_free(in);
          break;

	case PASO_SMOOTHER:
	  Paso_Preconditioner_Smoother_free((Paso_Preconditioner_Smoother*) in->solver_p);
	  break;
	  
        case PASO_MKL:
          Paso_MKL_free(in->mainBlock); 
          break;

        case PASO_UMFPACK:
          Paso_UMFPACK_free(in->mainBlock); 
          break;

   }
}
