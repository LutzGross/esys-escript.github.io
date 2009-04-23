
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: interface to the direct solvers                    */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "performance.h"
#include "Solver.h"

#ifdef SCSL
#include "SCSL.h"
#endif

#ifdef MKL
#include "MKL.h"
#endif

#ifdef UMFPACK
#include "UMFPACK.h"
#endif

/**************************************************************/

void Paso_solve(Paso_SystemMatrix* A,
                double* out,
                double* in,
                Paso_Options* options) {

  Paso_Performance pp;
  index_t package;
  Paso_resetError();
  if (Paso_SystemMatrix_getGlobalNumCols(A) != Paso_SystemMatrix_getGlobalNumRows(A)
                || A->col_block_size!=A->row_block_size) {
       Paso_setError(VALUE_ERROR,"Paso_solve: matrix has to be a square matrix.");
       return;
  }
  Performance_open(&pp,options->verbose);
  package=Paso_Options_getPackage(options->method,options->package,options->symmetric, A->mpi_info);
  if (Paso_noError()) {
     switch(package) {

        case PASO_PASO:
          Paso_Solver(A,out,in,options,&pp);
          A->solver_package=PASO_PASO;
          break;

        #ifdef SCSL
        case PASO_SCSL:
          if (A->mpi_info->size>1) {
              Paso_setError(VALUE_ERROR,"Paso_solve: SCSL package does not support MPI.");
              return;
          }
          Paso_SCSL(A,out,in,options,&pp);
          A->solver_package=PASO_SCSL;
          break;
        #endif

  
        #ifdef MKL
        case PASO_MKL:
          if (A->mpi_info->size>1) {
              Paso_setError(VALUE_ERROR,"Paso_solve: MKL package does not support MPI.");
              return;
          }
          Paso_MKL(A,out,in,options,&pp);
          A->solver_package=PASO_MKL;
          break;
        #endif

        #ifdef UMFPACK
        case PASO_UMFPACK:
          if (A->mpi_info->size>1) {
              Paso_setError(VALUE_ERROR,"Paso_solve: UMFPACK package does not support MPI.");
              return;
          }
          Paso_UMFPACK(A,out,in,options,&pp);
          A->solver_package=PASO_UMFPACK;
          break;
        #endif

        default:
           Paso_setError(VALUE_ERROR,"Paso_solve: unknown package code");
           break;
     }
  }
  Performance_close(&pp,options->verbose);
  return;
}

/*  free memory possibly resereved for a recall */

void Paso_solve_free(Paso_SystemMatrix* in) { 

     if (in==NULL) return;

     switch(in->solver_package) {

        case PASO_PASO:
          Paso_Solver_free(in);
          break;

        #ifdef SCSL
        case PASO_SCSL:
          Paso_SCSL_free(in);
          break;
        #endif

  
        #ifdef MKL
        case PASO_MKL:
          Paso_MKL_free(in); 
          break;
        #endif

        #ifdef UMFPACK
        case PASO_UMFPACK:
          Paso_UMFPACK_free(in); 
          break;
        #endif

   }
}
