
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
/* Author: artak@access.edu.au */

/**************************************************************/

#include "paso/Paso.h"
#include "paso/Solver.h"
#include "paso/SystemMatrix.h"
#include "paso/Options.h"
#include "Paso_tests.h"


/**************************************************************/

void Paso_test_run(Paso_SystemMatrix* A,double* b,dim_t level) {

  Paso_Options options;
  Paso_Options_setDefaults(&options);
   
 if(level==1) /* Solvers only*/
    {
      options.method=PASO_PCG;
      options.verbose=TRUE;
      options.preconditioner=PASO_JACOBI;
      fprintf(stdout,"Test solver: PCG with JACOBI\n");
      Paso_test_matrix(A,b,&options);
      
      fprintf(stdout,"Test solver: GMRES with JACOBI\n");
      A->solver=NULL;
      options.method=PASO_GMRES;
      Paso_test_matrix(A,b,&options);
      
      fprintf(stdout,"Test solver: PRES20 with JACOBI\n");
      A->solver=NULL;
      options.method=PASO_PRES20;
      Paso_test_matrix(A,b,&options); 
      
      fprintf(stdout,"Test solver: MINRES with JACOBI\n");
      A->solver=NULL;
      options.method=PASO_MINRES;
      Paso_test_matrix(A,b,&options); 
      
      fprintf(stdout,"Test solver: TFQMR with JACOBI\n");
      A->solver=NULL;
      options.method=PASO_TFQMR;
      Paso_test_matrix(A,b,&options); 
    }
 else if (level==2) /* Preconditiones only with default solver*/
    {
      Paso_Options_setDefaults(&options);
      options.method=PASO_DEFAULT;
      options.verbose=TRUE;
      options.preconditioner=PASO_JACOBI;
      fprintf(stdout,"Test preconditioner: PASO_DEFAULT with JACOBI\n");
      Paso_test_matrix(A,b,&options);
      
      Paso_Options_setDefaults(&options);
      A->solver=NULL;
      options.verbose=TRUE;
      fprintf(stdout,"Test preconditioner: PASO_DEFAULT with ILU\n");
      options.method=PASO_DEFAULT;
      options.preconditioner=PASO_ILU0;
      Paso_test_matrix(A,b,&options);
      
      Paso_Options_setDefaults(&options);
      A->solver=NULL;
      options.verbose=TRUE;
      fprintf(stdout,"Test preconditioner: PASO_DEFAULT with RILU\n");
      options.method=PASO_DEFAULT;
      options.preconditioner=PASO_RILU;
      Paso_test_matrix(A,b,&options); 

      Paso_Options_setDefaults(&options);
      A->solver=NULL;
      options.verbose=TRUE;
      fprintf(stdout,"Test preconditioner: PASO_DEFAULT with GS\n");
      options.method=PASO_DEFAULT;
      options.preconditioner=PASO_GS;
      Paso_test_matrix(A,b,&options); 

      Paso_Options_setDefaults(&options);
      A->solver=NULL;
      options.verbose=TRUE;
      fprintf(stdout,"Test preconditioner: PASO_DEFAULT with AMG\n");
      options.method=PASO_DEFAULT;
      options.preconditioner=PASO_AMG;
      Paso_test_matrix(A,b,&options);  
      
    }
}

void Paso_test_matrix(Paso_SystemMatrix* A, double* b, Paso_Options* options ) {
   
   dim_t n=Paso_SystemMatrix_getTotalNumRows(A);
   double *out=NULL;
   out=MEMALLOC(n,double);
   
   if (Paso_checkPtr(out)) {
      fprintf(stderr,"Cannot allocate memory\n");
      return;
   }
   Paso_solve(A,out,b,options);
   
   MEMFREE(out);
   
}

void Paso_test_data(char *fileName_p, double* b, Paso_Options* options ) {
   
   Paso_SystemMatrix* A=NULL;
   dim_t n=Paso_SystemMatrix_getTotalNumRows(A);
   double *out=MEMALLOC(n,double);
   A=Paso_SystemMatrix_loadMM_toCSR(fileName_p);

   if (Paso_checkPtr(out)) {
      return;
   }
   
   Paso_solve(A,out,b,options);
   
   Paso_SystemMatrix_free(A);
   MEMFREE(out);
}
