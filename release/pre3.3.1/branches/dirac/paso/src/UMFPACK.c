
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

/* Paso: interface to the Intel UMFPACK library */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "UMFPACK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

/*  free any extra stuff possibly used by the UMFPACK library */

void Paso_UMFPACK_free(Paso_SparseMatrix* A) {
     Paso_UMFPACK_Handler* pt =NULL;
     if ( (A->solver_p!=NULL) && (A->solver_package == PASO_UMFPACK) ) {
           pt=(Paso_UMFPACK_Handler*)(A->solver_p);
	   #ifdef UMFPACK
	   umfpack_di_free_symbolic(&(pt->symbolic));
	   umfpack_di_free_numeric(&(pt->numeric));
	   #endif
	   MEMFREE(pt);
           A->solver_p=NULL;
     }
}



/*  call the solver: */

void Paso_UMFPACK(Paso_SparseMatrix* A,
                          double* out,
                          double* in,
		          dim_t numRefinements,
		          bool_t verbose) 
{
#ifdef UMFPACK
     double time0;
     Paso_UMFPACK_Handler* pt = NULL;
     double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
     int error = UMFPACK_OK;
     
     if (! ( (A->type & MATRIX_FORMAT_BLK1) && (A->type & MATRIX_FORMAT_CSC)) ) {
        Esys_setError(TYPE_ERROR,"Paso_UMFPACK: UMFPACK requires CSC format with index offset 1 and block size 1.");
        return;
     }
     
     pt = (Paso_UMFPACK_Handler *)(A->solver_p);
     umfpack_di_defaults(control);
     
     if (pt==NULL) {
	int n = A->numRows;
	pt=(MEMALLOC(1,Paso_UMFPACK_Handler));
	time0=Esys_timer();
	if (Esys_checkPtr(pt)) return;
	A->solver_p = (void*) pt;
        A->solver_package=PASO_UMFPACK;
	/* call LDU symbolic factorization: */
	error=umfpack_di_symbolic(n,n,A->pattern->ptr,A->pattern->index,A->val,&(pt->symbolic),control,info);
	if (error != UMFPACK_OK) {
	   if (verbose) printf("UMFPACK: symbolic factorization factorization failed.\n");
	   Esys_setError(VALUE_ERROR,"UMFPACK: symbolic factorization failed.");
	   return;
	} else {
	   /* call LDU factorization: */
	   error= umfpack_di_numeric(A->pattern->ptr,A->pattern->index,A->val,pt->symbolic,&(pt->numeric),control,info);
	   if (error != UMFPACK_OK) {
	      if (verbose) printf("UMFPACK: LDU factorization failed.\n");
	      Esys_setError(ZERO_DIVISION_ERROR,"factorization failed. Most likely the matrix is singular.");
	      return;
	   }
	   if (verbose) printf("UMFPACK: LDU factorization completed (time = %e).\n",Esys_timer()-time0);
	}
     }
     if (Esys_noError())  {
	/* call forward backward substitution: */
	control[UMFPACK_IRSTEP]=numRefinements; /* number of refinement steps */
	time0=Esys_timer();
	error=umfpack_di_solve(UMFPACK_A,A->pattern->ptr,A->pattern->index,A->val,out,in,pt->numeric,control,info);
	if (error != UMFPACK_OK) {
	   if (verbose) printf("UMFPACK: forward/backward substitution failed.\n");
	   Esys_setError(VALUE_ERROR,"forward/backward substition failed. Most likely the matrix is singular.");
	   return;
	}
	if (verbose) printf("UMFPACK: forward/backward substitution completed (time = %e).\n",Esys_timer()-time0);
     }
     #else
	 Esys_setError(SYSTEM_ERROR,"Paso_UMFPACK:UMFPACK is not avialble.");
     #endif
}
