
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/* Paso: interface to the Intel UMFPACK library */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "UMFPACK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

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
	if (error == UMFPACK_ERROR_out_of_memory ) {
	    if (verbose) printf("UMFPACK: symbolic factorization failed because of memory overlow.\n");
	   Esys_setError(MEMORY_ERROR,"UMFPACK: symbolic factorization failed because of memory overlow.");
	   return;
	} else if (error == UMFPACK_WARNING_singular_matrix ) {
	   if (verbose) printf("UMFPACK: symbolic factorization failed because of singular matrix.\n");
	   Esys_setError(ZERO_DIVISION_ERROR,"UMFPACK: symbolic factorization failed because of singular matrix.");
	   return;
	} else if ( (error == UMFPACK_WARNING_determinant_underflow ) || ( error == UMFPACK_WARNING_determinant_overflow ) ) {
	   if (verbose) printf("UMFPACK: symbolic factorization failed because of under/overflow.\n");
	   Esys_setError(FLOATING_POINT_ERROR,"UMFPACK: symbolic factorization failed because of under/overflow.");
	   return;	 
       }  else if (error != UMFPACK_OK) {
	   if (verbose) printf("UMFPACK: symbolic factorization failed. Error code = %d.\n",error);
	   Esys_setError(SYSTEM_ERROR,"UMFPACK: symbolic factorization failed.");
	   return;
	} else {
	   /* call LDU factorization: */
	    error= umfpack_di_numeric(A->pattern->ptr,A->pattern->index,A->val,pt->symbolic,&(pt->numeric),control,info);
	    if (error == UMFPACK_ERROR_out_of_memory ) {
		if (verbose) printf("UMFPACK: LDU factorization failed because of memory overlow.\n");
	      Esys_setError(MEMORY_ERROR,"UMFPACK: LDU factorization failed because of memory overlow.");
	      return;
	    } else if (error == UMFPACK_WARNING_singular_matrix ) {
	      if (verbose) printf("UMFPACK: LDU factorization failed because of singular matrix.\n");
	      Esys_setError(ZERO_DIVISION_ERROR,"UMFPACK: LDU factorization failed because of singular matrix.");
	      return;
	    } else if ( (error == UMFPACK_WARNING_determinant_underflow ) || ( error == UMFPACK_WARNING_determinant_overflow ) ) {
	      if (verbose) printf("UMFPACK: symbolic factorization failed because of under/overflow.\n");
	      Esys_setError(FLOATING_POINT_ERROR,"UMFPACK: symbolic factorization failed because of under/overflow.");
	      return;	
	    } else if (error != UMFPACK_OK) {
	      if (verbose) printf("UMFPACK: LDU factorization failed. Error code = %d.\n",error);
	      Esys_setError(SYSTEM_ERROR,"UMFPACK: factorization failed.");
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
	if (error == UMFPACK_ERROR_out_of_memory ) {
		if (verbose) printf("UMFPACK: forward/backward substitution failed because of memory overlow.\n");
	      Esys_setError(MEMORY_ERROR,"UMFPACK: forward/backward substitution failed because of memory overlow.");
	      return;
	} else if (error == UMFPACK_WARNING_singular_matrix ) {
	      if (verbose) printf("UMFPACK: forward/backward substitution because of singular matrix.\n");
	      Esys_setError(ZERO_DIVISION_ERROR,"UMFPACK: forward/backward substitution failed because of singular matrix.");
	      return;
	} else if ( (error == UMFPACK_WARNING_determinant_underflow ) || ( error == UMFPACK_WARNING_determinant_overflow ) ) {
	      if (verbose) printf("UMFPACK: forward/backward substitution failed because of under/overflow.\n");
	      Esys_setError(FLOATING_POINT_ERROR,"UMFPACK: forward/backward substitution failed because of under/overflow.");
	      return;	
	} else if (error != UMFPACK_OK) {
	   if (verbose) printf("UMFPACK: forward/backward substitution failed. Error code = %d.\n", error);
	   Esys_setError(SYSTEM_ERROR,"UMFPACK: forward/backward substitution failed.");
	   return;
	}
	if (verbose) printf("UMFPACK: forward/backward substitution completed (time = %e).\n",Esys_timer()-time0);
     }
     #else
	 Esys_setError(SYSTEM_ERROR,"Paso_UMFPACK: UMFPACK is not available.");
     #endif
}




