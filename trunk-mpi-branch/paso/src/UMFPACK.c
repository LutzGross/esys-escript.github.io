/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: interface to the Intel UMFPACK library */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: gross@@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "UMFPACK.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

/*  free any extra stuff possibly used by the UMFPACK library */

void Paso_UMFPACK_free(Paso_SystemMatrix* A) {
     Paso_UMFPACK_Handler* pt =NULL;
#ifdef UMFPACK
      if (A->solver!=NULL) {
           pt=(Paso_UMFPACK_Handler*)(A->solver);
           umfpack_di_free_symbolic(&(pt->symbolic));
           umfpack_di_free_numeric(&(pt->numeric));
           MEMFREE(pt);
           A->solver=NULL;
     }
#endif
}
/*  call the solver: */

void Paso_UMFPACK(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options,
                          Paso_Performance* pp) {
#ifdef UMFPACK
     double time0;
     double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
     int error = UMFPACK_OK;
     Paso_UMFPACK_Handler* pt = NULL;

     if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        Paso_setError(TYPE_ERROR,"Paso_UMFPACK: UMFPACK requires CSR format with index offset 1 and block size 1.");
        return;
     }
     Performance_startMonitor(pp,PERFORMANCE_ALL);
     pt = (Paso_UMFPACK_Handler *)(A->solver);
     umfpack_di_defaults(control);

     if (pt==NULL) {
        int n = A->mainBlock->numRows;
        pt=MEMALLOC(1,Paso_UMFPACK_Handler);
        if (Paso_checkPtr(pt)) return;
        A->solver=(void*) pt;
        time0=Paso_timer();
        /* call LDU symbolic factorization: */
        error=umfpack_di_symbolic(n,n,A->mainBlock->pattern->ptr,A->mainBlock->pattern->index,A->mainBlock->val,&(pt->symbolic),control,info);
        if (error != UMFPACK_OK) {
             Paso_setError(VALUE_ERROR,"symbolic factorization failed.");
             Paso_UMFPACK_free(A);
        } else {
            /* call LDU factorization: */
            error= umfpack_di_numeric(A->mainBlock->pattern->ptr,A->mainBlock->pattern->index,A->mainBlock->val,pt->symbolic,&(pt->numeric),control,info);
           if (error != UMFPACK_OK) {
             Paso_setError(ZERO_DIVISION_ERROR,"factorization failed. Most likely the matrix is singular.");
             Paso_UMFPACK_free(A);
           }
           if (options->verbose) printf("timing UMFPACK: LDU factorization: %.4e sec.\n",Paso_timer()-time0);
        }
     }
     if (Paso_noError())  {
        time0=Paso_timer();
        /* call forward backward substitution: */
        control[UMFPACK_IRSTEP]=2; /* number of refinement steps */
        error=umfpack_di_solve(UMFPACK_A,A->mainBlock->pattern->ptr,A->mainBlock->pattern->index,A->mainBlock->val,out,in,pt->numeric,control,info);
        if (options->verbose) printf("timing UMFPACK: solve: %.4e sec\n",Paso_timer()-time0);
        if (error != UMFPACK_OK) {
              Paso_setError(VALUE_ERROR,"forward/backward substition failed. Most likely the matrix is singular.");
        }
     }
     Performance_stopMonitor(pp,PERFORMANCE_ALL);
#else
    Paso_setError(SYSTEM_ERROR,"Paso_UMFPACK:UMFPACK is not avialble.");
#endif

}
