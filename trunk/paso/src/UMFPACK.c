
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

#define MAX_BLOCK_SIZE 3

/**************************************************************/

/*  free any extra stuff possibly used by the UMFPACK library */

void Paso_UMFPACK_free(Paso_SystemMatrix* A) {
     Paso_UMFPACK_Handler* pt =NULL;
      if (A->solver!=NULL) {
           pt=(Paso_UMFPACK_Handler*)(A->solver);
           Paso_UMFPACK1_free((Paso_UMFPACK_Handler*)A->solver);
           A->solver=NULL;
     }
}
void Paso_UMFPACK1_free(Paso_UMFPACK_Handler* pt) {
    if (pt!=NULL) {
#ifdef UMFPACK
         umfpack_di_free_symbolic(&(pt->symbolic));
         umfpack_di_free_numeric(&(pt->numeric));
#endif
         MEMFREE(pt);
    }
}


/*  call the solver: */

void Paso_UMFPACK(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options,
                          Paso_Performance* pp) {
#ifdef UMFPACK
     double time0;
     Paso_UMFPACK_Handler* pt = NULL;
     options->converged=FALSE;

     if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        Paso_setError(TYPE_ERROR,"Paso_UMFPACK: UMFPACK requires CSR format with index offset 1 and block size 1.");
        return;
     }
     Performance_startMonitor(pp,PERFORMANCE_ALL);
     pt = (Paso_UMFPACK_Handler *)(A->solver);

     time0=Paso_timer();
     Paso_UMFPACK1(&pt, A->mainBlock, out, in, 2);
     options->set_up_time=0;
     options->time=Paso_timer()-time0;
     if (!Paso_noError()) {
         Paso_UMFPACK_free(A);
     } else {
        if (options->verbose) printf("UMFPACK: solve completed.\n");
        A->solver=(void*) pt;
        options->converged=TRUE;
        options->residual_norm=0;
        options->num_iter=1;
        options->num_level=0;
        options->num_inner_iter=0;
     }
     Performance_stopMonitor(pp,PERFORMANCE_ALL);
}

void Paso_UMFPACK1(Paso_UMFPACK_Handler** pt, Paso_SparseMatrix* A, double* out, double* in, const int refines) {
     double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
     int error = UMFPACK_OK;
     dim_t i,j;
     /*****/
     double **xx;
     double **bb;
     Paso_SparseMatrix *block[MAX_BLOCK_SIZE];
    
     xx=MEMALLOC(A->row_block_size,double*);
     bb=MEMALLOC(A->row_block_size,double*);
     if (Paso_checkPtr(xx) || Paso_checkPtr(bb)) return;
     /****/

     umfpack_di_defaults(control);
     
     for (i=0;i<A->row_block_size;i++) {
          xx[i]=MEMALLOC(A->numRows,double);
          bb[i]=MEMALLOC(A->numRows,double);
          if (Paso_checkPtr(xx[i]) && Paso_checkPtr(bb[i])) return;
     }
     
     /*#pragma omp parallel for private(i,j) schedule(static)*/
     for (i=0;i<A->numRows;i++) {
         for (j=0;j<A->row_block_size;j++) {
          bb[j][i]=in[A->row_block_size*i+j];
          xx[j][i]=0;  
         }
      }
     
     for (i=0;i<MAX_BLOCK_SIZE;++i) {
          block[i]=NULL;
     }
     
     for (i=0;i<A->row_block_size;++i) {
        block[i]=Paso_SparseMatrix_getBlock(A,i+1);

          if (*pt==NULL) {
             int n = block[i]->numRows;
             *pt=(MEMALLOC(1,Paso_UMFPACK_Handler));
             if (Paso_checkPtr(*pt)) return;
             /* call LDU symbolic factorization: */
             error=umfpack_di_symbolic(n,n,block[i]->pattern->ptr,block[i]->pattern->index,block[i]->val,&((*pt)->symbolic),control,info);
             if (error != UMFPACK_OK) {
                  Paso_setError(VALUE_ERROR,"symbolic factorization failed.");
                  return;
             } else {
                 /* call LDU factorization: */
                 error= umfpack_di_numeric(block[i]->pattern->ptr,block[i]->pattern->index,block[i]->val,(*pt)->symbolic,&((*pt)->numeric),control,info);
                if (error != UMFPACK_OK) {
                  Paso_setError(ZERO_DIVISION_ERROR,"factorization failed. Most likely the matrix is singular.");
                  return;
                }
             }
          }
          
          if (Paso_noError())  {
             /* call forward backward substitution: */
             control[UMFPACK_IRSTEP]=refines; /* number of refinement steps */
             error=umfpack_di_solve(UMFPACK_A,block[i]->pattern->ptr,block[i]->pattern->index,block[i]->val,xx[i],bb[i],(*pt)->numeric,control,info);
             if (error != UMFPACK_OK) {
                   Paso_setError(VALUE_ERROR,"forward/backward substition failed. Most likely the matrix is singular.");
                   return;
             }
          }
     
     }
     
     /*#pragma omp parallel for private(i,j) schedule(static)*/
      for (i=0;i<A->numRows;i++) {
          for (j=0;j<A->row_block_size;j++) {
          out[A->row_block_size*i+j]=xx[j][i];
          }
       }
       
     for (i=0;i<A->row_block_size;i++) {
                MEMFREE(xx[i]);
                MEMFREE(bb[i]);
                Paso_SparseMatrix_free(block[i]);
      }
      
#else
    Paso_setError(SYSTEM_ERROR,"Paso_UMFPACK:UMFPACK is not avialble.");
#endif

}
