
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

/* Paso: interface to the Intel MKL library */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: gross@@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "MKL.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

/*  free any extra stuff possibly used by the MKL library */

void Paso_MKL_free(Paso_SystemMatrix* A) {
#ifdef MKL
      index_t i;
      if (A->solver!=NULL) {
           _INTEGER_t mtype = MKL_MTYPE_UNSYM;
           if (A->type & MATRIX_FORMAT_SYM) mtype=MKL_MTYPE_SYM;
          _INTEGER_t n = A->mainBlock->numRows;
          _INTEGER_t maxfct=1; /* number of factorizations on the same pattern */
          _INTEGER_t mnum =1; /* factoriztion to be handeled in this call */
          _INTEGER_t msglvl=0; /* message level */
          _INTEGER_t nrhs=1; /* number of right hand sides */
          _INTEGER_t idum; /* dummy integer */
          _DOUBLE_PRECISION_t ddum; /* dummy float */
          _INTEGER_t error=MKL_ERROR_NO;  /* error code */
          _INTEGER_t iparm[64]; /* parameters */
          for (i=0;i<64;++i) iparm[i]=0;

          _INTEGER_t phase = MKL_PHASE_RELEASE_MEMORY;
          PARDISO ((_MKL_DSS_HANDLE_t *)(A->solver), &maxfct, &mnum, &mtype, &phase,
                   &n, A->mainBlock->val, A->mainBlock->pattern->ptr, A->mainBlock->pattern->index, &idum, &nrhs,
                   iparm, &msglvl,&ddum, &ddum, &error);
          MEMFREE(A->solver);
          if (error != MKL_ERROR_NO) Paso_setError(TYPE_ERROR,"memory release in paradiso library failed.");
     }
#endif
}
/*  call the solver: */

void Paso_MKL(Paso_SystemMatrix* A,
                          double* out,
                          double* in,
                          Paso_Options* options,
                          Paso_Performance* pp) {
#ifdef MKL
     double time0;
     index_t i;

     if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        Paso_setError(TYPE_ERROR,"Paso_MKL: MKL requires CSR format with index offset 1 and block size 1.");
        return;
     }
     Performance_startMonitor(pp,PERFORMANCE_ALL);
     _INTEGER_t mtype = MKL_MTYPE_UNSYM;
     if (A->type & MATRIX_FORMAT_SYM) mtype=MKL_MTYPE_SYM;
     _INTEGER_t n = A->mainBlock->numRows;
     _INTEGER_t maxfct=1; /* number of factorizations on the same pattern */
     _INTEGER_t mnum =1; /* factoriztion to be handeled in this call */
     _INTEGER_t msglvl=0; /* message level */
     _INTEGER_t nrhs=1; /* number of right hand sides */
     _INTEGER_t idum; /* dummy integer */
     _DOUBLE_PRECISION_t ddum; /* dummy float */
     _INTEGER_t phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;

     _INTEGER_t error=MKL_ERROR_NO;  /* error code */
     _INTEGER_t iparm[64]; /* parameters */
     _MKL_DSS_HANDLE_t* pt = (_MKL_DSS_HANDLE_t *)(A->solver);
     /* set iparm */
     for (i=0;i<64;++i) iparm[i]=0;
     iparm[0] = 1; /* no default settings*/
     switch (options->reordering) {
            case PASO_MINIMUM_FILL_IN:
               iparm[1]=MKL_REORDERING_MINIMUM_DEGREE;
               break;
            default:
                iparm[1]=MKL_REORDERING_NESTED_DISSECTION;
                break;
     }
     #ifdef _OPENMP
     iparm[2] = omp_get_max_threads();
     #else
     iparm[2] = 1;
     #endif
     iparm[5] = 0; /* store solution into output array */
     iparm[7] = 2; /* maximum number of refinements */
     iparm[9] = 13; /* 10**(-iparm[9]) preturbation of pivot elements */
     iparm[10] = 1; /* rescaling the matrix before factorization started */
     iparm[17] =0; /* =-1 report number of non-zeroes */
     iparm[18] =0; /* =-1 report flops */


     if (pt==NULL) {
        /* allocate address pointer */
        pt=MEMALLOC(64,_MKL_DSS_HANDLE_t);
        if (Paso_checkPtr(pt)) return;
        A->solver=(void*) pt;
        for (i=0;i<64;++i) pt[i]=NULL;
        time0=Paso_timer();
        /* symbolic factorization */
        phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, A->mainBlock->val, A->mainBlock->pattern->ptr, A->mainBlock->pattern->index, &idum, &nrhs,
                 iparm, &msglvl, in, out, &error);
        if (error != MKL_ERROR_NO) {
             Paso_setError(VALUE_ERROR,"symbolic factorization in paradiso library failed.");
             Paso_MKL_free(A);
        } else {
           /* LDU factorization */
           phase = MKL_PHASE_FACTORIZATION;
           PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A->mainBlock->val, A->mainBlock->pattern->ptr, A->mainBlock->pattern->index, &idum, &nrhs,
                iparm, &msglvl, in, out, &error);
           if (error != MKL_ERROR_NO) {
             Paso_setError(ZERO_DIVISION_ERROR,"factorization in paradiso library failed. Most likely the matrix is singular.");
             Paso_MKL_free(A);
           }
           if (options->verbose) printf("timing MKL: LDU factorization: %.4e sec.\n",Paso_timer()-time0);
        }
     }
     /* forward backward substitution\ */
     if (Paso_noError())  {
        time0=Paso_timer();
        phase = MKL_PHASE_SOLVE;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, A->mainBlock->val, A->mainBlock->pattern->ptr, A->mainBlock->pattern->index, &idum, &nrhs,
                 iparm, &msglvl, in, out, &error);
        if (options->verbose) printf("timing MKL: solve: %.4e sec\n",Paso_timer()-time0);
        if (error != MKL_ERROR_NO) {
              Paso_setError(VALUE_ERROR,"forward/backward substition in paradiso library failed. Most likely the matrix is singular.");
        }
     }
     Performance_stopMonitor(pp,PERFORMANCE_ALL);
#else
    Paso_setError(SYSTEM_ERROR,"Paso_MKL:MKL is not avialble.");
#endif
}
/*
 * $Log$
 *
 */
