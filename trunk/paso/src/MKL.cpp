
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

/* Paso: interface to the Intel MKL library */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "MKL.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

/*  free any extra stuff possibly used by the MKL library */

void Paso_MKL_free(Paso_SparseMatrix* A) {
#ifdef MKL
      index_t i;
      if (A!=NULL) {

         if ((A->solver_p!=NULL) && (A->solver_package==PASO_MKL) ) {
             _INTEGER_t mtype = MKL_MTYPE_UNSYM;
             _INTEGER_t n = A->numRows;
             _INTEGER_t maxfct=1; /* number of factorizations on the same pattern */
             _INTEGER_t mnum =1; /* factorization to be handled in this call */
             _INTEGER_t msglvl=0; /* message level */
             _INTEGER_t nrhs=1; /* number of right hand sides */
             _INTEGER_t idum; /* dummy integer */
             _DOUBLE_PRECISION_t ddum; /* dummy float */
             _INTEGER_t error=MKL_ERROR_NO;  /* error code */  
             _INTEGER_t iparm[64]; /* parameters */  
             _MKL_DSS_HANDLE_t* pt = (_MKL_DSS_HANDLE_t *)(A->solver_p);
             _INTEGER_t phase = MKL_PHASE_RELEASE_MEMORY;  
              for (i=0;i<64;++i) iparm[i]=0;  
  
            PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                   &n, A->val, A->pattern->ptr, A->pattern->index, &idum, &nrhs,  
                   iparm, &msglvl,&ddum, &ddum, &error);  
              delete[] A->solver_p;  
              A->solver_p=NULL;
	      if (error != MKL_ERROR_NO) Esys_setError(TYPE_ERROR,"memory release in PARDISO library failed.");  
        }
     }
#endif
}

/*  call the solver: */

void Paso_MKL(Paso_SparseMatrix* A,
              double* out,
              double* in,
	      index_t reordering,
	      dim_t numRefinements,
	      bool verbose)
{	      

#ifdef MKL
     double time0;
     index_t i;

     _INTEGER_t mtype = MKL_MTYPE_UNSYM;
     _INTEGER_t n = A->numRows;
     _INTEGER_t maxfct=1; /* number of factorizations on the same pattern */
     _INTEGER_t mnum =1; /* factorization to be handled in this call */
     _INTEGER_t msglvl=0; /* message level */
     _INTEGER_t nrhs=1; /* number of right hand sides */
     _INTEGER_t idum; /* dummy integer */
     _DOUBLE_PRECISION_t ddum; /* dummy float */
     _INTEGER_t phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;

     _INTEGER_t error=MKL_ERROR_NO;  /* error code */
     _INTEGER_t iparm[64]; /* parameters */
     _MKL_DSS_HANDLE_t* pt = (_MKL_DSS_HANDLE_t *)(A->solver_p);
     
     if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        Esys_setError(TYPE_ERROR,"Paso_MKL: MKL requires CSR format with index offset 1 and block size 1.");
        return;
     }

     /* set iparm */
     for (i=0;i<64;++i) iparm[i]=0;
     iparm[0] = 1; /* no default settings*/
     switch (reordering) {
            case PASO_MINIMUM_FILL_IN:
               iparm[1]=MKL_REORDERING_MINIMUM_DEGREE;
               break;
            default:
                iparm[1]=MKL_REORDERING_NESTED_DISSECTION;
                break;
     }
     #ifdef _OPENMP
     iparm[2] =omp_get_max_threads();
     #else
     iparm[2] = 1;
     #endif
     iparm[5] = 0; /* store solution into output array */
     iparm[7] = numRefinements; /* maximum number of refinements */
     iparm[9] = 13; /* 10**(-iparm[9]) perturbation of pivot elements */
     iparm[10] = 1; /* rescaling the matrix before factorization started */
     iparm[17] =0; /* =-1 report number of non-zeroes */
     iparm[18] =0; /* =-1 report flops */

     if (pt==NULL) {
        /* allocate address pointer */
        pt=new _MKL_DSS_HANDLE_t[64];
        if (Esys_checkPtr(pt)) return;
        for (i=0;i<64;++i) pt[i]=NULL;
        A->solver_p=(void*) pt;
	A->solver_package=PASO_MKL;
        /* symbolic factorization */
        phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;
        time0=Esys_timer();
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, A->val, A->pattern->ptr, A->pattern->index, &idum, &nrhs,
                 iparm, &msglvl, in, out, &error);
        if (error != MKL_ERROR_NO) {
             if (verbose) printf("MKL: symbolic factorization failed.\n");
             Esys_setError(VALUE_ERROR,"symbolic factorization in PARDISO library failed.");
             Paso_MKL_free(A);
        } else {
           /* LDU factorization */
           phase = MKL_PHASE_FACTORIZATION;
           PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A->val, A->pattern->ptr, A->pattern->index, &idum, &nrhs,
                iparm, &msglvl, in, out, &error);
           if (error != MKL_ERROR_NO) {
             if (verbose) printf("MKL: LDU factorization failed.\n");
             Esys_setError(ZERO_DIVISION_ERROR,"factorization in PARDISO library failed. Most likely the matrix is singular.");
             Paso_MKL_free(A);
           }
           if (verbose) printf("MKL: LDU factorization completed (time = %e).\n",Esys_timer()-time0);
        }
     }
     /* forward backward substitution\ */
     if (Esys_noError())  {
        time0=Esys_timer();
        phase = MKL_PHASE_SOLVE;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, A->val, A->pattern->ptr, A->pattern->index, &idum, &nrhs,
                 iparm, &msglvl, in, out, &error);
        if (verbose) printf("MKL: solve completed.\n");
        if (error != MKL_ERROR_NO) {
              if (verbose) printf("MKL: forward/backward substitution failed.\n");
              Esys_setError(VALUE_ERROR,"forward/backward substitution in PARDISO library failed. Most likely the matrix is singular.");
        } else {
	   if (verbose) printf("MKL: forward/backward substitution completed (time = %e).\n",Esys_timer()-time0);
        }
     }
#else
    Esys_setError(SYSTEM_ERROR,"Paso_MKL: MKL is not available.");
#endif
}

