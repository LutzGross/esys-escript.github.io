
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to the Intel MKL library */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: l.gross@uq.edu.au */

/****************************************************************************/

#include "MKL.h"
#include "Options.h"
#include "PasoException.h"

namespace paso {

void MKL_free(SparseMatrix<double>* A)
{
#ifdef ESYS_HAVE_MKL
    if (A && A->solver_p && A->solver_package==PASO_MKL) {
        ES_MKL_INT mtype = MKL_MTYPE_REAL_UNSYM;
        ES_MKL_INT n = A->numRows;
        ES_MKL_INT maxfct=1; // number of factorizations on the same pattern
        ES_MKL_INT mnum =1;  // factorization to be handled in this call
        ES_MKL_INT msglvl=0; // message level
        ES_MKL_INT nrhs=1;   // number of right hand sides
        ES_MKL_INT idum;     // dummy integer
        _DOUBLE_PRECISION_t ddum;      // dummy float
        ES_MKL_INT error=MKL_ERROR_NO; // error code
        ES_MKL_INT iparm[64];          // parameters
        _MKL_DSS_HANDLE_t* pt = (_MKL_DSS_HANDLE_t*)A->solver_p;
        ES_MKL_INT phase = MKL_PHASE_RELEASE_MEMORY;
        ES_MKL_INT* ptr = reinterpret_cast<ES_MKL_INT*>(A->pattern->ptr);
        ES_MKL_INT* index = reinterpret_cast<ES_MKL_INT*>(A->pattern->index);
        for (index_t i=0; i<64; ++i)
            iparm[i]=0;

        ES_PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A->val, ptr,
                   index, &idum, &nrhs, iparm, &msglvl,&ddum, &ddum, &error);
        delete[] A->solver_p;
        A->solver_p=NULL;
        if (error != MKL_ERROR_NO)
            throw PasoException("Memory release in MKL library failed.");
    }
#endif
}

void MKL_solve(SparseMatrix_ptr<double> A, double* out, double* in, index_t reordering,
               dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_MKL
    if (! (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
        throw PasoException("Paso: MKL requires CSR format with index offset 1 and block size 1.");
    }

    // MKL uses 'long long int' in 64-bit version, escript 'long'. Make sure
    // they are compatible
    if (sizeof(ES_MKL_INT) != sizeof(index_t)) {
        throw PasoException("Paso: MKL index type is not compatible with this escript build. Check compile options.");
    }

    ES_MKL_INT* ptr = reinterpret_cast<ES_MKL_INT*>(A->pattern->ptr);
    ES_MKL_INT* index = reinterpret_cast<ES_MKL_INT*>(A->pattern->index);

    ES_MKL_INT mtype = MKL_MTYPE_REAL_UNSYM;
    ES_MKL_INT n = A->numRows;
    ES_MKL_INT maxfct = 1; /* number of factorizations on the same pattern */
    ES_MKL_INT mnum = 1; /* factorization to be handled in this call */
    ES_MKL_INT msglvl = (verbose? 1 : 0); /* message level */
    ES_MKL_INT nrhs = 1; /* number of right hand sides */
    ES_MKL_INT idum; /* dummy integer */
    _DOUBLE_PRECISION_t ddum; /* dummy float */
    ES_MKL_INT phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;
    ES_MKL_INT error=MKL_ERROR_NO;  /* error code */
    ES_MKL_INT iparm[64]; /* parameters */
    _MKL_DSS_HANDLE_t* pt = (_MKL_DSS_HANDLE_t *)(A->solver_p);

    for (index_t i=0; i<64; ++i)
        iparm[i]=0;
    iparm[0] = 1; // no default settings
    switch (reordering) {
        case PASO_MINIMUM_FILL_IN:
            iparm[1] = MKL_REORDERING_MINIMUM_DEGREE;
            break;
        default:
#ifdef _OPENMP
            iparm[1] = MKL_REORDERING_NESTED_DISSECTION_OMP;
#else
            iparm[1] = MKL_REORDERING_NESTED_DISSECTION;
#endif
            break;
    }
    iparm[5] = 0; // store solution into output array
    iparm[7] = numRefinements; // maximum number of refinements
    iparm[9] = 13; // 10**(-iparm[9]) perturbation of pivot elements
    iparm[10] = 1; // rescaling the matrix before factorization started
    iparm[17] = 0; // =-1 report number of non-zeroes
    iparm[18] = 0; // =-1 report flops
#ifdef _OPENMP
    if (omp_get_max_threads() > 8)
        iparm[23] = 1;
#endif

    double time0;

    if (pt==NULL) {
        // allocate address pointer
        pt = new _MKL_DSS_HANDLE_t[64];
        for (index_t i=0; i<64; ++i)
            pt[i] = NULL;
        A->solver_p = (void*) pt;
        A->solver_package = PASO_MKL;
        // symbolic factorization
        phase = MKL_PHASE_SYMBOLIC_FACTORIZATION;
        time0 = escript::gettime();
        ES_PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A->val,
                   ptr, index, &idum, &nrhs, iparm, &msglvl, in, out, &error);
        if (error != MKL_ERROR_NO) {
             if (verbose)
                 printf("MKL: symbolic factorization failed.\n");
             MKL_free(A.get());
             throw PasoException("symbolic factorization in MKL library failed.");
        } else {
            // LDU factorization
            phase = MKL_PHASE_FACTORIZATION;
            ES_PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A->val, ptr,
                       index, &idum, &nrhs, iparm, &msglvl, in, out, &error);
            if (error != MKL_ERROR_NO) {
                if (verbose)
                    printf("MKL: LDU factorization failed.\n");
                MKL_free(A.get());
                throw PasoException("factorization in MKL library failed. Most likely the matrix is singular.");
           }
           if (verbose)
               printf("MKL: LDU factorization completed (time = %e).\n", escript::gettime()-time0);
        }
    }
    // forward backward substitution
    time0 = escript::gettime();
    phase = MKL_PHASE_SOLVE;
    ES_PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A->val,
               ptr, index, &idum, &nrhs, iparm, &msglvl, in, out, &error);
    if (verbose) printf("MKL: solve completed.\n");
    if (error != MKL_ERROR_NO) {
        if (verbose)
            printf("MKL: forward/backward substitution failed.\n");
        throw PasoException("forward/backward substitution in MKL library failed. Most likely the matrix is singular.");
    } else {
        if (verbose)
            printf("MKL: forward/backward substitution completed (time = %e).\n", escript::gettime()-time0);
    }
#else
    throw PasoException("Paso: MKL is not available.");
#endif
}

void MKL_free(SparseMatrix<cplx_t>* A)
{
#ifdef ESYS_HAVE_MKL
    throw PasoException("Paso MKL_free(): complex not implemented.");
#else
    throw PasoException("Paso: MKL is not available.");
#endif
}

void MKL_solve(SparseMatrix_ptr<cplx_t> A, cplx_t* out, cplx_t* in, index_t reordering,
               dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_MKL
    throw PasoException("Paso MKL_solve(): complex not implemented.");
#else
    throw PasoException("Paso: MKL is not available.");
#endif
}

} // namespace paso

