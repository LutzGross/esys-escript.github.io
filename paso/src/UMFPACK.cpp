
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to the UMFPACK library */

/****************************************************************************/

#include "Paso.h"
#include "UMFPACK.h"
#include "Options.h"
#include "PasoException.h"

#include <iostream>
#include <sstream>

namespace paso {

/// frees any UMFPACK related data from the matrix
void UMFPACK_free(SparseMatrix<double>* A)
{
    if (A && A->solver_p) {
        UMFPACK_Handler* pt = reinterpret_cast<UMFPACK_Handler*>(A->solver_p);
#ifdef ESYS_HAVE_UMFPACK
#ifdef ESYS_INDEXTYPE_LONG
        umfpack_dl_free_symbolic(&pt->symbolic);
        umfpack_dl_free_numeric(&pt->numeric);
#else
        umfpack_di_free_symbolic(&pt->symbolic);
        umfpack_di_free_numeric(&pt->numeric);
#endif // ESYS_INDEXTYPE_LONG
#endif
        delete pt;
        A->solver_p = NULL;
    }
}


/// calls the solver
void UMFPACK_solve(SparseMatrix_ptr<double> A, double* out, double* in,
                   dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_UMFPACK
    if (!( (A->type & MATRIX_FORMAT_BLK1) && (A->type & MATRIX_FORMAT_CSC)) ) {
        throw PasoException("Paso: UMFPACK requires CSC format with index offset 1 and block size 1.");
    }

    UMFPACK_Handler* pt = reinterpret_cast<UMFPACK_Handler*>(A->solver_p);
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
#ifdef ESYS_INDEXTYPE_LONG
    umfpack_dl_defaults(control);
#else
    umfpack_di_defaults(control);
#endif
#ifdef UMFPACK_ORDERING_METIS
    control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS; // in versions > 5.2
#endif
    double time0;
    int error;
    if (pt == NULL) {
        int n = A->numRows;
        pt = new UMFPACK_Handler;
        A->solver_p = (void*) pt;
        A->solver_package = PASO_UMFPACK;
        time0=escript::gettime();

        // call LDU symbolic factorization:
#ifdef ESYS_INDEXTYPE_LONG
        error = umfpack_dl_symbolic(n, n, A->pattern->ptr,
                                    A->pattern->index, A->val,
                                    &pt->symbolic, control, info);
#else
        error = umfpack_di_symbolic(n, n, A->pattern->ptr,
                                    A->pattern->index, A->val,
                                    &pt->symbolic, control, info);
#endif
        if (error != UMFPACK_OK) {
            std::string message;
            if (error == UMFPACK_ERROR_out_of_memory) {
                message = "UMFPACK: symbolic factorization failed because of "
                          "memory overflow.";
            } else if (error == UMFPACK_WARNING_singular_matrix) {
                message = "UMFPACK: symbolic factorization failed because of "
                          "singular matrix.";
            } else if (error == UMFPACK_WARNING_determinant_underflow ||
                       error == UMFPACK_WARNING_determinant_overflow) {
                message = "UMFPACK: symbolic factorization failed because of "
                          "under/overflow.";
            } else {
                std::stringstream ss;
                ss << "UMFPACK: symbolic factorization failed. UMFPACK "
                      "error code = " << error << ".";
                message = ss.str();
            }
            if (verbose)
                std::cout << message.c_str() << std::endl;
            throw PasoException(message);
        }

        // call LDU factorization:
#ifdef ESYS_INDEXTYPE_LONG
        error = umfpack_dl_numeric(A->pattern->ptr, A->pattern->index,
                                   A->val, pt->symbolic, &pt->numeric,
                                   control, info);
#else
        error = umfpack_di_numeric(A->pattern->ptr, A->pattern->index,
                                   A->val, pt->symbolic, &pt->numeric,
                                   control, info);
#endif
        if (error == UMFPACK_OK) {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization completed (time = "
                    << escript::gettime()-time0 << ")." << std::endl;
            }
        } else if (error == UMFPACK_ERROR_out_of_memory) {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed because of "
                    "memory overflow." << std::endl;
            }
            throw PasoException("UMFPACK: LDU factorization failed because of memory overflow.");
        } else if (error == UMFPACK_WARNING_singular_matrix) {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed because of "
                    "singular matrix." << std::endl;
            }
            throw PasoException("UMFPACK: LDU factorization failed because of singular matrix.");
        } else if (error == UMFPACK_WARNING_determinant_underflow
                   || error == UMFPACK_WARNING_determinant_overflow) {
            if (verbose) {
                std::cout << "UMFPACK: symbolic factorization failed because "
                    "of under/overflow." << std::endl;
            }
            throw PasoException("UMFPACK: symbolic factorization failed because of under/overflow.");
        } else {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed. UMFPACK "
                    "error code = " << error << "." << std::endl;
            }
            throw PasoException("UMFPACK: factorization failed.");
        }
    } // pt==NULL

    // call forward backward substitution
    control[UMFPACK_IRSTEP] = numRefinements; // number of refinement steps
    time0 = escript::gettime();
#ifdef ESYS_INDEXTYPE_LONG
    error = umfpack_dl_solve(UMFPACK_A, A->pattern->ptr, A->pattern->index,
                             A->val, out, in, pt->numeric, control, info);
#else
    error = umfpack_di_solve(UMFPACK_A, A->pattern->ptr, A->pattern->index,
                             A->val, out, in, pt->numeric, control, info);
#endif

    if (error == UMFPACK_OK) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution completed "
                "(time = " << escript::gettime()-time0 << ")." << std::endl;
        }
    } else if (error == UMFPACK_ERROR_out_of_memory) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed "
                "because of memory overflow." << std::endl;
        }
        throw PasoException("UMFPACK: forward/backward substitution failed because of memory overflow.");
    } else if (error == UMFPACK_WARNING_singular_matrix) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution because of "
                "singular matrix." << std::endl;
        }
        throw PasoException("UMFPACK: forward/backward substitution failed because of singular matrix.");
    } else if (error == UMFPACK_WARNING_determinant_underflow
                 || error == UMFPACK_WARNING_determinant_overflow) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed "
                "because of under/overflow." << std::endl;
        }
        throw PasoException("UMFPACK: forward/backward substitution failed because of under/overflow.");
    } else {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed. "
                "UMFPACK error code = " << error << "." << std::endl;
        }
        throw PasoException("UMFPACK: forward/backward substitution failed.");
    }
#else // ESYS_HAVE_UMFPACK
    throw PasoException("Paso: Not compiled with UMFPACK.");
#endif
}

/// frees any UMFPACK related data from the matrix
void UMFPACK_free(SparseMatrix<cplx_t>* A)
{
#ifdef ESYS_HAVE_UMFPACK
    throw PasoException("Paso UMFPACK_free(): complex not implemented.");
#else
    throw PasoException("Paso: Not compiled with UMFPACK.");
#endif
}

/// calls the solver
void UMFPACK_solve(SparseMatrix_ptr<cplx_t> A, cplx_t* out, cplx_t* in,
                   dim_t numRefinements, bool verbose)
{
#ifdef ESYS_HAVE_UMFPACK
    throw PasoException("Paso UMFPACK_solve(): complex not implemented.");
#else
    throw PasoException("Paso: Not compiled with UMFPACK.");
#endif
}

} // namespace paso

