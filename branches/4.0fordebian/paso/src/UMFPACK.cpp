
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: interface to the UMFPACK library */

/****************************************************************************/

#include "UMFPACK.h"
#include "Paso.h"
#include "Options.h"

#include <iostream>
#include <sstream>

namespace paso {

/// frees any UMFPACK related data from the matrix
void UMFPACK_free(SparseMatrix* A)
{
    if (A && A->solver_p) {
        UMFPACK_Handler* pt = reinterpret_cast<UMFPACK_Handler*>(A->solver_p);
#ifdef USE_UMFPACK
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
void UMFPACK_solve(SparseMatrix_ptr A, double* out, double* in,
                   dim_t numRefinements, bool verbose)
{
#ifdef USE_UMFPACK
    if (!( (A->type & MATRIX_FORMAT_BLK1) && (A->type & MATRIX_FORMAT_CSC)) ) {
        Esys_setError(TYPE_ERROR, "Paso: UMFPACK requires CSC format with index offset 1 and block size 1.");
        return;
    }

    UMFPACK_Handler* pt = reinterpret_cast<UMFPACK_Handler*>(A->solver_p);
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    umfpack_di_defaults(control);
    double time0;
    int error;

    if (pt == NULL) {
        int n = A->numRows;
        pt = new UMFPACK_Handler;
        A->solver_p = (void*) pt;
        A->solver_package = PASO_UMFPACK;
        time0=Esys_timer();

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
                Esys_setError(MEMORY_ERROR, message.c_str());
            } else if (error == UMFPACK_WARNING_singular_matrix) {
                message = "UMFPACK: symbolic factorization failed because of "
                          "singular matrix.";
                Esys_setError(ZERO_DIVISION_ERROR, message.c_str());
            } else if (error == UMFPACK_WARNING_determinant_underflow ||
                       error == UMFPACK_WARNING_determinant_overflow) {
                message = "UMFPACK: symbolic factorization failed because of "
                          "under/overflow.";
                Esys_setError(FLOATING_POINT_ERROR, message.c_str());
            } else {
                std::stringstream ss;
                ss << "UMFPACK: symbolic factorization failed. UMFPACK "
                      "error code = " << error << ".";
                message = ss.str();
                Esys_setError(SYSTEM_ERROR, message.c_str());
            }
            if (verbose)
                std::cout << message.c_str() << std::endl;
            return;
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
                    << Esys_timer()-time0 << ")." << std::endl;
            }
        } else if (error == UMFPACK_ERROR_out_of_memory) {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed because of "
                    "memory overflow." << std::endl;
            }
            Esys_setError(MEMORY_ERROR, "UMFPACK: LDU factorization failed because of memory overflow.");
            return;
        } else if (error == UMFPACK_WARNING_singular_matrix) {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed because of "
                    "singular matrix." << std::endl;
            }
            Esys_setError(ZERO_DIVISION_ERROR,"UMFPACK: LDU factorization failed because of singular matrix.");
            return;
        } else if (error == UMFPACK_WARNING_determinant_underflow
                   || error == UMFPACK_WARNING_determinant_overflow) {
            if (verbose) {
                std::cout << "UMFPACK: symbolic factorization failed because "
                    "of under/overflow." << std::endl;
            }
            Esys_setError(FLOATING_POINT_ERROR,"UMFPACK: symbolic factorization failed because of under/overflow.");
            return;
        } else {
            if (verbose) {
                std::cout << "UMFPACK: LDU factorization failed. UMFPACK "
                    "error code = " << error << "." << std::endl;
            }
            Esys_setError(SYSTEM_ERROR, "UMFPACK: factorization failed.");
            return;
        }
    } // pt==NULL

    // call forward backward substitution
    control[UMFPACK_IRSTEP] = numRefinements; // number of refinement steps
    time0 = Esys_timer();
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
                "(time = " << Esys_timer()-time0 << ")." << std::endl;
        }
    } else if (error == UMFPACK_ERROR_out_of_memory) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed "
                "because of memory overflow." << std::endl;
        }
        Esys_setError(MEMORY_ERROR, "UMFPACK: forward/backward substitution failed because of memory overflow.");
    } else if (error == UMFPACK_WARNING_singular_matrix) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution because of "
                "singular matrix." << std::endl;
        }
        Esys_setError(ZERO_DIVISION_ERROR, "UMFPACK: forward/backward substitution failed because of singular matrix.");
    } else if (error == UMFPACK_WARNING_determinant_underflow
                 || error == UMFPACK_WARNING_determinant_overflow) {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed "
                "because of under/overflow." << std::endl;
        }
        Esys_setError(FLOATING_POINT_ERROR, "UMFPACK: forward/backward substitution failed because of under/overflow.");
    } else {
        if (verbose) {
            std::cout << "UMFPACK: forward/backward substitution failed. "
                "UMFPACK error code = " << error << "." << std::endl;
        }
        Esys_setError(SYSTEM_ERROR, "UMFPACK: forward/backward substitution failed.");
    }
#else // USE_UMFPACK
    Esys_setError(SYSTEM_ERROR, "Paso: Not compiled with UMFPACK.");
#endif
}

} // namespace paso

