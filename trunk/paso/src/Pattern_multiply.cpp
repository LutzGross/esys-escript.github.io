
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

/* Paso: Pattern */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Pattern.h"

using esysUtils::IndexListArray;

namespace paso {

// computes the pattern coming from a matrix-matrix multiplication

Pattern* Pattern_multiply(int type, Pattern* A, Pattern* B)
{
    IndexListArray index_list(A->numOutput);

#pragma omp parallel for schedule(static)
    for (dim_t i = 0; i < A->numOutput; i++) {
        for (index_t iptrA = A->ptr[i]; iptrA < A->ptr[i+1]; ++iptrA) {
            const dim_t j = A->index[iptrA];
            for (index_t iptrB = B->ptr[j]; iptrB < B->ptr[j+1]; ++iptrB) {
                const dim_t k = B->index[iptrB];
                index_list[i].insertIndex(k);
            }
        }
    }
    Pattern* out=Pattern_fromIndexListArray(0, A->numOutput, index_list,
                                            0, B->numInput, 0);
    return out;
}

/*
 * Computes the pattern  of C = A binary operation B for CSR matrices A,B
 *
 * Note: we do not check whether A_ij(op)B_ij=0
 *
 */
Pattern* Pattern_binop(int type, Pattern* A, Pattern* B)
{
    IndexListArray index_list(A->numOutput);

#pragma omp parallel for schedule(static)
    for (dim_t i = 0; i < B->numOutput; i++) {
        index_t iptrA = A->ptr[i];
        index_t iptrB = B->ptr[i];

        while (iptrA < A->ptr[i+1] && iptrB < B->ptr[i+1]) {
            const dim_t j = A->index[iptrA];
            const dim_t k = B->index[iptrB];
            if (j < k) {
                index_list[i].insertIndex(j);
                iptrA++;
            } else if (j > k) {
                index_list[i].insertIndex(k);
                iptrB++;
            } else { // (j == k)
                index_list[i].insertIndex(j);
                iptrB++;
                iptrA++;
            }
        }
        while(iptrA < A->ptr[i+1]) {
            const dim_t j = A->index[iptrA];
            index_list[i].insertIndex(j);
            iptrA++;
        }
        while(iptrB < B->ptr[i+1]) {
            const dim_t k = B->index[iptrB];
            index_list[i].insertIndex(k);
            iptrB++;
        }
    }

    Pattern* out=Pattern_fromIndexListArray(0, A->numOutput, index_list, 0,
                                            A->numInput, 0);
    return out;
}

} // namespace paso

