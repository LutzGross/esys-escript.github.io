
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __PASO_BLOCKOPS_H__
#define __PASO_BLOCKOPS_H__

#include "Paso.h"
#include "PasoException.h"

#include <cstring> // memcpy

#ifdef ESYS_HAVE_LAPACK
   #ifdef ESYS_MKL_LAPACK
      #include <mkl_lapack.h>
      #include <mkl_cblas.h>
   #else
      #include <lapacke.h>
      #include <cblas.h>
   #endif
#endif

namespace paso {

inline void BlockOps_Cpy_N(dim_t N, double* R, const double* V)
{
    memcpy((void*)R, (void*)V, N*sizeof(double));
}

/// performs operation R=R-mat*V (V and R are not overlapping) - 2x2
inline void BlockOps_SMV_2(double* R, const double* mat, const double* V)
{
    const double S1 = V[0];
    const double S2 = V[1];
    const double A11 = mat[0];
    const double A12 = mat[2];
    const double A21 = mat[1];
    const double A22 = mat[3];
    R[0] -= A11 * S1 + A12 * S2;
    R[1] -= A21 * S1 + A22 * S2;
}

/// performs operation R=R-mat*V (V and R are not overlapping) - 3x3
inline void BlockOps_SMV_3(double* R, const double* mat, const double* V)
{
    const double S1 = V[0];
    const double S2 = V[1];
    const double S3 = V[2];
    const double A11 = mat[0];
    const double A21 = mat[1];
    const double A31 = mat[2];
    const double A12 = mat[3];
    const double A22 = mat[4];
    const double A32 = mat[5];
    const double A13 = mat[6];
    const double A23 = mat[7];
    const double A33 = mat[8];
    R[0] -= A11 * S1 + A12 * S2 + A13 * S3;
    R[1] -= A21 * S1 + A22 * S2 + A23 * S3;
    R[2] -= A31 * S1 + A32 * S2 + A33 * S3;
}

#define PASO_MISSING_CLAPACK throw PasoException("You need to install a LAPACK version to enable operations on block sizes > 3.")

/// performs operation R=R-mat*V (V and R are not overlapping) - NxN
inline void BlockOps_SMV_N(dim_t N, double* R, const double* mat, const double* V)
{
#ifdef ESYS_HAVE_LAPACK
    cblas_dgemv(CblasColMajor,CblasNoTrans, N, N, -1., mat, N, V, 1, 1., R, 1);
#else
    PASO_MISSING_CLAPACK;
#endif
}

inline void BlockOps_MV_N(dim_t N, double* R, const double* mat, const double* V)
{
#ifdef ESYS_HAVE_LAPACK
    cblas_dgemv(CblasColMajor,CblasNoTrans, N, N, 1., mat, N, V, 1, 0., R, 1);
#else
    PASO_MISSING_CLAPACK;
#endif
}

inline void BlockOps_invM_2(double* invA, const double* A, int* failed)
{
    const double A11 = A[0];
    const double A12 = A[2];
    const double A21 = A[1];
    const double A22 = A[3];
    double D = A11*A22-A12*A21;
    if (std::abs(D) > 0) {
        D = 1./D;
        invA[0] =  A22*D;
        invA[1] = -A21*D;
        invA[2] = -A12*D;
        invA[3] =  A11*D;
    } else {
        *failed = 1;
    }
}

inline void BlockOps_invM_3(double* invA, const double* A, int* failed)
{
    const double A11 = A[0];
    const double A21 = A[1];
    const double A31 = A[2];
    const double A12 = A[3];
    const double A22 = A[4];
    const double A32 = A[5];
    const double A13 = A[6];
    const double A23 = A[7];
    const double A33 = A[8];
    double D = A11*(A22*A33-A23*A32) +
               A12*(A31*A23-A21*A33) +
               A13*(A21*A32-A31*A22);
    if (std::abs(D) > 0) {
        D = 1./D;
        invA[0] = (A22*A33-A23*A32)*D;
        invA[1] = (A31*A23-A21*A33)*D;
        invA[2] = (A21*A32-A31*A22)*D;
        invA[3] = (A13*A32-A12*A33)*D;
        invA[4] = (A11*A33-A31*A13)*D;
        invA[5] = (A12*A31-A11*A32)*D;
        invA[6] = (A12*A23-A13*A22)*D;
        invA[7] = (A13*A21-A11*A23)*D;
        invA[8] = (A11*A22-A12*A21)*D;
    } else {
        *failed = 1;
    }
}

/// LU factorization of NxN matrix mat with partial pivoting
inline void BlockOps_invM_N(dim_t N, double* mat, index_t* pivot, int* failed)
{
#ifdef ESYS_HAVE_LAPACK
#ifdef ESYS_MKL_LAPACK
    int res = 0;
    dgetrf(&N, &N, mat, &N, pivot, &res);
    if (res != 0)
        *failed = 1;
#else
    int res = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, mat, N, pivot);
    if (res != 0)
        *failed = 1;
#endif // ESYS_MKL_LAPACK
#else
    PASO_MISSING_CLAPACK;
#endif
}

/// solves system of linear equations A*X=B
inline void BlockOps_solve_N(dim_t N, double* X, double* mat, index_t* pivot, int* failed)
{
#ifdef ESYS_HAVE_LAPACK
#ifdef ESYS_MKL_LAPACK
    int res = 0;
    int ONE = 1;
    dgetrs("N", &N, &ONE, mat, &N, pivot, X, &N, &res);
    if (res != 0)
        *failed = 1;
#else
    int res = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', N, 1, mat, N, pivot, X, N);
    if (res != 0)
        *failed = 1;
#endif // ESYS_MKL_LAPACK
#else
    PASO_MISSING_CLAPACK;
#endif
}

/// inplace matrix vector product - order 2
inline void BlockOps_MViP_2(const double* mat, double* V)
{
    const double S1 = V[0];
    const double S2 = V[1];
    const double A11 = mat[0];
    const double A12 = mat[2];
    const double A21 = mat[1];
    const double A22 = mat[3];
    V[0] = A11 * S1 + A12 * S2;
    V[1] = A21 * S1 + A22 * S2;
}

/// inplace matrix vector product - order 3
inline void BlockOps_MViP_3(const double* mat, double* V)
{
    const double S1 = V[0];
    const double S2 = V[1];
    const double S3 = V[2];
    const double A11 = mat[0];
    const double A21 = mat[1];
    const double A31 = mat[2];
    const double A12 = mat[3];
    const double A22 = mat[4];
    const double A32 = mat[5];
    const double A13 = mat[6];
    const double A23 = mat[7];
    const double A33 = mat[8];
    V[0] = A11 * S1 + A12 * S2 + A13 * S3;
    V[1] = A21 * S1 + A22 * S2 + A23 * S3;
    V[2] = A31 * S1 + A32 * S2 + A33 * S3;
}

inline void BlockOps_solveAll(dim_t n_block, dim_t n, double* D,
                              index_t* pivot, double* x)
{
    if (n_block == 1) {
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i)
            x[i] *= D[i];
    } else if (n_block == 2) {
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i)
            BlockOps_MViP_2(&D[4*i], &x[2*i]);
    } else if (n_block == 3) {
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i)
            BlockOps_MViP_3(&D[9*i], &x[3*i]);
    } else {
        int failed = 0;
#pragma omp parallel for
        for (dim_t i=0; i<n; ++i) {
            const dim_t block_size = n_block*n_block;
            BlockOps_solve_N(n_block, &x[n_block*i], &D[block_size*i], &pivot[n_block*i], &failed);
        }
        if (failed > 0) {
            throw PasoException("BlockOps_solveAll: solution failed.");
        }
    }
}

} // namespace paso

#endif // __PASO_BLOCKOPS_H__

