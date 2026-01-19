
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

/* Paso: matrix vector product with sparse matrix           */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SystemMatrix.h"

namespace paso {

template <>
void SystemMatrix<double>::MatrixVector_CSR_OFFSET0(double alpha, const double* in,
                                            double beta, double* out) const
{
    // start exchange
    startCollect(in);
    // process main block
    if (type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
        SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha, mainBlock, in, beta, out);
    } else {
        SparseMatrix_MatrixVector_CSR_OFFSET0(alpha, mainBlock, in, beta, out);
    }
    // finish exchange
    double* remote_values = finishCollect();
    // process couple block
    if (col_coupleBlock->pattern->ptr != NULL) {
        if (type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
            SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha, col_coupleBlock, remote_values, 1., out);
        } else {
            SparseMatrix_MatrixVector_CSR_OFFSET0(alpha, col_coupleBlock, remote_values, 1., out);
        }
    }
}

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */
template <>
void SystemMatrix<double>::MatrixVector(double alpha, const double* in, double beta,
                                double* out) const
{
    if (is_balanced) {
        throw PasoException("MatrixVector: balanced matrix is not supported.");
    }
    if (type & MATRIX_FORMAT_CSC) {
        if (mpi_info->size > 1) {
            throw PasoException("MatrixVector: CSC is not supported by MPI.");
        } else {
            if (type & MATRIX_FORMAT_OFFSET1) {
                SparseMatrix_MatrixVector_CSC_OFFSET1(alpha, mainBlock, in, beta, out);
            } else {
                SparseMatrix_MatrixVector_CSC_OFFSET0(alpha, mainBlock, in, beta, out);
            }
        }
    } else {
        if (type & MATRIX_FORMAT_OFFSET1) {
            if (mpi_info->size > 1) {
                throw PasoException("MatrixVector: CSR with offset 1 is not supported in MPI.");
            } else {
                SparseMatrix_MatrixVector_CSR_OFFSET1<double>(alpha, mainBlock, in, beta, out);
            }
        } else {
            MatrixVector_CSR_OFFSET0(alpha, in, beta, out);
        }
    }
}

template <>
void SystemMatrix<cplx_t>::MatrixVector(double alpha, const cplx_t* in, double beta,
                                cplx_t* out) const
{
#if defined(ESYS_HAVE_MUMPS)
    if (is_balanced) {
        throw PasoException("MatrixVector: balanced matrix is not supported.");
    }
    if (type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) {
        SparseMatrix_MatrixVector_CSR_OFFSET1<cplx_t>(alpha, mainBlock, in, beta, out);
        return;
    } else {
        throw PasoException("MatrixVector: MUMPS requires CSR format with "
                            "index offset 1 and block size 1.");
    }
#endif
    throw PasoException("MatrixVector: require MUMPS for complex matrices.");
}

} // namespace paso

