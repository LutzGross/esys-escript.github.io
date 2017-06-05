
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
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

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */
void SystemMatrix::MatrixVector(double alpha, const double* in, double beta,
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
                SparseMatrix_MatrixVector_CSR_OFFSET1(alpha, mainBlock, in, beta, out);
            }
        } else {
            MatrixVector_CSR_OFFSET0(alpha, in, beta, out);
        }
    }
}

void SystemMatrix::MatrixVector_CSR_OFFSET0(double alpha, const double* in,
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

} // namespace paso

