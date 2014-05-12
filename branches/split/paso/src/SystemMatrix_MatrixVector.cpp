
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

/* Paso: matrix vector product with sparse matrix           */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

namespace paso {

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */
void SystemMatrix_MatrixVector(double alpha, SystemMatrix_ptr A,
                               const double* in, double beta, double* out)
{
    if (A->is_balanced) {
        Esys_setError(VALUE_ERROR, "SystemMatrix_MatrixVector: balanced matrix is not supported.");
        return;
    }
    if (A->type & MATRIX_FORMAT_CSC) {
        if (A->mpi_info->size > 1) {
            Esys_setError(SYSTEM_ERROR,"SystemMatrix_MatrixVector: CSC is not supported by MPI.");
            return;
        } else {
            if (A->type & MATRIX_FORMAT_OFFSET1) {
                SparseMatrix_MatrixVector_CSC_OFFSET1(alpha,A->mainBlock,in,beta,out);
            } else {
                SparseMatrix_MatrixVector_CSC_OFFSET0(alpha,A->mainBlock,in,beta,out);
            }
        }
    } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
        Esys_setError(SYSTEM_ERROR,"SystemMatrix_MatrixVector: TRILINOS is not supported with MPI.");
        return;
    } else {
        if (A->type & MATRIX_FORMAT_OFFSET1) {
            if (A->mpi_info->size > 1) {
                Esys_setError(SYSTEM_ERROR,"SystemMatrix_MatrixVector: CSR with offset 1 is not supported in MPI.");
                return;
            } else {
                SparseMatrix_MatrixVector_CSR_OFFSET1(alpha,A->mainBlock,in,beta,out);
            }
        } else {
            if (Esys_noError()) {
                SystemMatrix_MatrixVector_CSR_OFFSET0(alpha,A,in,beta,out);
            }
        }
    }
}

void SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, SystemMatrix_ptr A,
                                           const double* in, const double beta,
                                           double* out)
{
    // start exchange
    A->startCollect(in);
    // process main block
    if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
        SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha,A->mainBlock,in,beta,out);
    } else {
        SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->mainBlock,in,beta,out);
    }
    // finish exchange
    double* remote_values = A->finishCollect();
    // process couple block
    if (A->col_coupleBlock->pattern->ptr != NULL) {
        if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
            SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha,A->col_coupleBlock,remote_values,1.,out);
        } else {
            SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->col_coupleBlock,remote_values,1.,out);
        }
    }
}

} // namespace paso

