
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: matrix vector product with sparse matrix           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/**************************************************************************/

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */


void  Paso_SystemMatrix_MatrixVector(const double alpha,
                                     Paso_SystemMatrix* A,
                                     const double* in,
                                     const double beta,
                                     double* out) {

  /*double *snd_buffer=NULL, *rcv_buffer=NULL;*/
  Esys_MPIInfo *mpi_info=A->mpi_info;

  if (A->is_balanced) {
     Esys_setError(VALUE_ERROR,"Paso_SystemMatrix_MatrixVector: balanced matrix is not supported.");
     return;
  }
  if (A->type & MATRIX_FORMAT_CSC) {
     if ( mpi_info->size>1) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSC is not supported by MPI.");
           return;
     } else {
       if (A->type & MATRIX_FORMAT_OFFSET1) {
         Paso_SparseMatrix_MatrixVector_CSC_OFFSET1(alpha,A->mainBlock,in,beta,out);
       } else {
         Paso_SparseMatrix_MatrixVector_CSC_OFFSET0(alpha,A->mainBlock,in,beta,out);
       }
     }
  } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: TRILINOS is not supported with MPI.");
           return;
  } else {
     if (A->type & MATRIX_FORMAT_OFFSET1) {
           if ( mpi_info->size>1) {
              Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSR with index 1 is not supported by MPI.");
              return;
           } else {
              Paso_SparseMatrix_MatrixVector_CSR_OFFSET1(alpha,A->mainBlock,in,beta,out);
           }
     } else {
         if (Esys_noError()) {
            Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(alpha,A,in,beta,out);
         }
     }
  }
}

void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha,
                                                 Paso_SystemMatrix* A,
                                                 const double* in,
                                                 const double beta,
                                                 double* out)
{
  double *remote_values=NULL;
  /* start exchange */
  Paso_SystemMatrix_startCollect(A,in);
  /* process main block */
  if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha,A->mainBlock,in,beta,out);
  } else {
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->mainBlock,in,beta,out);
  }
  /* finish exchange */
  remote_values=Paso_SystemMatrix_finishCollect(A);
  /* process couple block */
  if (A->col_coupleBlock->pattern->ptr!=NULL) {
      if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(alpha,A->col_coupleBlock,remote_values,1.,out); 
      } else {
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->col_coupleBlock,remote_values,1.,out); 
      }
  }
}
