
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: matrix vector product with sparse matrix           */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/**************************************************************************/

/*  raw scaled vector update operation: out = alpha * A * in + beta * out */

/* has to be called within a parallel region                              */
/* barrier synconization is performed to make sure that the input vector available */

void  Paso_SystemMatrix_MatrixVector(double alpha,
                                     Paso_SystemMatrix* A,
                                     double* in,
                                     double beta,
                                     double* out) {

  double *snd_buffer=NULL, *rcv_buffer=NULL;
  dim_t N;
  Paso_MPIInfo *mpi_info=A->mpi_info;

  if (A->type & MATRIX_FORMAT_CSC) {
     if ( mpi_info->size>1) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSC is not supported by MPI.");
           return;
     } else {
       if (A->type & MATRIX_FORMAT_OFFSET1) {
         Paso_SparseMatrix_MatrixVector_CSC_OFFSET1(alpha,A->mainBlock,in,beta,out);
       } else {
         Paso_SparseMatrix_MatrixVector_CSC_OFFSET0(alpha,A->mainBlock,in,beta,out);
       }
     }
  } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: TRILINOS is not supported with MPI.");
           return;
  } else {
     if (A->type & MATRIX_FORMAT_OFFSET1) {
           if ( mpi_info->size>1) {
              Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_MatrixVector: CSR with index 1 is not supported by MPI.");
              return;
           } else {
              Paso_SparseMatrix_MatrixVector_CSR_OFFSET1(alpha,A->mainBlock,in,beta,out);
           }
     } else {
         Paso_SystemMatrix_allocBuffer(A);
         if (Paso_noError()) {
            Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(alpha,A,in,beta,out);
         }
         Paso_SystemMatrix_freeBuffer(A);
     }
  }
}

void  Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha,
                                                 Paso_SystemMatrix* A,
                                                 double* in,
                                                 double beta,
                                                 double* out)
{
  double *remote_values=NULL;
  /* start exchange */
  Paso_SystemMatrix_startCollect(A,in);
  /* process main block */
  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->mainBlock,in,beta,out);
  /* finish exchange */
  remote_values=Paso_SystemMatrix_finishCollect(A);
  /* process couple block */
  Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(alpha,A->coupleBlock,remote_values,1.,out);
}
