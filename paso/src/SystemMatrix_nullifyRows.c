
/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: SystemMatrix                                       */

/*  nullify rows in the matrix                             */

/*  the rows are marked by positive values in     */
/*  mask_row. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row       */
/*  are set to main_diagonal_value                   */


/**************************************************************/

/* Author: l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_nullifyRows(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value) {

  Paso_MPIInfo *mpi_info=A->mpi_info;
  if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Paso_noError()) {
            Paso_SparseMatrix_nullifyRows_CSR_BLK1(A->mainBlock,mask_row,main_diagonal_value);
            Paso_SparseMatrix_nullifyRows_CSR_BLK1(A->coupleBlock,mask_row,0.); 
         }
       }
  } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: TRILINOS is not supported with MPI.");
           return;
       } else {
         Paso_SystemMatrix_allocBuffer(A);
         if (Paso_noError()) {
            Paso_SparseMatrix_nullifyRows_CSR(A->mainBlock,mask_row,main_diagonal_value);
            Paso_SparseMatrix_nullifyRows_CSR(A->coupleBlock,mask_row,0.);
         }
         Paso_SystemMatrix_freeBuffer(A);
       }
  } 
  return;
}
