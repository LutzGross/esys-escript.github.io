
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

/* Paso: SystemMatrix                                       */

/*  nullify rows and columns in the matrix                    */
/*
/*  the rows and columns are marked by positive values in     */
/*  mask_row and mask_col. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row and      */
/*  mask_col are set to main_diagonal_value                   */


/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value) {

  double *remote_values=NULL;
  Paso_MPIInfo *mpi_info=A->mpi_info;
  if (mpi_info->size>1) {
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Paso_noError()) {
            Paso_SystemMatrix_startColCollect(A,mask_col);
            Paso_SystemMatrix_startRowCollect(A,mask_row);
            Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(A->mainBlock,mask_row,mask_col,main_diagonal_value);
            remote_values=Paso_SystemMatrix_finishColCollect(A);
            Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(A->col_coupleBlock,mask_row,remote_values,0.); 
            /* remote_values=Paso_SystemMatrix_finishRowCollect(A);
            Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(A->row_coupleBlock,remote_values,mask_col,0.); 
            */
         }
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Paso_noError()) {
            Paso_SystemMatrix_startColCollect(A,mask_col) ;
            Paso_SystemMatrix_startRowCollect(A,mask_row) ;
            Paso_SparseMatrix_nullifyRowsAndCols_CSR(A->mainBlock,mask_row,mask_col,main_diagonal_value);
            remote_values=Paso_SystemMatrix_finishColCollect(A);
            Paso_SparseMatrix_nullifyRowsAndCols_CSR(A->col_coupleBlock,mask_row,remote_values,0.);
            /*
            remote_values=Paso_SystemMatrix_finishRowCollect(A);
            Paso_SparseMatrix_nullifyRowsAndCols_CSR(A->row_coupleBlock,remote_values,mask_col,0.); 
            */
         }
       }
     } 
  } else { 
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_SparseMatrix_nullifyRowsAndCols_CSC_BLK1(A->mainBlock,mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported.");
       } else {
           Paso_SparseMatrix_nullifyRowsAndCols_CSR_BLK1(A->mainBlock,mask_row,mask_col,main_diagonal_value);
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Paso_SparseMatrix_nullifyRowsAndCols_CSC(A->mainBlock,mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
       } else {
           Paso_SparseMatrix_nullifyRowsAndCols_CSR(A->mainBlock,mask_row,mask_col,main_diagonal_value);
       }
     }
  }
  return;
}
