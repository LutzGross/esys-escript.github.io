
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

/* Paso: SystemMatrix                                         */

/*  nullify rows and columns in the matrix                    */
/*                                                            */
/*  The rows and columns are marked by positive values in     */
/*  mask_row and mask_col. Values on the main diagonal        */
/*  which are marked to set to zero by both mask_row and      */
/*  mask_col are set to main_diagonal_value.                  */


/****************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_nullifyRowsAndCols(Paso_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value)
{
  double *remote_values=NULL;
  Esys_MPIInfo *mpi_info=A->mpi_info;
  if (mpi_info->size>1) {
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Esys_noError()) {
            Paso_SystemMatrix_startColCollect(A, mask_col);
            Paso_SystemMatrix_startRowCollect(A, mask_row);
            A->mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row,mask_col,main_diagonal_value);
            remote_values=Paso_SystemMatrix_finishColCollect(A);
            A->col_coupleBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, remote_values, 0.);
            remote_values=Paso_SystemMatrix_finishRowCollect(A);
            A->row_coupleBlock->nullifyRowsAndCols_CSR_BLK1(remote_values, mask_col, 0.);
         }
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Esys_noError()) {
            Paso_SystemMatrix_startColCollect(A,mask_col) ;
            Paso_SystemMatrix_startRowCollect(A,mask_row) ;
            A->mainBlock->nullifyRowsAndCols_CSR(mask_row,mask_col,main_diagonal_value);
            remote_values=Paso_SystemMatrix_finishColCollect(A);
            A->col_coupleBlock->nullifyRowsAndCols_CSR(mask_row,remote_values,0.);
            remote_values=Paso_SystemMatrix_finishRowCollect(A);
            A->row_coupleBlock->nullifyRowsAndCols_CSR(remote_values,mask_col,0.); 
         }
       }
     } 
  } else { 
     if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           A->mainBlock->nullifyRowsAndCols_CSC_BLK1(mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported.");
       } else {
           A->mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row,mask_col,main_diagonal_value);
       }
     } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           A->mainBlock->nullifyRowsAndCols_CSC(mask_row,mask_col,main_diagonal_value);
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRowsAndCols: TRILINOS is not supported with MPI.");
       } else {
           A->mainBlock->nullifyRowsAndCols_CSR(mask_row,mask_col,main_diagonal_value);
       }
     }
  }
}

