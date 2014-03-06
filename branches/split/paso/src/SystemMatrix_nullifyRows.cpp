
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


/************************************************************************************/

/* Paso: SystemMatrix                                         */

/*  nullify rows in the matrix                                */

/*  The rows are marked by positive values in                 */
/*  mask_row. Values on the main diagonal                     */
/*  which are marked to set to zero by both mask_row          */
/*  are set to main_diagonal_value                            */


/************************************************************************************/

/* Author: l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_nullifyRows(Paso_SystemMatrix* A, double* mask_row, double main_diagonal_value) {
 
  double* remote_values;
  /*Esys_MPIInfo *mpi_info=A->mpi_info;*/
  if (A ->col_block_size==1 && A ->row_block_size ==1) {
       if (A->type & MATRIX_FORMAT_CSC) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Esys_noError()) {
            Paso_SystemMatrix_startRowCollect(A,mask_row);
            Paso_SparseMatrix_nullifyRows_CSR_BLK1(A->mainBlock,mask_row,main_diagonal_value);
            Paso_SparseMatrix_nullifyRows_CSR_BLK1(A->col_coupleBlock,mask_row,0.); 
            remote_values=Paso_SystemMatrix_finishRowCollect(A);
            Paso_SparseMatrix_nullifyRows_CSR_BLK1(A->row_coupleBlock,remote_values,0.);
         }
       }
  } else {
       if (A->type & MATRIX_FORMAT_CSC) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: CSC is not supported by MPI.");
           return;
       } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
           Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_nullifyRows: TRILINOS is not supported with MPI.");
           return;
       } else {
         if (Esys_noError()) {
            Paso_SystemMatrix_startRowCollect(A,mask_row);
            Paso_SparseMatrix_nullifyRows_CSR(A->mainBlock,mask_row,main_diagonal_value);
            Paso_SparseMatrix_nullifyRows_CSR(A->col_coupleBlock,mask_row,0.);
            remote_values=Paso_SystemMatrix_finishRowCollect(A);
            Paso_SparseMatrix_nullifyRows_CSR(A->row_coupleBlock,remote_values,0.);
         }
       }
  } 
}

