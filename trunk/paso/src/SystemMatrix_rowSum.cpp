
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

/* Paso: SystemMatrix: calculates row sum                     */


/************************************************************************************/

/* Author: l.gross@auq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_rowSum(Paso_SystemMatrix* A, double* row_sum) {
   dim_t irow, nrow;
   if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1)) {
        Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_rowSum: No normalization available for compressed sparse column or index offset 1.");
      } else {
         nrow=A->mainBlock->numRows*A->row_block_size;
         #pragma omp parallel for private(irow) schedule(static)
         for (irow=0; irow<nrow ; ++irow) {
               row_sum[irow]=0.;
         }
         paso::SparseMatrix_addRow_CSR_OFFSET0(A->mainBlock,row_sum);
         paso::SparseMatrix_addRow_CSR_OFFSET0(A->col_coupleBlock,row_sum);
   }
}

