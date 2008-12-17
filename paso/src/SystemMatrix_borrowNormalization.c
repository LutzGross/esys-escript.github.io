
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix                                           */

/*  returns a borrowed reference to the matrix normaliztion vector */

/*  if the vector is in valid a new vector is calculate as the inverse */
/*  of the sum of the absolute value in each row/column                */

/**************************************************************/

/* Copyrights by ACcESS Australia 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

double* Paso_SystemMatrix_borrowNormalization(Paso_SystemMatrix* A) {
   dim_t irow, nrow;
   index_t irow_failed, irow_failed_local;
   register double fac;
   if (!A->normalizer_is_valid) {
      if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_SYM) || (A->type & MATRIX_FORMAT_OFFSET1)) {
        Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_borrowNormalization: No normalization available for compressed sparse column, symmetric storage scheme or index offset 1.");
      } else {
          if (Paso_checkPtr(A->normalizer)) {
              Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrix_borrowNormalization: no memory alloced for normalizer.");

          } else {
             nrow=A->mainBlock->numRows*A->row_block_size;
             #pragma omp parallel for private(irow) schedule(static)
             for (irow=0; irow<nrow ; ++irow) {
                  A->normalizer[irow]=0;
             }
             Paso_SparseMatrix_addAbsRow_CSR_OFFSET0(A->mainBlock,A->normalizer);
             Paso_SparseMatrix_addAbsRow_CSR_OFFSET0(A->col_coupleBlock,A->normalizer);
   
             #pragma omp parallel
             {
                irow_failed_local=-1;
                #pragma omp for private(irow,fac) schedule(static)
                for (irow=0; irow<nrow ; ++irow) {
                    fac=A->normalizer[irow];
                    if (ABS(fac)>0) {
                       A->normalizer[irow]=1./fac;
                    } else {
                       A->normalizer[irow]=1.;
                       irow_failed=irow;
                    }
                }
                #pragma omp critical
                irow_failed=irow_failed_local;
             }
             if (irow_failed>=0) {
                Paso_setError(ZERO_DIVISION_ERROR,"There is a row containing zero entries only.");
             }
             A->normalizer_is_valid=TRUE;
          }
      }
      Paso_MPIInfo_noError(A->mpi_info );
   }
   return A->normalizer;
}
