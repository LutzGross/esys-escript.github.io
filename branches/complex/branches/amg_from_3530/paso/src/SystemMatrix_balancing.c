
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix                                         */

/*  Returns a borrowed reference to the matrix normalization vector. */

/*  If the vector is in valid a new vector is calculated as the  */
/*  inverse of the sum of the absolute value in each row/column. */

/**************************************************************/

/* Copyrights by ACcESS Australia 2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_applyBalanceInPlace(const Paso_SystemMatrix* A, double* x, const bool_t RHS)
{
   const dim_t nrow=Paso_SystemMatrix_getTotalNumRows(A);
   const dim_t ncol=Paso_SystemMatrix_getTotalNumCols(A);
   index_t i;
   
   if (A->is_balanced) {
      if (RHS) {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0; i<nrow ; ++i) {
	    x[i]*=A->balance_vector[i];
	 }
      } else {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0; i<ncol ; ++i) {
	    x[i]*=A->balance_vector[i];
	 }
      }
   }
   return;
}

void Paso_SystemMatrix_applyBalance(const Paso_SystemMatrix* A, double* x_out, const double* x, const bool_t RHS)
{
   const dim_t nrow=Paso_SystemMatrix_getTotalNumRows(A);
   const dim_t ncol=Paso_SystemMatrix_getTotalNumCols(A);
   index_t i;
   
   if (A->is_balanced) {
      if (RHS) {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0; i<nrow ; ++i) {
	    x_out[i]=x[i]*A->balance_vector[i];
	 }
      } else {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0; i<ncol ; ++i) {
	    x_out[i]=x[i]*A->balance_vector[i];
	 }
      }
   }
   return;
}

void Paso_SystemMatrix_balance(Paso_SystemMatrix* A) {
   dim_t irow;
   const dim_t nrow=Paso_SystemMatrix_getTotalNumRows(A);
   register double fac;

   if (!A->is_balanced) {
      if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_OFFSET1)) {
        Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_balance: No normalization available for compressed sparse column or index offset 1.");
      } 
      if (Esys_checkPtr(A->balance_vector)) {
	 Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_balance: no memory allocated for balance vector.");
      } 
      if ( ! ( (Paso_SystemMatrix_getGlobalNumRows(A) == Paso_SystemMatrix_getGlobalNumCols(A)) && (A->row_block_size == A->col_block_size) ) ) {
	 Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_balance: matrix needs to be a square matrix.");
      }
      if (Esys_noError() ) {

             /* calculate absolute max value over each row */
             #pragma omp parallel for private(irow) schedule(static)
             for (irow=0; irow<nrow ; ++irow) {
		A->balance_vector[irow]=0;
             }
             Paso_SparseMatrix_maxAbsRow_CSR_OFFSET0(A->mainBlock,A->balance_vector);
             if(A->col_coupleBlock->pattern->ptr!=NULL) {
		Paso_SparseMatrix_maxAbsRow_CSR_OFFSET0(A->col_coupleBlock,A->balance_vector);  
             }

             /* set balancing vector */
             #pragma omp parallel
             {
                #pragma omp for private(irow,fac) schedule(static)
                for (irow=0; irow<nrow ; ++irow) {
		    fac=A->balance_vector[irow];
                    if (fac>0) {
		       A->balance_vector[irow]=sqrt(1./fac);
                    } else {
		       A->balance_vector[irow]=1.;
                    }
                }
	    }
            {
		  /* rescale matrix: */
		  double *remote_values=NULL;
		  /* start exchange */
		  Paso_SystemMatrix_startCollect(A, A->balance_vector);
		  /* process main block */
		  Paso_SparseMatrix_applyDiagonal_CSR_OFFSET0(A->mainBlock,A->balance_vector, A->balance_vector);
		  /* finish exchange */
		  remote_values=Paso_SystemMatrix_finishCollect(A);
		  /* process couple block */
		  if (A->col_coupleBlock->pattern->ptr!=NULL) {
		     Paso_SparseMatrix_applyDiagonal_CSR_OFFSET0(A->col_coupleBlock,A->balance_vector, remote_values); 
		  }
		  if (A->row_coupleBlock->pattern->ptr!=NULL) {
		     Paso_SparseMatrix_applyDiagonal_CSR_OFFSET0(A->row_coupleBlock, remote_values, A->balance_vector); 
		  }
	     }
             A->is_balanced=TRUE;
          }
      }
}
