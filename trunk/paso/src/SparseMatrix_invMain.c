
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/* Paso: inverts the main diagonal of a SparseMatrix: */

/************************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "Solver.h"
#include "BlockOps.h"
#include "PasoUtil.h"

void Paso_SparseMatrix_invMain(Paso_SparseMatrix * A_p, double* inv_diag, int* pivot) {
   index_t failed=0;
   register double A11;
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   const dim_t m_block=A_p->col_block_size;
   const dim_t block_size=A_p->block_size;
   dim_t i;
   register index_t iPtr;
   index_t* main_ptr=Paso_Pattern_borrowMainDiagonalPointer(A_p->pattern);
   /* check matrix is square */
   if (m_block != n_block) {
      Esys_setError(TYPE_ERROR, "Paso_SparseMatrix_invMain: square block size expected.");
   }
   if (Esys_noError()) {
         
      if (n_block==1) {
          #pragma omp parallel for private(i, iPtr, A11) schedule(static)
         for (i = 0; i < n; i++) {
	    iPtr= main_ptr[i];
	    A11=A_p->val[iPtr];
	    if ( ABS(A11) > 0.) { 
	        inv_diag[i]=1./A11;
	    } else {
	       failed=1;
	    }
         }
      } else if (n_block==2) {
	 #pragma omp parallel for private(i, iPtr) schedule(static)
         for (i = 0; i < n; i++) {
	       iPtr= main_ptr[i];
	       Paso_BlockOps_invM_2(&inv_diag[i*4], &A_p->val[iPtr*4], &failed);
	}
      } else if (n_block==3) {
          #pragma omp parallel for private(i, iPtr) schedule(static)
          for (i = 0; i < n; i++) {
	      iPtr= main_ptr[i];
	      Paso_BlockOps_invM_3(&inv_diag[i*9], &A_p->val[iPtr*9], &failed);
          }
      } else {    
	 #pragma omp parallel for private(i, iPtr) schedule(static)
	 for (i = 0; i < n; i++) {
	    iPtr= main_ptr[i];
	    Paso_BlockOps_Cpy_N(block_size, &inv_diag[i*block_size], &A_p->val[iPtr*block_size]);
	    Paso_BlockOps_invM_N(n_block, &inv_diag[i*block_size], &pivot[i*n_block], &failed);
	 }
      }
   }
   if (failed > 0) {
      Esys_setError(ZERO_DIVISION_ERROR, "Paso_SparseMatrix_invMain: non-regular main diagonal block.");
   }
}
void Paso_SparseMatrix_applyBlockMatrix(Paso_SparseMatrix * A_p, double* block_diag, int* pivot, double*x, double *b) {
   
   /*   inv_diag=MEMALLOC( A->numRows * A_p-> block_size,double);
   pivot=MEMALLOC( A->numRows * A->row_block_size */
   dim_t n=A_p->numRows;
   dim_t n_block=A_p->row_block_size;
   Paso_Copy(n_block*n, x,b);
   Paso_BlockOps_solveAll(n_block,n,block_diag,pivot,x);
}

