
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

/* Paso: inverts the main diagonal of a SparseMatrix: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2010 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "Solver.h"

void Paso_SparseMatrix_invMain(Paso_SparseMatrix * A_p, double* inv_diag, int* pivot) {
   /*   inv_diag=MEMALLOC( A->numRows * A_p-> block_size,double);
   pivot=MEMALLOC( A->numRows * A->row_block_size */
   dim_t n=A_p->numRows;
   dim_t n_block=A_p->row_block_size;
   dim_t m_block=A_p->col_block_size;
   /* dim_t block_size=A_p->block_size; */
   dim_t i;
   register index_t iPtr;
   register double A11,A12,A13,A21,A22,A23,A31,A32,A33,D;
   index_t* main_ptr=Paso_Pattern_borrowMainDiagonalPointer(A_p->pattern);
   /* check matrix is square */
   if (m_block != n_block) {
      Paso_setError(TYPE_ERROR, "Paso_SparseMatrix_invMain: square block size expected.");
   }
   if (Paso_noError()) {
         
      if (n_block==1) {
          #pragma omp parallel for private(i, iPtr, A11) schedule(static)
         for (i = 0; i < n; i++) {
	    iPtr= main_ptr[i];
	    A11=A_p->val[iPtr];
	    if ( ABS(A11) > 0.) { 
	        inv_diag[i]=1./A11;
	    } else {
	       Paso_setError(ZERO_DIVISION_ERROR, "Paso_SparseMatrix_invMain: non-regular main diagonal block.");
	    }
         }
      } else if (n_block==2) {
	 #pragma omp parallel for private(i, iPtr, A11,A12,A21,A22,D) schedule(static)
         for (i = 0; i < n; i++) {
	       iPtr= main_ptr[i];
	       A11=A_p->val[iPtr*4];
	       A12=A_p->val[iPtr*4+2];
	       A21=A_p->val[iPtr*4+1];
	       A22=A_p->val[iPtr*4+3];
	       D = A11*A22-A12*A21;
	       if (ABS(D) > 0 ){
		 D=1./D;
		 inv_diag[i*4]= A22*D;
		 inv_diag[i*4+1]=-A21*D;
		 inv_diag[i*4+2]=-A12*D;
		 inv_diag[i*4+3]= A11*D;
	       } else {
		  Paso_setError(ZERO_DIVISION_ERROR, "Paso_SparseMatrix_invMain: non-regular main diagonal block.");
	       }
	}
      } else if (n_block==3) {
          #pragma omp parallel for private(i, iPtr,A11,A12,A13,A21,A22,A23,A31,A32,A33,D) schedule(static)
          for (i = 0; i < n; i++) {
	      iPtr= main_ptr[i];
	      A11=A_p->val[iPtr*9  ];
	      A21=A_p->val[iPtr*9+1];
	      A31=A_p->val[iPtr*9+2];
	      A12=A_p->val[iPtr*9+3];
	      A22=A_p->val[iPtr*9+4];
	      A32=A_p->val[iPtr*9+5];
	      A13=A_p->val[iPtr*9+6];
	      A23=A_p->val[iPtr*9+7];
	      A33=A_p->val[iPtr*9+8];
	      D  =  A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
	      if (ABS(D) > 0 ){
		 D=1./D;
		 inv_diag[i*9  ]=(A22*A33-A23*A32)*D;
		 inv_diag[i*9+1]=(A31*A23-A21*A33)*D;
		 inv_diag[i*9+2]=(A21*A32-A31*A22)*D;
		 inv_diag[i*9+3]=(A13*A32-A12*A33)*D;
		 inv_diag[i*9+4]=(A11*A33-A31*A13)*D;
		 inv_diag[i*9+5]=(A12*A31-A11*A32)*D;
		 inv_diag[i*9+6]=(A12*A23-A13*A22)*D;
		 inv_diag[i*9+7]=(A13*A21-A11*A23)*D;
		 inv_diag[i*9+8]=(A11*A22-A12*A21)*D;
	      } else {
		 Paso_setError(ZERO_DIVISION_ERROR, "Paso_SparseMatrix_invMain: non-regular main diagonal block.");
	      }
          }
      } else {     
	 /*
	 #pragma omp parallel for private(i, INFO) schedule(static)
	 for (i = 0; i < n; i++) {
	    iPtr= main_ptr[i];
	    memcpy((void*)&(diag[i*block_size], (void*)(&A_p->val[iPtr*block_size]), sizeof(double)*(size_t)block_size);
	    
	    int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	    lda, integer *ipiv, integer *info);
	    
	    dgetrf_(n_block, m_block, &(diag[i*block_size], n_block, pivot[i*n_block], info );
	    if (INFO >0 ) {
	       Paso_setError(ZERO_DIVISION_ERROR, "Paso_SparseMatrix_invMain: non-regular main diagonal block.");
          }
	 */
	 Paso_setError(TYPE_ERROR, "Paso_SparseMatrix_invMain: Right now there is support block size less than 4 only");
      }
   }
}
void Paso_SparseMatrix_applyBlockMatrix(Paso_SparseMatrix * A_p, double* block_diag, int* pivot, double*x, double *b) {
   
   /*   inv_diag=MEMALLOC( A->numRows * A_p-> block_size,double);
   pivot=MEMALLOC( A->numRows * A->row_block_size */
   dim_t n=A_p->numRows;
   dim_t n_block=A_p->row_block_size;
  
   Paso_Solver_applyBlockDiagonalMatrix(n_block,n,block_diag,pivot,x,b);
}

