/* $Id$ */

/**************************************************************/

/* Finley: jacobi preconditioner: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"
#include "Solver.h"

/**************************************************************/

/* inverse preconditioner setup */

void Finley_Solver_setJacobi(Finley_SystemMatrix * A_p) {
#if ITERATIVE_SOLVER == NO_LIB
  Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
  int n = A_p->num_cols * A_p->col_block_size;
  int iCol, iRow, iPtr;
  double rowSum;
  int rowSign;
  /* check matrix is square */
  if (A_p->col_block_size !=1) {
    Finley_ErrorCode = TYPE_ERROR;
    sprintf(Finley_ErrorMsg, "Jacobi preconditioner requires block size 1.");
    return;
  }
  #ifdef FINLEY_SOLVER_TRACE
  printf("Jacobi preconditioner is used.\n");
  #endif
  /* allocate vector to hold main diagonal entries: */
  prec->values = (double *) MEMALLOC(sizeof(double) * n);
  if (Finley_checkPtr(prec->values)) return;

  /* TODO: validate for CSC */
  /* TODO: block_size>1 */
  switch(A_p->type) {
     case CSR:
   
       #pragma omp parallel for private(iRow, iPtr,rowSum,rowSign) schedule(static)
       for (iRow = 0; iRow < A_p->num_rows; iRow++) {
         rowSum = 0.0;
         rowSign = 1;
         for (iPtr = A_p->ptr[iRow]; iPtr < A_p->ptr[iRow + 1]; iPtr++) {
   	        rowSum += ABS(A_p->val[iPtr]);
   	        if (iRow == A_p->index[iPtr]) {
   	          rowSign = rowSign - 2 * (A_p->val[iPtr] < 0.0);
   	        }
         }	/* for iPtr */
         if (rowSum>0) {
   	        prec->values[iRow] = ((double) rowSign) / rowSum;
         } else {
   	        prec->values[iRow] = (double) rowSign;
         }
       } /* for iRow */
       break;

     case CSC:
       for (iRow = 0; iRow < A_p->num_rows; iRow++) {
           rowSum = 0.0;
           rowSign = 1;
           for (iCol = 0; iCol < A_p->num_cols; iCol++) {
       	     for (iPtr = A_p->ptr[iCol] ; (A_p->index[iPtr] < iRow ) && (iPtr < A_p->ptr[iCol + 1]); iPtr++);
   	        if (iRow  == A_p->index[iPtr]) {
   	          rowSum += ABS(A_p->val[iPtr]);
   	          if (iCol == A_p->index[iPtr]) {
   	            rowSign = rowSign - 2 * (A_p->val[iPtr] < 0.0);
   	          }
           	}
           } /* for iCol */
           if (rowSum>0)
   	      prec->values[iRow] = ((double) rowSign) / rowSum;
           else
   	      prec->values[iRow] = ((double) rowSign);
       } /* for iRow */
       break;
     default:
       Finley_ErrorCode = TYPE_ERROR;
       sprintf(Finley_ErrorMsg, "Unknown matrix type.");
       return;
  } /* switch A_p->type */
#endif
} 

/**************************************************************/

/* inverse preconditioner solution using raw pointers */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Finley_Solver_solveJacobi(Finley_SystemMatrix * A_p, double * x, double * b) {
#if ITERATIVE_SOLVER == NO_LIB
  Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
  int n = A_p->num_cols * A_p->col_block_size, i;
  #pragma omp for private(i) schedule(static)
  for (i = 0; i < n; i++) {
     x[i] = prec->values[i] * b[i]; 
  }
  return;
#endif
}

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:58  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
