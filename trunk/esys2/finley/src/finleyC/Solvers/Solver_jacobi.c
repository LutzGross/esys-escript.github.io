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
#include "Util.h"
#include "Common.h"

/**************************************************************/

/* inverse preconditioner setup */

void Finley_Solver_setJacobi(Finley_SystemMatrix * A_p) {
#if ITERATIVE_SOLVER == NO_LIB
  Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
  int n = A_p->num_cols;
  int n_block=A_p->row_block_size;
  int block_size=A_p->block_size;
  int i,iPtr,info;
  /* check matrix is square */
  if (A_p->col_block_size !=A_p->row_block_size) {
    Finley_ErrorCode = TYPE_ERROR;
    sprintf(Finley_ErrorMsg, "Jacobi preconditioner square block size.");
    return;
  }
  /* check matrix is square */
  if (n_block>3) {
    Finley_ErrorCode = TYPE_ERROR;
    sprintf(Finley_ErrorMsg, "Right now the Jacobi preconditioner supports block size less than 4 only");
    return;
  }
  /* allocate vector to hold main diagonal entries: */
  prec->values = MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
  prec->pivot = MEMALLOC(  ((size_t) n) * ((size_t) n_block),int);
  if (! (Finley_checkPtr(prec->values) || Finley_checkPtr(prec->pivot)) ) {

     if (n_block==1) {
        #pragma omp parallel for private(i, iPtr) schedule(static)
        for (i = 0; i < A_p->pattern->n_ptr; i++) {
           /* find main diagonal */
           for (iPtr = A_p->pattern->ptr[i]; iPtr < A_p->pattern->ptr[i + 1]; iPtr++) {
               if (A_p->pattern->index[iPtr]==i+INDEX_OFFSET) {
                   if (ABS(A_p->val[iPtr])>0.) {
                      prec->values[i]=1./A_p->val[iPtr];
                   } else {
                      prec->values[i]=1.;
                   }
                   break;
               }
           }
        }
     } else {
        #pragma omp parallel for private(i, iPtr,info) schedule(static)
        for (i = 0; i < A_p->pattern->n_ptr; i++) {
           /* find main diagonal */
           for (iPtr = A_p->pattern->ptr[i]; iPtr < A_p->pattern->ptr[i + 1]; iPtr++) {
               if (A_p->pattern->index[iPtr]==i+INDEX_OFFSET) {
                   info=Finley_Util_SmallMatLU(n_block,&A_p->val[iPtr*block_size],&prec->values[i*block_size],&prec->pivot[i*n_block]);
                   if (info>0) {
                        Finley_ErrorCode = ZERO_DIVISION_ERROR;
                        sprintf(Finley_ErrorMsg, "non-regular main diagonal block in row %d",i);
                   }
                   break;
               }
           }
        }
     }
  }
#endif
} 

/**************************************************************/

/* inverse preconditioner solution using raw pointers */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Finley_Solver_solveJacobi(Finley_SystemMatrix * A_p, double * x, double * b) {
#if ITERATIVE_SOLVER == NO_LIB
  Finley_Solver_Preconditioner* prec=(Finley_Solver_Preconditioner*) (A_p->iterative);
  int n_block=A_p->row_block_size;
  int block_size=A_p->block_size;
  int i;

  if (n_block==1) {
     #pragma omp for private(i) schedule(static)
     for (i = 0; i < A_p->pattern->n_ptr; i++) {
        x[i] = prec->values[i] * b[i]; 
     }
  } else {
     #pragma omp for private(i) schedule(static)
     for (i = 0; i < A_p->pattern->n_ptr; i++) {
        Finley_Util_SmallMatForwardBackwardSolve(n_block,1,&prec->values[i*block_size],&prec->pivot[i*n_block],&x[i*n_block],&b[i*n_block]);
     }
  }
  return;
#endif
}

/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:32  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.3  2004/12/07 10:12:06  gross
 * GMRES added
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:17  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:21  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:58  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
