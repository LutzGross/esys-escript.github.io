/* $Id$ */

/**********************************************************************/

/* Finley: Solver: set the pointer mainDiag of the main diagonal in A_p->value */
/*                 if mainDiag[i]<0 row has no main diagonal entry. */
/*                 matrix has to be given in CSR format             */

/**********************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "System.h"
#include "Solver.h"

/**************************************************************/

void Finley_Solver_getMainDiagonal(Finley_SystemMatrix* A_p,maybelong* mainDiag) {
#if ITERATIVE_SOLVER == NO_LIB
  int irow,iptr;
  maybelong n=A_p->pattern->n_ptr;
  #pragma omp parallel for private(irow,iptr) schedule(static)
  for (irow=0;irow<n;irow++) {
       mainDiag[irow]=-1;
       for (iptr=A_p->pattern->ptr[irow];iptr<A_p->pattern->ptr[irow+1]; iptr++) {
            if (A_p->pattern->index[iptr]==irow) {
               mainDiag[irow]=iptr;
               break;
            }
       }
  }
#endif
}

