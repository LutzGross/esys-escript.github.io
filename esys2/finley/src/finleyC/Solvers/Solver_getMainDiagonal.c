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
  maybelong n=A_p->num_rows;
  #pragma omp parallel for private(irow,iptr) schedule(static)
  for (irow=0;irow<n;irow++) {
       mainDiag[irow]=-1;
       for (iptr=A_p->ptr[irow];iptr<A_p->ptr[irow+1]; iptr++) {
            if (A_p->index[iptr]==irow) {
               mainDiag[irow]=iptr;
               break;
            }
       }
  }
#endif
}

/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:47  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:58  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:14  gross
 * Finley C code has been included
 *
 *
 */
