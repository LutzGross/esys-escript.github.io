/* $Id$ */

/**************************************************************/

/* Finley: returns the matrix format requested by a particular linear solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "System.h"

/**************************************************************/

int Finley_SystemMatrix_getSystemMatrixTypeId(int solver,int symmetry) {
  int out=CSR;
  if (solver==ESCRIPT_DIRECT || solver==ESCRIPT_CHOLEVSKY) {
      #if DIRECT_SOLVER == SGI_SCSL
         out= symmetry ? CSC_BLK1_SYM : CSC_BLK1;
      #else 
         out= CSR;
      #endif
  }  else {
     #if ITERATIVE_SOLVER == SGI_SCSL
        out= CSR_BLK1;
     #else
        out= CSR;
     #endif
  }
  return out;
}

/*
 * $Log$
 * Revision 1.3  2005/06/09 05:38:02  jgs
 * Merge of development branch back to main trunk on 2005-06-09
 *
 * Revision 1.2  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.2  2005/05/13 05:48:17  gross
 * some changes to get this running under gcc
 *
 * Revision 1.1.2.1  2004/11/15 00:59:05  gross
 * and anotther missing file
 *
 */
