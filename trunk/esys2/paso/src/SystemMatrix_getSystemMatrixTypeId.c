/* $Id$ */

/**************************************************************/

/* Paso: returns the matrix format requested by a particular linear solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"
#include "SystemMatrix.h"

/**************************************************************/

index_t Paso_SystemMatrix_getSystemMatrixTypeId(index_t solver,index_t package, bool_t symmetry) {
  index_t out=CSR;
  package=Paso_Options_getPackage(solver,package,symmetry);

  switch(package)  {

     case PASO_PASO:
       out=CSR;
       break;

     case PASO_SCSL:
       out= symmetry ? CSC_BLK1_SYM : CSC_BLK1;
       break;
/*
     case PASO_MKL:
       out= CSR_BLK1;
       break;
*/

/*
     case PASO_UMFPACK:
       out= CSR_BLK1;
      break;
*/

     default:
        Paso_setError(VALUE_ERROR,"unknown package code");
  }
  return out;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
