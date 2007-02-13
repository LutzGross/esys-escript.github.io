/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/
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
  index_t out=MATRIX_FORMAT_DEFAULT;
  package=Paso_Options_getPackage(solver,package,symmetry);

  switch(package)  {

     case PASO_PASO:
       out=MATRIX_FORMAT_DEFAULT;
       break;

     case PASO_SCSL:
       out=MATRIX_FORMAT_CSC + MATRIX_FORMAT_BLK1;
       /* if (solver == PASO_CHOLEVSKY) out+=MATRIX_FORMAT_SYM */
       break;

     case PASO_MKL:
       out=MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1;
       /* if (solver == PASO_CHOLEVSKY) out+=MATRIX_FORMAT_SYM */
       break;

     case PASO_UMFPACK:
       out=MATRIX_FORMAT_CSC + MATRIX_FORMAT_BLK1;
      break;

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
