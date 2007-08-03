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

/* Paso: SystemMatrix :                           */
/*  sets the values of the system matrix to a value */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005, 2006, 2007 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/**************************************************************/

void  Paso_SystemMatrix_setValues(Paso_SystemMatrix* in,double value) {

  if (in!=NULL) {
      Paso_SparseMatrix_setValues(in->mainBlock, value);
      Paso_SparseMatrix_setValues(in->coupleBlock, value);
      in->normalizer_is_valid=FALSE;
  }
}
