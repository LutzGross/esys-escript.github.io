
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix :                           */
/*  sets the values of the system matrix to a value */

/**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

/**************************************************************/

void  Paso_SystemMatrix_setValues(Paso_SystemMatrix* in,double value) {

  if (in!=NULL) {
      Paso_SparseMatrix_setValues(in->mainBlock, value);
      Paso_SparseMatrix_setValues(in->col_coupleBlock, value);
      Paso_SparseMatrix_setValues(in->row_coupleBlock, value);
      in->normalizer_is_valid=FALSE;
  }
}
