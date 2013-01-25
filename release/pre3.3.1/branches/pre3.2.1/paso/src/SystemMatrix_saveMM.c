
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "SystemMatrix.h"
#include "mmio.h"

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix * A_p, char * fileName_p) {
   if (A_p->mpi_info->size > 1) {
      Esys_setError(IO_ERROR, "Paso_SystemMatrix_saveMM: support single processor only");
   } else {
       Paso_SparseMatrix_saveMM(A_p->mainBlock,fileName_p);
   }
   return;
}
