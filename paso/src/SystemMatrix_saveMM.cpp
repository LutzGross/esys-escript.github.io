
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/* Paso: SystemMatrix is saved to Matrix Market format      */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/************************************************************************************/

#include "SystemMatrix.h"
#include "mmio.h"

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix * A_p, char * fileName_p) {
   if (A_p->mpi_info->size > 1) {
      Esys_setError(IO_ERROR, "Paso_SystemMatrix_saveMM: Only single processor supported.");
   } else {
       paso::SparseMatrix_saveMM(A_p->mainBlock,fileName_p);
   }
}

