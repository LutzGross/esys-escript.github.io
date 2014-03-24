
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

/* Paso: SystemMatrix saving to Harwell-Boeing format         */

/************************************************************************************/

/* Copyright: ACcESS Australia 2005                           */
/* Author: imran@esscc.uq.edu.au                              */

/************************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_saveHB( Paso_SystemMatrix *A_p, char *filename_p ) {
	FILE *fileHandle_p = NULL;
        if (A_p->mpi_info->size > 1) {
	    Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: Only single processor runs are supported.\n");
    	    return;
        }
	fileHandle_p = fopen( filename_p, "w" );
	if( fileHandle_p == NULL ) {
		Esys_setError(IO_ERROR,"Paso_SystemMatrix_saveHB: File could not be opened for writing.");
		return;
	}

	if ( A_p->type & MATRIX_FORMAT_CSC) {
        paso::SparseMatrix_saveHB_CSC( A_p->mainBlock,fileHandle_p);
        } else {
      	      Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: Only CSC is currently supported.\n");
	}

	/* close the file */
	fclose( fileHandle_p );
}

