
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix saving to Harwell-Boeing format         */

/**************************************************************/

/* Copyright: ACcESS Australia 2005                           */
/* Author: imran@esscc.uq.edu.au                              */

/**************************************************************/

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
             Paso_SparseMatrix_saveHB_CSC( A_p->mainBlock,fileHandle_p);
        } else {
      	      Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: Only CSC is currently supported.\n");
	}

	/* close the file */
	fclose( fileHandle_p );
}

