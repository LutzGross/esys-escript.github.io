
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: SystemMatrix is saved to Harwell-Boeing format     */

/**************************************************************/

/* Copyright: ACcESS Australia 2005                           */
/* Author: imran@esscc.uq.edu.au                              */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_saveHB( Paso_SystemMatrix *A_p, char *filename_p ) {
	FILE *fileHandle_p = NULL;
        if (A_p->mpi_info->size > 1) {
	    Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: currently single processor runs are supported.\n");
    	    return;
        }
	fileHandle_p = fopen( filename_p, "w" );
	if( fileHandle_p == NULL ) {
		Paso_setError(IO_ERROR,"File could not be opened for writing.");
		return;
	}

	if ( A_p->type & MATRIX_FORMAT_CSC) {
             Paso_SparseMatrix_saveHB_CSC( A_p->mainBlock,fileHandle_p);
        } else {
      	      Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: only CSC is currently supported.\n");
	}

	/* close the file */
	fclose( fileHandle_p );

	return;
}
