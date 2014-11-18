
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

/* Paso: SystemMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "SystemMatrix.h"

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix * A_p, char * fileName_p) {
  FILE * fileHandle_p = NULL;

  if (A_p->mpi_info->size > 1) {
       Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveHB: currently single processor runs are supported.\n");
       return;
  }
  if (A_p->type & MATRIX_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM does not support symmetric storage scheme");
    return;
  }
  /* open the file */
  fileHandle_p = fopen(fileName_p, "w");
  if (fileHandle_p==NULL) {
    Paso_setError(IO_ERROR,"file could not be opened for writing");
    return;
  }

  if (A_p->type & MATRIX_FORMAT_CSC) {
    Paso_SparseMatrix_saveHB_CSC( A_p->mainBlock,fileHandle_p);
  } else { 
    /* Paso_SparseMatrix_saveHB_CSR( A_p->mainBlock,fileHandle_p); */
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM does not support CSR yet.");
    return;
  }

  /* close the file */
  fclose(fileHandle_p);
  
  return;
}