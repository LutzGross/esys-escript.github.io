
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

/* Paso: returns the matrix format requested by a particular linear solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"
#include "SystemMatrix.h"

/**************************************************************/

index_t Paso_SystemMatrix_getSystemMatrixTypeId(const index_t solver,const index_t preconditioner, const index_t package,const  bool_t symmetry, Paso_MPIInfo *mpi_info) {
  index_t true_package;
  index_t out=MATRIX_FORMAT_DEFAULT;
  true_package=Paso_Options_getPackage(solver,package,symmetry, mpi_info);

  switch(true_package)  {

     case PASO_PASO:
      /* AMG on blocksize > 1 not supported yet */
  /*     if (preconditioner == PASO_AMG) {
          out=MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
       } else {
  */
          out=MATRIX_FORMAT_DEFAULT;
   /*    } */
       break;


     case PASO_MKL:
       out=MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1;
       /* if (solver == PASO_CHOLEVSKY) out+=MATRIX_FORMAT_SYM */
       break;

     case PASO_UMFPACK:
       if (mpi_info->size > 1) {
           Paso_setError(VALUE_ERROR,"The selected solver UMFPACK requires CSC format which is not supported with MPI.");
           return -1;
       }
       out=MATRIX_FORMAT_CSC + MATRIX_FORMAT_BLK1;
      break;

    case PASO_TRILINOS:
      out=MATRIX_FORMAT_TRILINOS_CRS + MATRIX_FORMAT_BLK1; /* Distributed CRS */
      break;

     default:
        Paso_setError(VALUE_ERROR,"unknown package code");
  }
  return out;
}
