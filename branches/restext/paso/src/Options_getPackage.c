
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

/* Paso: returns the package to be used                   */

/**************************************************************/

/* Copyrights by ACcESS Australia 2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Options.h"

/**************************************************************/

index_t Paso_Options_getPackage(index_t solver,index_t package, bool_t symmetry, Paso_MPIInfo *mpi_info) {
  index_t out=PASO_PASO;
  if (package==PASO_DEFAULT) {
      if (solver==PASO_DIRECT) {
         #ifdef MKL
            out=PASO_MKL;
         #else
            #ifdef UMFPACK
              out=PASO_UMFPACK;
            #else
              #ifdef PASTIX
               out=PASO_PASTIX
              #endif
            #endif
         #endif
         if (( (out == PASO_MKL) || (out==PASO_UMFPACK) || (out == PASO_PASTIX) ) && (mpi_info->size>1) ) {  /* these packages require CSC  which is not supported with MPI */
              out= PASO_PASO;
         }
      } else {
         out=PASO_PASO;
      }
  } else if (package==PASO_PASO) {
      out=PASO_PASO;
  } else if (package==PASO_PASTIX) {
      out=PASO_PASTIX;
  } else if (package==PASO_MKL) {
      out=PASO_MKL;
  } else if (package==PASO_UMFPACK) {
      out=PASO_UMFPACK;
  } else if (package==PASO_TRILINOS) {
      out=PASO_TRILINOS;
  } else {
      Paso_setError(VALUE_ERROR,"Unidentified package.");
  }
  return out;
}
