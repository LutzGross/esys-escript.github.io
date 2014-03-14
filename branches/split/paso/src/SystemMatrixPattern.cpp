
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

/* Paso: SystemMatrixPattern */

/************************************************************************************/
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/************************************************************************************/

/* allocates a SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type,
                                                         Paso_Distribution *output_distribution,
                                                         Paso_Distribution *input_distribution,
                                                         Paso_Pattern* mainPattern,
                                                         Paso_Pattern* col_couplePattern,
                                                         Paso_Pattern* row_couplePattern,
                                                         Paso_Connector* col_connector,
                                                         Paso_Connector* row_connector) 
{
  Paso_SystemMatrixPattern*out=NULL;
  Esys_resetError();

  if (output_distribution->mpi_info != input_distribution->mpi_info ) {
     Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and input_distribution MPI communicators don't match.");
     return NULL;
  }
  if (output_distribution->mpi_info != col_connector->mpi_info ) {
     Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and col_connector MPI communicators don't match.");
     return NULL;
  }
  if (output_distribution->mpi_info != row_connector->mpi_info ) {
     Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and row_connector MPI communicators don't match.");
     return NULL;
  }


  if (mainPattern->type != type)  {
      Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of mainPattern does not match expected type.");
  }
  if (col_couplePattern->type != type)  {
      Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of col_couplePattern does not match expected type.");
  }
  if (row_couplePattern->type != type)  {
      Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of row_couplePattern does not match expected type.");
  }
  if (col_couplePattern->numOutput != mainPattern->numOutput) {
            Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of outputs for couple and main pattern don't match.");
  }
  if (mainPattern->numOutput !=  Paso_Distribution_getMyNumComponents(output_distribution)) {
      Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of outputs and given distribution don't match.");
  }
  if (mainPattern->numInput != Paso_Distribution_getMyNumComponents(input_distribution)) {
     Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of input for main pattern and number of send components in connector don't match.");
  }
  if (col_couplePattern->numInput != col_connector->recv->numSharedComponents) {
     Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of inputs for column couple pattern and number of received components in connector don't match.");
  }
  if (row_couplePattern->numOutput != row_connector->recv->numSharedComponents) {
     Esys_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of inputs for row couple pattern and number of received components in connector don't match.");
  }

  out=new Paso_SystemMatrixPattern;
  if (Esys_checkPtr(out)) return NULL;
  out->type=type;
  out->reference_counter=1;
  out->mainPattern=Paso_Pattern_getReference(mainPattern);
  out->row_couplePattern=Paso_Pattern_getReference(row_couplePattern);
  out->col_couplePattern=Paso_Pattern_getReference(col_couplePattern);
  out->row_connector=Paso_Connector_getReference(row_connector);
  out->col_connector=Paso_Connector_getReference(col_connector);
  out->output_distribution=Paso_Distribution_getReference(output_distribution);
  out->input_distribution=Paso_Distribution_getReference(input_distribution);
  out->mpi_info= output_distribution->mpi_info;
  #ifdef Paso_TRACE
  printf("Paso_SystemMatrixPattern_alloc: system matrix pattern has been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_getReference(Paso_SystemMatrixPattern* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a SystemMatrixPattern: */

void Paso_SystemMatrixPattern_free(Paso_SystemMatrixPattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_Pattern_free(in->mainPattern);
        Paso_Pattern_free(in->row_couplePattern);
        Paso_Pattern_free(in->col_couplePattern);
        Paso_Connector_free(in->row_connector);
        Paso_Connector_free(in->col_connector);
        Paso_Distribution_free(in->output_distribution);
        Paso_Distribution_free(in->input_distribution);
        delete in;
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrixPattern_free: system matrix pattern has been deallocated.\n");
        #endif
     }
   }
}
dim_t Paso_SystemMatrixPattern_getNumOutput(Paso_SystemMatrixPattern* in) {
    if (in!=NULL) {
       return 0;
    } else {
       return in->mainPattern->numOutput;
    }
}

