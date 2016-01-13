
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/
 
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

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
  Paso_resetError();

  if (output_distribution->mpi_info != input_distribution->mpi_info ) {
     Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and input_distribution mpi communicator don't match.");
     return NULL;
  }
  if (output_distribution->mpi_info != col_connector->mpi_info ) {
     Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and col_connector mpi communicator don't match.");
     return NULL;
  }
  if (output_distribution->mpi_info != row_connector->mpi_info ) {
     Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: output_distribution and row_connector mpi communicator don't match.");
     return NULL;
  }


  if (mainPattern->type != type)  {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of mainPattern does not match expected type.");
  }
  if (col_couplePattern->type != type)  {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of col_couplePattern does not match expected type.");
  }
  if (row_couplePattern->type != type)  {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: type of row_couplePattern does not match expected type.");
  }
  if (col_couplePattern->numOutput != mainPattern->numOutput) {
            Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of output for couple and main pattern don't match.");
  }
  if (mainPattern->numOutput !=  Paso_Distribution_getMyNumComponents(output_distribution)) {
      Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of output and given distribution don't match.");
  }
  if (mainPattern->numInput != Paso_Distribution_getMyNumComponents(input_distribution)) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of input for main pattern and number of send components in connector don't match.");
  }
  if (col_couplePattern->numInput != col_connector->recv->numSharedComponents) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of inputs for column couple pattern and number of received components in connector don't match.");
  }
  if (row_couplePattern->numOutput != row_connector->recv->numSharedComponents) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: number of inputs for row couple pattern and number of received components in connector don't match.");
  }
  if (mainPattern->output_block_size != col_couplePattern->output_block_size) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: output block sizes of main and column couple pattern do not match.");
  }
  if (mainPattern->input_block_size != col_couplePattern->input_block_size) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: input block sizes of main and column couple pattern do not match.");
  }
  if (mainPattern->output_block_size != row_couplePattern->output_block_size) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: output block sizes of main and row couple pattern do not match.");
  }
  if (mainPattern->input_block_size != col_couplePattern->input_block_size) {
     Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: input block sizes of main and row couple pattern do not match.");
  }

  out=MEMALLOC(1,Paso_SystemMatrixPattern);
  if (Paso_checkPtr(out)) return NULL;
  out->type=type;
  out->reference_counter=1;
  out->mainPattern=Paso_Pattern_getReference(mainPattern);
  out->row_couplePattern=Paso_Pattern_getReference(row_couplePattern);
  out->col_couplePattern=Paso_Pattern_getReference(col_couplePattern);
  out->row_connector=Paso_Connector_getReference(row_connector);
  out->col_connector=Paso_Connector_getReference(col_connector);
  out->output_distribution=Paso_Distribution_getReference(output_distribution);
  out->input_distribution=Paso_Distribution_getReference(input_distribution);
  out->mpi_info= Paso_MPIInfo_getReference(output_distribution->mpi_info);
  #ifdef Paso_TRACE
  printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern* in) {
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
        Paso_MPIInfo_free(in->mpi_info);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrixPattern_free: system matrix pattern as been deallocated.\n");
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
