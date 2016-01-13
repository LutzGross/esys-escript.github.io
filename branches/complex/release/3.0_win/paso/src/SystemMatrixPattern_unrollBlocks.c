
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

/* Paso: SystemMatrixPattern_unrollBlocks */

/**************************************************************/

/* Author: gross@access.edu.au */

/**************************************************************/

#include "SystemMatrixPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */


Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern, 
                                           int type, dim_t output_block_size,dim_t input_block_size) {
  Paso_SystemMatrixPattern*out=NULL;
  Paso_Pattern *new_mainPattern=NULL, *new_col_couplePattern=NULL, *new_row_couplePattern=NULL;
  Paso_Distribution* new_output_distribution=NULL, *new_input_distribution=NULL;
  Paso_Connector *new_col_connector=NULL, *new_row_connector=NULL;

  if ( ( output_block_size == 1 ) && (input_block_size == 1) && ((pattern->type & PATTERN_FORMAT_OFFSET1) == (type & PATTERN_FORMAT_OFFSET1) ) ) {
     out = Paso_SystemMatrixPattern_getReference(pattern);
  } else {
     new_mainPattern=Paso_Pattern_unrollBlocks(pattern->mainPattern,type,output_block_size,input_block_size);
     new_col_couplePattern=Paso_Pattern_unrollBlocks(pattern->col_couplePattern,type,output_block_size,input_block_size);
     new_row_couplePattern=Paso_Pattern_unrollBlocks(pattern->row_couplePattern,type,output_block_size,input_block_size);
     if (output_block_size>1) {
          new_output_distribution=Paso_Distribution_alloc(pattern->output_distribution->mpi_info,
                                                          pattern->output_distribution->first_component,
                                                          output_block_size,0);
          new_row_connector=Paso_Connector_unroll(pattern->row_connector,output_block_size);
     } else {
          new_output_distribution=Paso_Distribution_getReference(pattern->output_distribution);
          new_row_connector= Paso_Connector_getReference(pattern->row_connector);
     }
     if (input_block_size>1) {
          new_input_distribution=Paso_Distribution_alloc(pattern->input_distribution->mpi_info,
                                                          pattern->input_distribution->first_component,
                                                          input_block_size,0);
          new_col_connector=Paso_Connector_unroll(pattern->col_connector,input_block_size);
     } else {
          new_input_distribution=Paso_Distribution_getReference(pattern->input_distribution);
          new_col_connector=Paso_Connector_getReference(pattern->col_connector);
     }
     if (Paso_noError()) {
        out=Paso_SystemMatrixPattern_alloc(type,
                                           new_output_distribution,
                                           new_input_distribution,
                                           new_mainPattern,
                                           new_col_couplePattern,
                                           new_row_couplePattern,
                                           new_col_connector,
                                           new_row_connector);
     }
     Paso_Pattern_free(new_mainPattern);
     Paso_Pattern_free(new_col_couplePattern);
     Paso_Pattern_free(new_row_couplePattern);
     Paso_Distribution_free(new_output_distribution);
     Paso_Distribution_free(new_input_distribution);
     Paso_Connector_free(new_row_connector);
     Paso_Connector_free(new_col_connector);
  }
  if (Paso_noError()) {
      return out;
  } else {
     Paso_SystemMatrixPattern_free(out);
     return NULL;
  }
}
