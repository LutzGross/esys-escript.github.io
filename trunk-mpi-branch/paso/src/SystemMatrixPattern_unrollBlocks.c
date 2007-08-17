/* $Id: SystemMatrixPattern_unrollBlocks.c 1140 2007-05-15 03:23:17Z ksteube $ */

/*
********************************************************************************
*               Copyright 2006, 2007 by ACcESS MNRF                            *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

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
  Paso_Pattern *new_mainPattern=NULL, *new_couplePattern=NULL;
  Paso_Distribution* new_output_distribution=NULL, *new_input_distribution=NULL;
  Paso_Coupler *new_coupler=NULL;

  new_mainPattern=Paso_Pattern_unrollBlocks(pattern->mainPattern,type,output_block_size,input_block_size);
  new_couplePattern=Paso_Pattern_unrollBlocks(pattern->couplePattern,type,output_block_size,input_block_size);
  new_output_distribution=Paso_Distribution_alloc(pattern->output_distribution->mpi_info,
                                                  pattern->output_distribution->first_component,
                                                  output_block_size,0);
  new_input_distribution=Paso_Distribution_alloc(pattern->input_distribution->mpi_info,
                                                  pattern->input_distribution->first_component,
                                                  input_block_size,0);
  new_coupler=Paso_Coupler_unroll(pattern->coupler,input_block_size);
  if (Paso_noError()) {
     out=Paso_SystemMatrixPattern_alloc(type,
                                        new_output_distribution,
                                        new_input_distribution,
                                        new_mainPattern,
                                        new_couplePattern,
                                        new_coupler);
  }
  Paso_Pattern_free(new_mainPattern);
  Paso_Pattern_free(new_couplePattern);
  Paso_Distribution_free(new_output_distribution);
  Paso_Distribution_free(new_input_distribution);
  Paso_Coupler_free(new_coupler);

  if (Paso_noError()) {
     return out;
  } else {
     Paso_SystemMatrixPattern_free(out);
     return NULL;
  }
}
