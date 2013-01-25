
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

/*   Paso: system matrix pattern                            */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIXPATTERN
#define INC_PASO_SYSTEMMATRIXPATTERN

#include "Distribution.h"
#include "Pattern.h"
#include "Coupler.h"

/**************************************************************/

typedef struct Paso_SystemMatrixPattern {
  int type;

  Paso_MPIInfo *mpi_info;

  Paso_Pattern* mainPattern;
  Paso_Pattern* couplePattern;
  Paso_Coupler* coupler;
  Paso_Distribution *output_distribution; 
  Paso_Distribution *input_distribution; 

  dim_t reference_counter;
  

} Paso_SystemMatrixPattern;


/*  interfaces: */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type, Paso_Distribution* output_distribution, Paso_Distribution* input_distribution, Paso_Pattern* mainPattern, Paso_Pattern* couplePattern, Paso_Coupler* coupler);
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern*);
void Paso_SystemMatrixPattern_free(Paso_SystemMatrixPattern*);
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern,
                                           int type, dim_t output_block_size,dim_t input_block_size);

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */
