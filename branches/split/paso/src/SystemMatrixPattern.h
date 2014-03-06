
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

/*   Paso: system matrix pattern                              */

/************************************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#ifndef INC_PASO_SYSTEMMATRIXPATTERN
#define INC_PASO_SYSTEMMATRIXPATTERN

#include "Distribution.h"
#include "Pattern.h"
#include "Coupler.h"

/************************************************************************************/

typedef struct Paso_SystemMatrixPattern {
  int type;

  Esys_MPIInfo *mpi_info;

  
  Paso_Pattern* mainPattern;
  Paso_Pattern* col_couplePattern;
  Paso_Pattern* row_couplePattern;
  Paso_Connector* col_connector;
  Paso_Connector* row_connector;
  Paso_Distribution *output_distribution; 
  Paso_Distribution *input_distribution; 

  dim_t reference_counter;
  

} Paso_SystemMatrixPattern;


/*  interfaces: */


PASO_DLL_API
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type, Paso_Distribution* output_distribution, Paso_Distribution* input_distribution, Paso_Pattern* mainPattern, Paso_Pattern* col_couplePattern, Paso_Pattern* row_couplePattern, Paso_Connector* col_connector, Paso_Connector* row_connector);

PASO_DLL_API
Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_getReference(Paso_SystemMatrixPattern*);

PASO_DLL_API
void Paso_SystemMatrixPattern_free(Paso_SystemMatrixPattern*);

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern* pattern,
                                           int type, dim_t output_block_size,dim_t input_block_size);
index_t Paso_SystemMatrixPattern_getNumOutput(Paso_SystemMatrixPattern*);

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */
