
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*   Paso: coupler                                            */ 

/************************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#ifndef INC_PASO_SHAREDCOMPONENTS
#define INC_PASO_SHAREDCOMPONENTS

#include "Common.h"
#include "esysUtils/Esys_MPI.h"

/************************************************************************************/

typedef struct Paso_SharedComponents {

  dim_t local_length;        /* local array length shared */

  dim_t numNeighbors;        /* number of processor sharing values with this processor */

  index_t* offsetInShared; /* offsetInSharedInput[i] points to the first input value in array shared
                              for processor i. Has length numNeighbors+1 */

  Esys_MPI_rank* neighbor;  /* list of the processor sharing values with this processor */

  index_t* shared;           /* list of the (local) components which are shared with other 
                                processors. Has length numSharedComponents */
                          
  dim_t numSharedComponents; /* = offsetInShared[numNeighbors] */

  Esys_MPIInfo *mpi_info;
  dim_t reference_counter;

} Paso_SharedComponents;



PASO_DLL_API
Paso_SharedComponents* Paso_SharedComponents_alloc(dim_t local_length,
                                                   dim_t numNeighbors,
                                                   Esys_MPI_rank* neighbor,
                                                   index_t* shared,
                                                   index_t* offsetInShared,
                                                   index_t m, index_t b,
                                                   Esys_MPIInfo *mpi_info);
                                 

PASO_DLL_API
Paso_SharedComponents* Paso_SharedComponents_getReference(Paso_SharedComponents*);

PASO_DLL_API
void Paso_SharedComponents_free(Paso_SharedComponents*);

#endif 
