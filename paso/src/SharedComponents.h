
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

/*   Paso: coupler                                            */ 

/**************************************************************/

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SHAREDCOMPONENTS
#define INC_PASO_SHAREDCOMPONENTS

#include "Paso_MPI.h"
#include "Common.h"

/**************************************************************/

typedef struct Paso_SharedComponents {

  dim_t local_length;        /* local array length shared */

  dim_t numNeighbors;        /* number of processor sharing values with this processor */

  index_t* offsetInShared; /* offsetInSharedInput[i] points to the first input value in array shared
                              for processor i. Has length numNeighbors+1 */

  Paso_MPI_rank* neighbor;  /* list of the processor sharing values with this processor */

  index_t* shared;           /* list of the (local) componets which is shared with other 
                                processors. Has length numSharedComponents */
                          
  dim_t numSharedComponents; /* = offsetInShared[numNeighbors] */

  Paso_MPIInfo *mpi_info;
  dim_t reference_counter;

} Paso_SharedComponents;



PASO_DLL_API
Paso_SharedComponents* Paso_SharedComponents_alloc(dim_t local_length,
                                                   dim_t numNeighbors,
                                                   Paso_MPI_rank* neighbor,
                                                   index_t* shared,
                                                   index_t* offsetInShared,
                                                   index_t m, index_t b,
                                                   Paso_MPIInfo *mpi_info);
                                 

PASO_DLL_API
Paso_SharedComponents* Paso_SharedComponents_getReference(Paso_SharedComponents*);

PASO_DLL_API
void Paso_SharedComponents_free(Paso_SharedComponents*);

#endif 
