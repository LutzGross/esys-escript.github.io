/* $Id:$ */

/*
********************************************************************************
*               Copyright  2007 by ACcESS MNRF                                 *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/*   Paso: coupler                                            */ 

/**************************************************************/

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SHAREDCOMPONENTS
#define INC_PASO_SHAREDCOMPONENTS

#include "Distribution.h"

/**************************************************************/

typedef struct Paso_SharedComponents {
  dim_t numComponents;       /* number of components used by this process */
  index_t* globalComponent;  /* assigns a global id to each component     
                                globalComponent[i] is on this processor for i=0,..,myNumComponents,
                                             and
                                globalComponent[i] is on a remote processor for i=myNumComponents,..,numComponents
                                             more accurate 
                                globalComponent[i] is on processor neighbours[p]
                                          for i=offsetInShared[p]+myNumComponents ... offsetInShared[p+1]+myNumComponents */
  index_t* ordering;         /* provides the reordering of the input globalComponent 
                                   Paso_SharedComponents.globalComponent[ordering[i]] = globalComponent[i]  */

  Paso_Distribution *distribution; /* global distribution of the components */

  dim_t myNumComponents;     /* number of components used by this process */
  dim_t numNeighbors;        /* number of processor sharing values with this processor */
  Paso_MPI_rank* neighbor;  /* list of the processor sharing values with this processor */

  dim_t numSharedComponents; /* = numComponents-myNumComponents */
  index_t* shared;           /* list of the (local) componets which is shared with other 
                                processors. Has length numSharedComponents */
                          
  index_t* offsetInShared; /* offsetInSharedInput[i] points to the first input value in array shared
                              for processor i. Has length numSharedComponents */

  Paso_MPIInfo *mpi_info;
  dim_t reference_counter;
} Paso_SharedComponents;


Paso_SharedComponents* Paso_SharedComponents_alloc(dim_t numComponents,
                                                   index_t* globalComponent, 
                                                   index_t m, index_t b,
                                                   Paso_Distribution *distribution);
                                 
Paso_SharedComponents* Paso_SharedComponents_getReference(Paso_SharedComponents*);
void Paso_SharedComponents_free(Paso_SharedComponents*);

#endif 
