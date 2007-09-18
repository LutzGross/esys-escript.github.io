
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

/*                                                                                                                     */
/* NodeMapping provides a mapping from the local nodes typically to the degrees of freedom,                            */
/*    the reduced degrees of freedom or the reduced node set                                                           */
/*                                                                                                                     */

#ifndef INC_FINLEY_NODEMAPPING
#define INC_FINLEY_NODEMAPPING

#include "paso/Paso_MPI.h"


struct Finley_NodeMapping {
  dim_t numNodes; /* number of FEM nodes */
  index_t *target; /* target[i] defines the target if FEM  node i =0,...,numNodes */
  index_t unused;  /* target[i]=unused defines that no target is defined for FEM  node i */
  dim_t numTargets; /* number of targets */
  index_t *map;  /* maps the target nodes back to the FEM nodes: target[map[i]]=i */
  dim_t reference_counter;
};
typedef struct Finley_NodeMapping Finley_NodeMapping;

Finley_NodeMapping* Finley_NodeMapping_alloc(dim_t numNodes, index_t* target, index_t unused);
void Finley_NodeMapping_free(Finley_NodeMapping*);
Finley_NodeMapping*  NodeMapping_getReference(Finley_NodeMapping *in );

#endif
