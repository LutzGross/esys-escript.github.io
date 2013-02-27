
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

/*                                                                                                                     */
/* NodeMapping provides a mapping from the local nodes typically to the degrees of freedom,                            */
/*    the reduced degrees of freedom or the reduced node set                                                           */
/*                                                                                                                     */

#ifndef INC_DUDLEY_NODEMAPPING
#define INC_DUDLEY_NODEMAPPING

#include "esysUtils/Esys_MPI.h"

struct Dudley_NodeMapping {
    dim_t numNodes;		/* number of FEM nodes */
    index_t *target;		/* target[i] defines the target if FEM  node i =0,...,numNodes */
    index_t unused;		/* target[i]=unused defines that no target is defined for FEM  node i */
    dim_t numTargets;		/* number of targets */
    index_t *map;		/* maps the target nodes back to the FEM nodes: target[map[i]]=i */
    dim_t reference_counter;
};
typedef struct Dudley_NodeMapping Dudley_NodeMapping;

Dudley_NodeMapping *Dudley_NodeMapping_alloc(dim_t numNodes, index_t *target, index_t unused);
void Dudley_NodeMapping_free(Dudley_NodeMapping *);
Dudley_NodeMapping *Dudley_NodeMapping_getReference(Dudley_NodeMapping *in);

#endif
