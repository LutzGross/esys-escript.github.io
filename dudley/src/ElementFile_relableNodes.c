
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

/*   Dudley: ElementFile */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/************************************************************************************/

#include "ElementFile.h"

/************************************************************************************/

void Dudley_ElementFile_relableNodes(index_t * newNode, index_t offset, Dudley_ElementFile * in)
{
    dim_t i, j, NN;

    if (in != NULL)
    {
	NN = in->numNodes;
#pragma omp parallel for private(j,i) schedule(static)
	for (j = 0; j < in->numElements; j++)
	{
	    for (i = 0; i < NN; i++)
	    {
		in->Nodes[INDEX2(i, j, NN)] = newNode[in->Nodes[INDEX2(i, j, NN)] - offset];
	    }
	}
    }
}
