
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

/*   Dudley: Mesh: NodeFile */

/* copies node file in into node file out starting from offset          */
/* the nodes offset to in->numNodes+offset-1 in out will be overwritten */

/************************************************************************************/

#include "NodeFile.h"

/************************************************************************************/

void Dudley_NodeFile_copyTable(int offset, Dudley_NodeFile * out, int idOffset, int dofOffset, Dudley_NodeFile * in)
{
    int i, n;
    /* check dimension and file size */
    if (out->numDim != in->numDim)
    {
	Dudley_setError(TYPE_ERROR, "Dudley_NodeFile_copyTable: dimensions of node files don't match");
    }
    if (out->numNodes < in->numNodes + offset)
    {
	Dudley_setError(MEMORY_ERROR, "Dudley_NodeFile_copyTable: node table is too small.");
    }
    if (Dudley_noError())
    {
#pragma omp parallel for private(i,n) schedule(static)
	for (n = 0; n < in->numNodes; n++)
	{
	    out->Id[offset + n] = in->Id[n] + idOffset;
	    out->Tag[offset + n] = in->Tag[n];
	    out->globalDegreesOfFreedom[offset + n] = in->globalDegreesOfFreedom[n] + dofOffset;
	    for (i = 0; i < out->numDim; i++)
		out->Coordinates[INDEX2(i, offset + n, out->numDim)] = in->Coordinates[INDEX2(i, n, in->numDim)];
	}
    }
}
