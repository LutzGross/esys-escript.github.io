
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

/****************************************************************************/

/*   Dudley: Mesh: NodeFile */

/* copies node file in into node file out starting from offset          */
/* the nodes offset to in->numNodes+offset-1 in out will be overwritten */

/****************************************************************************/

#include "NodeFile.h"

namespace dudley {

void Dudley_NodeFile_copyTable(int offset, Dudley_NodeFile * out, int idOffset, int dofOffset, Dudley_NodeFile * in)
{
    /* check dimension and file size */
    if (out->numDim != in->numDim)
        throw DudleyException("Dudley_NodeFile_copyTable: dimensions of node files don't match");

    if (out->numNodes < in->numNodes + offset)
        throw DudleyException("Dudley_NodeFile_copyTable: node table is too small.");

    int i, n;
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

} // namespace dudley

