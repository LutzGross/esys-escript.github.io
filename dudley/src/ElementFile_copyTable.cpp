
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

/*****************************************************************************

  Dudley: ElementFile

  copies element file in into element file out starting from offset
  the elements offset to in->numElements+offset-1 in out will be overwritten

*****************************************************************************/
#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"

namespace dudley {

void Dudley_ElementFile_copyTable(index_t offset, Dudley_ElementFile* out,
        index_t node_offset, index_t idOffset, Dudley_ElementFile* in)
{
    dim_t i, n;
    dim_t NN, NN_in;
    if (in == NULL)
        return;
    NN = out->numNodes;
    NN_in = in->numNodes;
    if (NN_in > NN)
    {
        throw DudleyException("Dudley_ElementFile_copyTable: dimensions of element files don't match.");
    }
    if (out->MPIInfo->comm != in->MPIInfo->comm)
    {
        throw DudleyException("Dudley_ElementFile_copyTable: MPI communicators of element files don't match.");
    }

#pragma omp parallel for private(i,n) schedule(static)
    for (n = 0; n < in->numElements; n++)
    {
        out->Owner[offset + n] = out->Owner[n];
        out->Id[offset + n] = in->Id[n] + idOffset;
        out->Tag[offset + n] = in->Tag[n];
        for (i = 0; i < NN; i++)
            out->Nodes[INDEX2(i, offset + n, NN)] = in->Nodes[INDEX2(i, n, NN_in)] + node_offset;
    }
}

} // namespace dudley

