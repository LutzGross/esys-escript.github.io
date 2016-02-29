
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

/*   Dudley: ElementFile */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"

namespace dudley {

void Dudley_ElementFile_relableNodes(index_t* newNode, index_t offset, Dudley_ElementFile* in)
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

} // namespace dudley

