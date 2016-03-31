
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************/

/*   Dudley: Mesh: NodeFile */

/*  set tags to newTag where mask>0 */

/****************************************************************************/

#include "NodeFile.h"
#include "Util.h"

namespace dudley {

void Dudley_NodeFile_setTags(Dudley_NodeFile * self, const int newTag, const escript::Data* mask)
{
    if (self == NULL)
        return;
    dim_t n;
    const double *mask_array;

    dim_t numNodes = self->numNodes;
    if (1 != mask->getDataPointSize()) {
        throw DudleyException("Dudley_NodeFile_setTags: number of components of mask is 1.");
    } else if (mask->getNumDataPointsPerSample() != 1 ||
            mask->getNumSamples() != numNodes) {
        throw DudleyException("Dudley_NodeFile_setTags: illegal number of samples of mask Data object");
    }

#pragma omp parallel for private(n,mask_array)
    for (n = 0; n < numNodes; n++)
    {
        mask_array = mask->getSampleDataRO(n);
        if (mask_array[0] > 0)
            self->Tag[n] = newTag;
    }
    Dudley_NodeFile_setTagsInUse(self);
}

} // namespace dudley

