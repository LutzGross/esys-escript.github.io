
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

/*       Dudley: Mesh: ElementFile */

/*      set tags to newTag where mask>0 */

/****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"
#include "Assemble.h"
#include "Util.h"

namespace dudley {

void Dudley_ElementFile_setTags(Dudley_ElementFile* self, int newTag, const escript::Data* mask)
{
    const dim_t numElements = self->numElements;
    const int numQuad = dudley::Assemble_reducedIntegrationOrder(mask) ? 1 : (self->numDim + 1);

    if (1 != mask->getDataPointSize()) {
        throw DudleyException("Dudley_ElementFile_setTags: number of components of mask is 1.");
    } else if (!mask->numSamplesEqual(numQuad, numElements)) {
        throw DudleyException("Dudley_ElementFile_setTags: illegal number of samples of mask Data object");
    }

    if (mask->actsExpanded()) {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            const double* mask_array = mask->getSampleDataRO(n);
            if (mask_array[0] > 0)
                self->Tag[n] = newTag;
        }
    } else {
#pragma omp parallel for
        for (index_t n = 0; n < numElements; n++) {
            const double* mask_array = mask->getSampleDataRO(n);
            bool check = false;
            for (int q = 0; q < numQuad; q++)
                check = check || mask_array[q];
            if (check)
                self->Tag[n] = newTag;
        }
    }
    Dudley_ElementFile_setTagsInUse(self);
}

} // namespace dudley

