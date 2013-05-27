
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


/****************************************************************************

  Finley Mesh: ElementFile

  Sets tags to newTag where mask>0

*****************************************************************************/

#include "ElementFile.h"
#include "Util.h"
#include "Assemble.h"

void Finley_ElementFile_setTags(Finley_ElementFile* self, const int newTag,
                                escriptDataC* mask)
{
    Finley_resetError();
    if (!self)
        return;

    const int numElements=self->numElements;
    const int numQuad=Finley_ReferenceElementSet_borrowReferenceElement(
        self->referenceElementSet, Finley_Assemble_reducedIntegrationOrder(mask))
        ->Parametrization->numQuadNodes;

    if (1 != getDataPointSize(mask)) {
        Finley_setError(TYPE_ERROR, "Finley_ElementFile_setTags: number of components of mask must be 1.");
        return;
    } else if (!numSamplesEqual(mask, numQuad, numElements)) {
        Finley_setError(TYPE_ERROR, "Finley_ElementFile_setTags: illegal number of samples of mask Data object");
        return;
    }

    if (isExpanded(mask)) {
#pragma omp parallel for
        for (int n=0; n<numElements; n++) {
            const double *mask_array=getSampleDataRO(mask,n);
            if (mask_array[0]>0)
                self->Tag[n]=newTag;
        }
    } else {
#pragma omp parallel for
        for (int n=0; n<numElements; n++) {
            const double *mask_array=getSampleDataRO(mask,n);
            bool check=false;
            for (int q=0; q<numQuad; q++)
                check = (check || mask_array[q]);
            if (check)
                self->Tag[n]=newTag;
        }
    }
    Finley_ElementFile_setTagsInUse(self);
}

