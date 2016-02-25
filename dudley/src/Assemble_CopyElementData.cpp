
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

/****************************************************************************

  Assemblage routines: copies data between elements

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

namespace dudley {

void Assemble_CopyElementData(Dudley_ElementFile* elements, escript::Data* out, const escript::Data* in)
{
    Dudley_resetError();
    if (!elements)
        return;

    dim_t numQuad;
    if (Assemble_reducedIntegrationOrder(in)) {
        numQuad = QuadNums[elements->numDim][0];
    } else {
        numQuad = QuadNums[elements->numDim][1];
    }

    const dim_t numElements = elements->numElements;
    const int numComps = out->getDataPointSize();

    // check out and in
    if (numComps != in->getDataPointSize())
    {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!in->numSamplesEqual(numQuad, numElements)) {
        Dudley_setError(TYPE_ERROR, "Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (!out->numSamplesEqual(numQuad, numElements)) {
        Dudley_setError(TYPE_ERROR, "Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!out->actsExpanded()) {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_CopyElementData: expanded Data object is expected for output data.");
    } else {
        out->requireWrite();
        if (in->actsExpanded()) {
            const size_t len_size = numComps * numQuad * sizeof(double);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++)
                memcpy(out->getSampleDataRW(n), in->getSampleDataRO(n), len_size);
        } else {
            const size_t len_size = numComps * sizeof(double);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const double* in_array = in->getSampleDataRO(n);
                double* out_array = out->getSampleDataRW(n);
                for (int q = 0; q < numQuad; q++)
                    memcpy(out_array + q * numComps, in_array, len_size);
            }
        }
    }
}

} // namespace dudley

