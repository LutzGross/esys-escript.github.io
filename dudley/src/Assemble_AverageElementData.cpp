
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

  Assemblage routines: averages data

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

namespace dudley {

void Assemble_AverageElementData(Dudley_ElementFile* elements, escript::Data* out, const escript::Data* in)
{
    Dudley_resetError();
    if (!elements)
        return;

    double wq;
    int numQuad_in, numQuad_out;
    if (Assemble_reducedIntegrationOrder(in)) {
        numQuad_in = QuadNums[elements->numDim][0];
        wq = QuadWeight[elements->numDim][0];
    } else {
        numQuad_in = QuadNums[elements->numDim][1];
        wq = QuadWeight[elements->numDim][1];
    }
    if (Assemble_reducedIntegrationOrder(out)) {
        numQuad_out = QuadNums[elements->numDim][0];
    } else {
        numQuad_out = QuadNums[elements->numDim][1];
    }

    const int numComps = out->getDataPointSize();
    const dim_t numElements = elements->numElements;

    // check out and in
    if (numComps != in->getDataPointSize()) {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_AverageElementData: number of components of input and output Data do not match.");
    } else if (!in->numSamplesEqual(numQuad_in, numElements)) {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!out->numSamplesEqual(numQuad_out, numElements)) {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!out->actsExpanded()) {
        Dudley_setError(TYPE_ERROR,
                        "Assemble_AverageElementData: expanded Data object is expected for output data.");
    } else {
        out->requireWrite();
        if (in->actsExpanded()) {
            const double vol = wq * numQuad_in;
            const double volinv = 1. / vol;
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const double* in_array = in->getSampleDataRO(n);
                double* out_array = out->getSampleDataRW(n);
                for (int i = 0; i < numComps; ++i) {
                    double rtmp = 0.;
                    for (int q = 0; q < numQuad_in; ++q)
                        rtmp += in_array[INDEX2(i, q, numComps)] * wq;
                    rtmp *= volinv;
                    for (int q = 0; q < numQuad_out; ++q)
                        out_array[INDEX2(i, q, numComps)] = rtmp;
                }
            }
        } else { // constant data
            const size_t numComps_size = numComps * sizeof(double);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const double* in_array = in->getSampleDataRO(n);
                double* out_array = out->getSampleDataRW(n);
                for (int q = 0; q < numQuad_out; q++)
                    memcpy(out_array + q * numComps, in_array, numComps_size);
            }
        }
    }
}

} // namespace dudley

