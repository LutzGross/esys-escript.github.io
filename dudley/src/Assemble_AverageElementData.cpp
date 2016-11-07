
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

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in)
{
    if (!elements)
        return;

    if (out.isComplex() || in.isComplex())
    {
        throw DudleyException("Assemble_AverageElementData: complex arguments not supported.");      
    }
    double wq;
    int numQuad_in, numQuad_out;
    if (hasReducedIntegrationOrder(in)) {
        numQuad_in = QuadNums[elements->numDim][0];
        wq = QuadWeight[elements->numDim][0];
    } else {
        numQuad_in = QuadNums[elements->numDim][1];
        wq = QuadWeight[elements->numDim][1];
    }
    if (hasReducedIntegrationOrder(out)) {
        numQuad_out = QuadNums[elements->numDim][0];
    } else {
        numQuad_out = QuadNums[elements->numDim][1];
    }

    // check out and in
    const dim_t numElements = elements->numElements;
    const int numComps = out.getDataPointSize();

    if (numComps != in.getDataPointSize()) {
        throw DudleyException("Assemble_AverageElementData: number of components of input and output Data do not match.");
    } else if (!in.numSamplesEqual(numQuad_in, numElements)) {
        throw DudleyException("Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!out.numSamplesEqual(numQuad_out, numElements)) {
        throw DudleyException("Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!out.actsExpanded()) {
        throw DudleyException("Assemble_AverageElementData: expanded Data object is expected for output data.");
    } else {
        out.requireWrite();
        if (in.actsExpanded()) {
            const double vol = wq * numQuad_in;
            const double volinv = 1. / vol;
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const double* in_array = in.getSampleDataRO(n, static_cast<escript::DataTypes::real_t>(0));
                double* out_array = out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0));
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
                const double* in_array = in.getSampleDataRO(n, static_cast<escript::DataTypes::real_t>(0));
                double* out_array = out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0));
                for (int q = 0; q < numQuad_out; q++)
                    memcpy(out_array + q * numComps, in_array, numComps_size);
            }
        }
    }
}

} // namespace dudley

