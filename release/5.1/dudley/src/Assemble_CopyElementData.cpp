
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

namespace dudley {

using escript::ValueError;

template<typename Scalar>
void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in)
{
    if (!elements)
        return;

    dim_t numQuad = (hasReducedIntegrationOrder(in) ?
            QuadNums[elements->numDim][0] : QuadNums[elements->numDim][1]);

    // check out and in
    const dim_t numElements = elements->numElements;
    const int numComps = out.getDataPointSize();

    if (numComps != in.getDataPointSize()) {
        throw ValueError("Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!in.numSamplesEqual(numQuad, numElements)) {
        throw ValueError("Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (!out.numSamplesEqual(numQuad, numElements)) {
        throw ValueError("Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!out.actsExpanded()) {
        throw ValueError("Assemble_CopyElementData: expanded Data object is expected for output data.");
    } else if (in.isComplex() != out.isComplex()) {
        throw ValueError("Assemble_CopyElementData: complexity of input and output Data must match.");
    } else {
        const Scalar zero = static_cast<Scalar>(0);
        out.requireWrite();
        if (in.actsExpanded()) {
            const size_t len_size = numComps * numQuad * sizeof(Scalar);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++)
                memcpy(out.getSampleDataRW(n, zero),
                       in.getSampleDataRO(n, zero), len_size);
        } else {
            const size_t len_size = numComps * sizeof(Scalar);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const Scalar* in_array = in.getSampleDataRO(n, zero);
                Scalar* out_array = out.getSampleDataRW(n, zero);
                for (int q = 0; q < numQuad; q++)
                    memcpy(out_array + q * numComps, in_array, len_size);
            }
        }
    }
}

// instantiate our two supported versions
template void Assemble_CopyElementData<escript::DataTypes::real_t>(
                const ElementFile* elements, escript::Data& out,
                const escript::Data& in);
template void Assemble_CopyElementData<escript::DataTypes::cplx_t>(
                const ElementFile* elements, escript::Data& out,
                const escript::Data& in);

} // namespace dudley

