
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Assemblage routines: averages data between elements.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace finley {

template<typename Scalar>
void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in)
{
    if (!elements)
        return;

    const double* wq;
    int numQuad_in, numQuad_out;
    if (util::hasReducedIntegrationOrder(in)) {
        numQuad_in = elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
        wq = &elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->QuadWeights[0];
    } else {
        numQuad_in = elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
        wq = &elements->referenceElementSet->referenceElement->Parametrization->QuadWeights[0];
    }
    if (util::hasReducedIntegrationOrder(out)) {
        numQuad_out = elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad_out = elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    // check out and in
    const dim_t numElements = elements->numElements;
    const int numComps = out.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);

    if (numComps != in.getDataPointSize()) {
        throw escript::ValueError("Assemble_AverageElementData: number of components of input and output data do not match.");
    } else if (!in.numSamplesEqual(numQuad_in,numElements)) {
        throw escript::ValueError("Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!out.numSamplesEqual(numQuad_out,numElements)) {
        throw escript::ValueError("Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!out.actsExpanded()) {
        throw escript::ValueError("Assemble_AverageElementData: expanded Data object is expected for output data.");
    } else if (in.isComplex() != out.isComplex()) {
        throw escript::ValueError("Assemble_AverageElementData: complexity of input and output data must match.");
    } else {
        out.requireWrite();
        if (in.actsExpanded()) {
            double vol = 0.;
            for (int q = 0; q < numQuad_in; ++q)
                vol += wq[q];
            const double volinv = 1./vol;
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const Scalar* in_array = in.getSampleDataRO(n, zero);
                Scalar* out_array = out.getSampleDataRW(n, zero);
                for (int i = 0; i < numComps; ++i) {
                    Scalar rtmp = zero;
                    for (int q = 0; q < numQuad_in; ++q)
                        rtmp += in_array[INDEX2(i,q,numComps)] * wq[q];
                    rtmp *= volinv;
                    for (int q = 0; q < numQuad_out; ++q)
                        out_array[INDEX2(i,q,numComps)] = rtmp;
                }
            }
        } else { // constant data
            const size_t numComps_size = numComps * sizeof(Scalar);
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const Scalar* in_array = in.getSampleDataRO(n, zero);
                Scalar* out_array = out.getSampleDataRW(n, zero);
                for (int q = 0; q < numQuad_out; q++)
                    memcpy(out_array+q*numComps, in_array, numComps_size);
            }
        }
    }
}

// instantiate our two supported versions
template void Assemble_AverageElementData<escript::DataTypes::real_t>(
                const ElementFile* elements,
                escript::Data& out, const escript::Data& in);
template void Assemble_AverageElementData<escript::DataTypes::cplx_t>(
                const ElementFile* elements,
                escript::Data& out, const escript::Data& in);

} // namespace finley

