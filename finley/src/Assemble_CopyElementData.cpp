
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Assemblage routines: copies data between elements.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

template<typename Scalar>
void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in)
{
    if (!elements)
        return;

    int numQuad_out, numQuad_in;
    if (util::hasReducedIntegrationOrder(out)) {
        numQuad_out = elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad_out = elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }
    if (util::hasReducedIntegrationOrder(in)) {
        numQuad_in = elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad_in = elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    const dim_t numElements = elements->numElements;
    const int numComps = out.getDataPointSize();

    if (numComps != in.getDataPointSize()) {
        throw escript::ValueError("Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!out.actsExpanded()) {
        throw escript::ValueError("Assemble_CopyElementData: expanded Data object is expected for output data.");
    } else if (!out.numSamplesEqual(numQuad_out,numElements)) {
        throw escript::ValueError("Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!in.numSamplesEqual(numQuad_in,numElements)) {
        throw escript::ValueError("Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (in.isComplex() != out.isComplex()) {
        throw escript::ValueError("Assemble_CopyElementData: complexity of input and output Data must match.");
    } else {
        const Scalar zero = static_cast<Scalar>(0);
        if (numQuad_in == 1) {
            const size_t len_size = numComps*sizeof(Scalar);
            out.requireWrite();
#pragma omp parallel for
            for (index_t n = 0; n < numElements; n++) {
                const Scalar* in_array = in.getSampleDataRO(n, zero);
                Scalar* out_array = out.getSampleDataRW(n, zero);
                for (int q = 0; q < numQuad_out; q++)
                    memcpy(out_array+q*numComps, in_array, len_size);
            }
        } else if (numQuad_in == numQuad_out) {
            out.requireWrite();
            if (in.actsExpanded()) {
                const size_t len_size = numComps*numQuad_in*sizeof(Scalar);
#ifndef _WIN32 // TODO: why fatal error C1001: An internal error has occurred in the compiler.?
#pragma omp parallel for
#endif
                for (index_t n = 0; n < numElements; n++) 
                    memcpy(out.getSampleDataRW(n, zero),
                            in.getSampleDataRO(n, zero), len_size);
            } else {
                const size_t len_size = numComps*sizeof(Scalar);
#pragma omp parallel for
                for (index_t n = 0; n < numElements; n++) {
                    const Scalar* in_array = in.getSampleDataRO(n, zero);
                    Scalar* out_array = out.getSampleDataRW(n, zero);
                    for (int q = 0; q < numQuad_in; q++)
                        memcpy(out_array+q*numComps, in_array, len_size);
                }
            }
        } else {
            throw escript::ValueError("Assemble_CopyElementData: unable to process given number of data points.");
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

} // namespace finley

