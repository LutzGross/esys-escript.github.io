
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

  Assemblage routines: copies data between elements.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in)
{
    resetError();
    if (!elements)
        return;

    int numQuad;
    if (util::hasReducedIntegrationOrder(in)) {
        numQuad=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    const int numElements=elements->numElements;
    const int numComps=out.getDataPointSize();

    if (numComps != in.getDataPointSize()) {
        setError(TYPE_ERROR,"Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!in.numSamplesEqual(numQuad,numElements)) {
        setError(TYPE_ERROR,"Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (!out.numSamplesEqual(numQuad,numElements)) {
        setError(TYPE_ERROR,"Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!out.actsExpanded()) {
        setError(TYPE_ERROR,"Assemble_CopyElementData: expanded Data object is expected for output data.");
    } else {
        if (in.actsExpanded()) {
            const size_t len_size=numComps*numQuad*sizeof(double);
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<numElements; n++) 
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(n), len_size);
        } else {
            const size_t len_size=numComps*sizeof(double);
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<numElements; n++) {
                const double *in_array = in.getSampleDataRO(n);
                double *out_array = out.getSampleDataRW(n);
                for (int q=0; q<numQuad; q++)
                    memcpy(out_array+q*numComps, in_array, len_size);
            }
        }
    }
}

} // namespace finley

