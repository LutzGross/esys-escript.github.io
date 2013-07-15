
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

  Assemblage routines: averages data between elements.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in)
{
    resetError();
    if (!elements)
        return;

    double *wq;
    int numQuad_in, numQuad_out;
    if (util::hasReducedIntegrationOrder(in)) {
        numQuad_in=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
        wq=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->QuadWeights;
    } else {
        numQuad_in=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
        wq=elements->referenceElementSet->referenceElement->Parametrization->QuadWeights;
    }
    if (util::hasReducedIntegrationOrder(out)) {
        numQuad_out=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad_out=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    // check out and in
    const int numElements=elements->numElements;
    const int numComps=out.getDataPointSize();

    if (numComps != in.getDataPointSize()) {
        setError(TYPE_ERROR, "Assemble_AverageElementData: number of components of input and output data do not match.");
    } else if (!in.numSamplesEqual(numQuad_in,numElements)) {
        setError(TYPE_ERROR, "Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!out.numSamplesEqual(numQuad_out,numElements)) {
        setError(TYPE_ERROR, "Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!out.actsExpanded()) {
        setError(TYPE_ERROR, "Assemble_AverageElementData: expanded Data object is expected for output data.");
    } else {
        escript::Data& _in(*const_cast<escript::Data*>(&in));
        if (in.actsExpanded()) {
            double vol=0.;
            for (int q=0; q<numQuad_in;++q) vol+=wq[q];
            const double volinv=1./vol;
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<numElements; n++) {
                const double *in_array = _in.getSampleDataRO(n);
                double *out_array = out.getSampleDataRW(n);
                for (int i=0; i<numComps; ++i) {
                    double rtmp=0.;
                    for (int q=0; q<numQuad_in; ++q)
                        rtmp+=in_array[INDEX2(i,q,numComps)]*wq[q];
                    rtmp*=volinv;
                    for (int q=0; q<numQuad_out; ++q)
                        out_array[INDEX2(i,q,numComps)]=rtmp;
                }
            }
        } else { // constant data
            const size_t numComps_size=numComps*sizeof(double);
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<numElements; n++) {
                const double *in_array = _in.getSampleDataRO(n);
                double *out_array = out.getSampleDataRW(n);
                for (int q=0; q<numQuad_out; q++)
                    memcpy(out_array+q*numComps, in_array, numComps_size);
            }
        }
    }
}

} // namespace finley

