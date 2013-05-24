
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

void Finley_Assemble_AverageElementData(Finley_ElementFile* elements,
                                        escriptDataC* out, escriptDataC* in)
{
    Finley_resetError();
    if (!elements)
        return;

    double *wq;
    dim_t numQuad_in, numQuad_out;
    if (Finley_Assemble_reducedIntegrationOrder(in)) {
        numQuad_in=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
        wq=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->QuadWeights;
    } else {
        numQuad_in=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
        wq=elements->referenceElementSet->referenceElement->Parametrization->QuadWeights;
    }
    if (Finley_Assemble_reducedIntegrationOrder(out)) {
        numQuad_out=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad_out=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    // check out and in
    const dim_t numElements=elements->numElements;
    const dim_t numComps=getDataPointSize(out);

    if (numComps != getDataPointSize(in)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_AverageElementData: number of components of input and output data do not match.");
    } else if (!numSamplesEqual(in,numQuad_in,numElements)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!numSamplesEqual(out,numQuad_out,numElements)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!isExpanded(out)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_AverageElementData: expanded Data object is expected for output data.");
    } else {
        // now we can start
        if (isExpanded(in)) {
            double vol=0.;
            for (dim_t q=0; q< numQuad_in;++q) vol+=wq[q];
            const double volinv=1./vol;
            requireWrite(out);
#pragma omp parallel for
            for (dim_t n=0; n<numElements; n++) {
                const double *in_array = getSampleDataRO(in,n);
                double *out_array = getSampleDataRW(out,n);
                for (dim_t i=0; i<numComps; ++i) {
                    double rtmp=0.;
                    for (dim_t q=0; q<numQuad_in; ++q)
                        rtmp+=in_array[INDEX2(i,q,numComps)]*wq[q];
                    rtmp*=volinv;
                    for (dim_t q=0; q<numQuad_out; ++q)
                        out_array[INDEX2(i,q,numComps)]=rtmp;
                }
            }
        } else { // constant data
            const size_t numComps_size=numComps*sizeof(double);
            requireWrite(out);
#pragma omp parallel for
            for (dim_t n=0; n<numElements; n++) {
                const double *in_array = getSampleDataRO(in,n);
                double *out_array = getSampleDataRW(out,n);
                for (dim_t q=0; q<numQuad_out; q++)
                    memcpy(out_array+q*numComps, in_array, numComps_size);
            }
        }
    }
}

