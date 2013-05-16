
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

void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,
                                     escriptDataC* out, escriptDataC* in)
{
    Finley_resetError();
    if (!elements)
        return;

    dim_t numQuad;
    if (Finley_Assemble_reducedIntegrationOrder(in)) {
        numQuad=elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
    } else {
        numQuad=elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
    }

    // check out and in
    const dim_t numElements=elements->numElements;
    const dim_t numComps=getDataPointSize(out);

    if (numComps!=getDataPointSize(in)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!numSamplesEqual(in,numQuad,numElements)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (!numSamplesEqual(out,numQuad,numElements)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!isExpanded(out)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: expanded Data object is expected for output data.");
    } else {
        // now we can start
        if (isExpanded(in)) {
            const size_t len_size=numComps*numQuad*sizeof(double);
            requireWrite(out);
#pragma omp parallel for schedule(static)
            for (dim_t n=0; n<numElements; n++) 
                memcpy(getSampleDataRW(out,n), getSampleDataRO(in,n), len_size);
         } else {
            const size_t len_size=numComps*sizeof(double);
            requireWrite(out);
#pragma omp parallel for schedule(static)
            for (dim_t n=0; n<numElements; n++) {
                const double *in_array = getSampleDataRO(in,n);
                double *out_array = getSampleDataRW(out,n);
                for (dim_t q=0; q<numQuad; q++)
                    memcpy(out_array+q*numComps, in_array, len_size);
            }
         }
    }
    return;
}
