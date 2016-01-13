
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assemblage routines: copies data between elements       */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/******************************************************************************************************/


void Finley_Assemble_AverageElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in) {
    dim_t n,q, numElements, numQuad_in, numQuad_out, i;
    __const double *in_array;
    double *out_array, vol, volinv, *wq;
    register double rtmp;
    dim_t numComps=getDataPointSize(out);
    size_t numComps_size;

    Finley_resetError();
    if( elements == NULL )
    {
       return;
    }

    numElements=elements->numElements;
    if (Finley_Assemble_reducedIntegrationOrder(in)) {
       numQuad_in=elements->ReferenceElementReducedOrder->numQuadNodes;
       wq=elements->ReferenceElementReducedOrder->QuadWeights;
    } else {
       numQuad_in=elements->ReferenceElement->numQuadNodes;
       wq=elements->ReferenceElement->QuadWeights;
    }
    if (Finley_Assemble_reducedIntegrationOrder(out)) {
       numQuad_out=elements->ReferenceElementReducedOrder->numQuadNodes;
    } else {
       numQuad_out=elements->ReferenceElement->numQuadNodes;
    }

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_AverageElementData: number of components of input and output Data do not match.");
    } else if (!numSamplesEqual(in,numQuad_in,numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_AverageElementData: illegal number of samples of input Data object");
    } else if (!numSamplesEqual(out,numQuad_out,numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_AverageElementData: illegal number of samples of output Data object");
    } else if (!isExpanded(out)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_AverageElementData: expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Finley_noError()) {
         if (isExpanded(in)) {
	     void* buffer=allocSampleBuffer(in);
             vol=0;
             for (q=0; q< numQuad_in;++q) vol+=wq[q];
             volinv=1./vol;
	     requireWrite(out);
	     #pragma omp parallel private(n, i, rtmp, q, in_array, out_array)
	     {
               # pragma omp for schedule(static)
               for (n=0;n<numElements;n++) {
                 in_array=getSampleDataRO(in,n,buffer);
                 out_array=getSampleDataRW(out,n);
                 for (i=0; i<numComps; ++i) {
                     rtmp=0;
                     for (q=0; q< numQuad_in;++q) rtmp+=in_array[INDEX2(i,q,numComps)]*wq[q];
                     rtmp*=volinv;
                     for (q=0; q< numQuad_out;++q) out_array[INDEX2(i,q,numComps)]=rtmp;
                 }
               }
	     }
	     freeSampleBuffer(buffer);
         } else {
	     void* buffer=allocSampleBuffer(in);
             numComps_size=numComps*sizeof(double);
	     requireWrite(out);
	     #pragma omp parallel private(q,n,out_array,in_array)
	     {
               # pragma omp for schedule(static)
               for (n=0;n<numElements;n++) {
                 in_array=getSampleDataRO(in,n,buffer);
                 out_array=getSampleDataRW(out,n);
                 for (q=0;q<numQuad_out;q++) memcpy(out_array+q*numComps,in_array,numComps_size);
               }
	     }
	     freeSampleBuffer(buffer);
         }
    }
    return;
}
