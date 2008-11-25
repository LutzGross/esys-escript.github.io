
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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


void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in) {
    dim_t n,q, numElements, numQuad;
    __const double *in_array;
    double *out_array;
    dim_t numComps=getDataPointSize(out);
    size_t len_size;

    Finley_resetError();
    if( elements == NULL )
    {
       return;
    }

    numElements=elements->numElements;
    if (Finley_Assemble_reducedIntegrationOrder(in)) {
       numQuad=elements->ReferenceElementReducedOrder->numQuadNodes;
    } else {
       numQuad=elements->ReferenceElement->numQuadNodes;
    }

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: number of components of input and output Data do not match.");
    } else if (!numSamplesEqual(in,numQuad,numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: illegal number of samples of input Data object");
    } else if (!numSamplesEqual(out,numQuad,numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: illegal number of samples of output Data object");
    } else if (!isExpanded(out)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyElementData: expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Finley_noError()) {
         if (isExpanded(in)) {
             len_size=numComps*numQuad*sizeof(double);
             # pragma omp parallel for private(n) schedule(static)
             for (n=0;n<numElements;n++) 
                 memcpy(getSampleDataRW(out,n),getSampleDataRO(in,n), len_size);
         } else {
             len_size=numComps*sizeof(double);
             # pragma omp parallel for private(q,n,out_array,in_array) schedule(static)
             for (n=0;n<numElements;n++) {
                 in_array=getSampleDataRO(in,n);
                 out_array=getSampleDataRW(out,n);
                 for (q=0;q<numQuad;q++) memcpy(out_array+q*numComps,in_array,len_size);
             }
         }
    }
    return;
}
