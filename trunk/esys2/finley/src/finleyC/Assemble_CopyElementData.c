/* $Id$ */

/**************************************************************/

/*    assemblage routines: copies data between elements       */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Util.h"
#include "Finley.h"
#include "Assemble.h"
#include "ElementFile.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/******************************************************************************************************/


void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in) {
    if (elements==NULL) return;
    dim_t n,q;
    dim_t numElements=elements->numElements;
    dim_t numQuad=elements->ReferenceElement->numQuadNodes;
    dim_t numComps=getDataPointSize(out);
    double *in_array,*out_array;

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"number of components of input and output Data do not match.");
    } else if (numSamplesEqual(out,numQuad,numElements)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of output Data object");
    } else if (numSamplesEqual(in,numQuad,numElements)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of input Data object");
    } else if (!isExpanded(out)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Finley_ErrorCode==NO_ERROR) {
         if (isExpanded(in)) {
             # pragma omp parallel for private(n) schedule(static)
             for (n=0;n<numElements;n++) 
                 Finley_copyDouble(numComps*numQuad,getSampleData(in,n),getSampleData(out,n));
         } else {
             # pragma omp parallel for private(q,n,out_array,in_array) schedule(static)
             for (n=0;n<numElements;n++) {
                 in_array=getSampleData(in,n);
                 out_array=getSampleData(out,n);
                 for (q=0;q<numQuad;q++) Finley_copyDouble(numComps,in_array,out_array+q*numComps);
             }
         }
    }
    return;
}
/*
 * $Log$
 * Revision 1.2  2005/07/08 04:07:45  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:46  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:56  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
