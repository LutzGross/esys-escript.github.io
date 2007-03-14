/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/*    assemblage routines: copies data between elements       */

/**************************************************************/

/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/******************************************************************************************************/


void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in) {
    dim_t n,q, numElements, numQuad;
    dim_t numComps=getDataPointSize(out);
    double *in_array,*out_array;
    Finley_resetError();

    if (elements==NULL) return;
    numElements=elements->numElements;
    numQuad=elements->ReferenceElement->numQuadNodes;

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_setError(TYPE_ERROR,"__FILE__: number of components of input and output Data do not match.");
    } else if (!numSamplesEqual(in,numQuad,numElements)) {
       Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of input Data object");
    } else if (!numSamplesEqual(out,numQuad,numElements)) {
       Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of output Data object");
    } else if (!isExpanded(out)) {
       Finley_setError(TYPE_ERROR,"__FILE__: expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Finley_noError()) {
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
 * Revision 1.4  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.3  2005/08/12 01:45:42  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.2.2.2  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2.2.1  2005/08/02 05:29:11  gross
 * bug in finley/src/Assemble_CopyElementData fixed
 *
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
