/* $Id$ */
/**************************************************************/

/*   Finley: ElementFile                                                      */

/* copies element file in into element file out starting from offset          */
/* the elements offset to in->numElements+offset-1 in out will be overwritten */
                                                                                                                                                   
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "ElementFile.h"

/****************************************************************************/


void Finley_ElementFile_copyTable(index_t offset,Finley_ElementFile* out,index_t node_offset, index_t idOffset,Finley_ElementFile* in) {
    dim_t i,n;
    dim_t NN=out->ReferenceElement->Type->numNodes;
    if (in==NULL) return;
    /* check dimension and file size */
    if (out->ReferenceElement->Type->TypeId!=in->ReferenceElement->Type->TypeId) {
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"dimensions of element files don't match");
    }
    if (out->numElements<in->numElements+offset) {
        Finley_ErrorCode=MEMORY_ERROR;
        sprintf(Finley_ErrorMsg,"element table is too small.");
    }
    if (Finley_ErrorCode==NO_ERROR) {
       #pragma omp parallel for private(i,n) schedule(static)
       for(n=0;n<in->numElements;n++) {
          out->Id[offset+n]=in->Id[n]+idOffset;
          out->Tag[offset+n]=in->Tag[n];
          out->Color[offset+n]=in->Color[n]+out->maxColor+1;
          for(i=0;i<NN;i++) out->Nodes[INDEX2(i,offset+n,NN)]=in->Nodes[INDEX2(i,n,in->ReferenceElement->Type->numNodes)]+node_offset;
       }
       out->minColor=MIN(out->minColor,in->minColor+out->maxColor+1);
       out->maxColor=MAX(out->maxColor,in->maxColor+out->maxColor+1);
    }
}
/* 
* $Log$
* Revision 1.2  2005/07/08 04:07:49  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/06/30 01:53:55  gross
* a bug in coloring fixed
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:49  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
