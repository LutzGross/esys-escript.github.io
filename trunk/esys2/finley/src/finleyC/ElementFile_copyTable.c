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


void Finley_ElementFile_copyTable(int offset,Finley_ElementFile* out,int node_offset, int idOffset,Finley_ElementFile* in) {
    maybelong i,n;
    maybelong NN=out->ReferenceElement->Type->numNodes;
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
          out->Color[offset+n]=in->Color[n]+out->numColors;
          for(i=0;i<NN;i++) out->Nodes[INDEX2(i,offset+n,NN)]=in->Nodes[INDEX2(i,n,in->ReferenceElement->Type->numNodes)]+node_offset;
       }
       out->numColors+=in->numColors;
    }
}
/* 
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
