/* $Id$ */
/**************************************************************/

/*   Finley: ElementFile */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_relableNodes(int* newNode,int offset,Finley_ElementFile* in) {
   maybelong i,j,NN;
   if (in!=NULL) {
     NN=in->ReferenceElement->Type->numNodes;
     #pragma omp parallel for private(j,i) schedule(static)
     for (j=0;j<in->numElements;j++) {
       for (i=0;i<NN;i++) {
         in->Nodes[INDEX2(i,j,NN)]=newNode[in->Nodes[INDEX2(i,j,NN)]-offset];
       }
     }
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
