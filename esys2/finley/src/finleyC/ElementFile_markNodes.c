/* $Id$ */
/**************************************************************/

/*   Finley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_markNodes(int* mask,int offset,Finley_ElementFile* in,int useLinear) {
   int i,NN,NN2,e,color,*lin_node;
   if (in!=NULL) {
     int id[in->ReferenceElement->Type->numNodes];
     for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
     if (useLinear) {
        NN=in->LinearReferenceElement->Type->numNodes;
        lin_node=in->ReferenceElement->Type->linearNodes;
     } else {
        NN=in->ReferenceElement->Type->numNodes;
        lin_node=id;
     }
     NN2=in->ReferenceElement->Type->numNodes;
     #pragma omp parallel private(color)
     {
        for (color=0;color<in->numColors;color++) {
          #pragma omp for private(e,i) schedule(static)
          for (e=0;e<in->numElements;e++) {
            if (in->Color[e]==color) {
               for (i=0;i<NN;i++) mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
            }
          }
        }
        #pragma omp barrier
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
