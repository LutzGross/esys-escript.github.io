/* $Id$ */
/**************************************************************/
/*                                                                                            */
/*   Finley: ElementFile                                                                      */
/*                                                                                            */
/*   scatter the ElementFile in into the  ElementFile out using index[0:out->numElements-1].  */
/*   index has to be between 0 and in->numElements-1.                                         */
/*   a conservative assumtion on the coloring is made                                         */
/*                                                                                            */
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_scatter(index_t* index, Finley_ElementFile* in, Finley_ElementFile* out) {
   index_t k;
   dim_t e,j;
   if (in!=NULL) {
     dim_t NN_in=in->ReferenceElement->Type->numNodes;
     dim_t NN_out=out->ReferenceElement->Type->numNodes;
     /*OMP */
     #pragma omp parallel for private(e,k,j) schedule(static)
     for (e=0;e<in->numElements;e++) {
        k=index[e];
        out->Id[k]=in->Id[e];
        out->Tag[k]=in->Tag[e];
        out->Color[k]=in->Color[e]+out->maxColor+1;
        for(j=0;j<MIN(NN_out,NN_in);j++) out->Nodes[INDEX2(j,k,NN_out)]=in->Nodes[INDEX2(j,e,NN_in)];
     }
     out->minColor=MIN(out->minColor,in->minColor+out->maxColor+1);
     out->maxColor=MAX(out->maxColor,in->maxColor+out->maxColor+1);
   }
}
/* 
* $Log$
* Revision 1.2  2005/07/08 04:07:50  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/06/30 01:53:55  gross
* a bug in coloring fixed
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:50  gross
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
