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

void Finley_ElementFile_scatter(int* index, Finley_ElementFile* in, Finley_ElementFile* out) {
   maybelong e,k,j;
   if (in!=NULL) {
     int NN_in=in->ReferenceElement->Type->numNodes;
     int NN_out=out->ReferenceElement->Type->numNodes;
     /*OMP */
     #pragma omp parallel for private(e,k,j) schedule(static)
     for (e=0;e<in->numElements;e++) {
        k=index[e];
        out->Id[k]=in->Id[e];
        out->Tag[k]=in->Tag[e];
        out->Color[k]=in->Color[e]+out->numColors;
        for(j=0;j<MIN(NN_out,NN_in);j++) out->Nodes[INDEX2(j,k,NN_out)]=in->Nodes[INDEX2(j,e,NN_in)];
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
