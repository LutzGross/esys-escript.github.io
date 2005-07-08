/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile                                   */

/*   gathers the NodeFile out from the  NodeFile in using index[0:out->numNodes-1].  */
/*   index has to be between 0 and in->numNodes-1. */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_gather(int* index, Finley_NodeFile* in, Finley_NodeFile* out) {
   dim_t i,j;
   index_t k;
   if (in->Id!=NULL) {
     #pragma omp parallel for private(i,j,k) schedule(static)
     for (i=0;i<out->numNodes;i++) {
        k=index[i];
        out->Id[i]=in->Id[k];
        out->Tag[i]=in->Tag[k];
        out->degreeOfFreedom[i]=in->degreeOfFreedom[k];
        out->reducedDegreeOfFreedom[i]=in->reducedDegreeOfFreedom[k];
        out->toReduced[i]=in->toReduced[k];
        for(j=0;j<in->numDim;j++) out->Coordinates[INDEX2(j,i,in->numDim)]=in->Coordinates[INDEX2(j,k,in->numDim)];
     }
   }
}
/* 
* $Log$
* Revision 1.2  2005/07/08 04:07:55  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:54  gross
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
