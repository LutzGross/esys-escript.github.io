/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile */

/*   scatters the NodeFile in into NodeFile out using index[0:in->numNodes-1].  */
/*   index has to be between 0 and out->numNodes-1. */
/*   coloring is choosen for the worst case */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_scatter(maybelong* index, Finley_NodeFile* in, Finley_NodeFile* out) {
   maybelong i,k,j;
   if (in->Id!=NULL) {
     #pragma omp parallel for private(i,j,k) schedule(static)
     for (i=0;i<in->numNodes;i++) {
        k=index[i];
        out->Id[k]=in->Id[i];
        out->Tag[k]=in->Tag[i];
        out->degreeOfFreedom[k]=in->degreeOfFreedom[i];
        out->reducedDegreeOfFreedom[k]=in->reducedDegreeOfFreedom[i];
        out->toReduced[k]=in->toReduced[i];
        for(j=0;j<in->numDim;j++) out->Coordinates[INDEX2(j,k,in->numDim)]=in->Coordinates[INDEX2(j,i,in->numDim)];
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
