/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile */

/* copies node file in into node file out starting from offset          */
/* the nodes offset to in->numNodes+offset-1 in out will be overwritten */


/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_copyTable(int offset,Finley_NodeFile* out,int idOffset,int dofOffset,Finley_NodeFile* in) {
    int i,n;
    /* check dimension and file size */
    if (out->numDim!=in->numDim) {
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"dimensions of node files don't match");
    }
    if (out->numNodes<in->numNodes+offset) {
        Finley_ErrorCode=MEMORY_ERROR;
        sprintf(Finley_ErrorMsg,"node table is too small.");
    }
    if (Finley_ErrorCode==NO_ERROR) {
       #pragma omp parallel for private(i,n) schedule(static)
       for(n=0;n<in->numNodes;n++) {
          out->Id[offset+n]=in->Id[n]+idOffset;
          out->Tag[offset+n]=in->Tag[n];
          out->degreeOfFreedom[offset+n]=in->degreeOfFreedom[n]+dofOffset;
          out->reducedDegreeOfFreedom[offset+n]=in->reducedDegreeOfFreedom[n]+dofOffset;
          out->toReduced[offset+n]=in->toReduced[n]+dofOffset;
          for(i=0;i<out->numDim;i++) out->Coordinates[INDEX2(i,offset+n,out->numDim)]=in->Coordinates[INDEX2(i,n,in->numDim)];
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
