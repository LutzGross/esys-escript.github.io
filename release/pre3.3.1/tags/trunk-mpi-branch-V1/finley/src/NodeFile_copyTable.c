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

/*   Finley: Mesh: NodeFile */

/* copies node file in into node file out starting from offset          */
/* the nodes offset to in->numNodes+offset-1 in out will be overwritten */


/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_copyTable(int offset,Finley_NodeFile* out,int idOffset,int dofOffset,Finley_NodeFile* in) {
    int i,n;
    /* check dimension and file size */
    if (out->numDim!=in->numDim) {
        Finley_setError(TYPE_ERROR,"__FILE__: dimensions of node files don't match");
    }
    if (out->numNodes<in->numNodes+offset) {
        Finley_setError(MEMORY_ERROR,"__FILE__: node table is too small.");
    }
    if (Finley_noError()) {
       #pragma omp parallel for private(i,n) schedule(static)
       for(n=0;n<in->numNodes;n++) {
          out->Id[offset+n]=in->Id[n]+idOffset;
          out->Tag[offset+n]=in->Tag[n];
          out->degreeOfFreedom[offset+n]=in->degreeOfFreedom[n]+dofOffset;
          out->reducedDegreeOfFreedom[offset+n]=in->reducedDegreeOfFreedom[n]+dofOffset;
          out->toReduced[offset+n]=in->toReduced[n]+dofOffset;
#ifdef PASO_MPI
					out->Dom[offset+n]=in->Dom[n];
#endif
          for(i=0;i<out->numDim;i++) out->Coordinates[INDEX2(i,offset+n,out->numDim)]=in->Coordinates[INDEX2(i,n,in->numDim)];
       }
    }
}
/* 
* $Log$
* Revision 1.2  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.1.1.1.6.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
