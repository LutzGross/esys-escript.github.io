/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: NodeFile */

/*   scatters the NodeFile in into NodeFile out using index[0:in->numNodes-1].  */
/*   index has to be between 0 and out->numNodes-1. */
/*   coloring is choosen for the worst case */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_scatter(index_t* index, Finley_NodeFile* in, Finley_NodeFile* out) {
   dim_t i,j;
   index_t k;
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
* Revision 1.3  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
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
