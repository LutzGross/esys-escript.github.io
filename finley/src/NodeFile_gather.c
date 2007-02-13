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

/*   Finley: Mesh: NodeFile                                   */

/*   gathers the NodeFile out from the  NodeFile in using index[0:out->numNodes-1].  */
/*   index has to be between 0 and in->numNodes-1. */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_gather(int* index, Finley_NodeFile* in, Finley_NodeFile* out) {
   dim_t i,j;
   index_t k;
   if (in->Id!=NULL) {
     #pragma omp parallel for private(i,j,k) schedule(static)
     for (i=0;i<out->numNodes;i++) {
        k=index[i];
#ifdef PASO_MPI
        out->Dom[i]=in->Dom[k];
#endif
        out->Id[i]=in->Id[k];
        out->Tag[i]=in->Tag[k];
        out->degreeOfFreedom[i]=in->degreeOfFreedom[k];
        for(j=0;j<in->numDim;j++) out->Coordinates[INDEX2(j,i,in->numDim)]=in->Coordinates[INDEX2(j,k,in->numDim)];
     }
     out->isPrepared=FINLEY_UNPREPARED;
   }
}
