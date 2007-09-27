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
#ifdef PASO_MPI
        out->Dom[k]=in->Dom[i];
#endif
        out->Id[k]=in->Id[i];
        out->Tag[k]=in->Tag[i];
        out->degreeOfFreedom[k]=in->degreeOfFreedom[i];
        for(j=0;j<in->numDim;j++) out->Coordinates[INDEX2(j,k,in->numDim)]=in->Coordinates[INDEX2(j,i,in->numDim)];
     }
     out->isPrepared=FINLEY_UNPREPARED;
   }
}
