
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: Mesh: NodeFile */

/*   scatters the NodeFile in into NodeFile out using index[0:in->numNodes-1].  */
/*   index has to be between 0 and out->numNodes-1. */
/*   coloring is choosen for the worst case */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_scatterEntries(dim_t n, index_t* index, index_t min_index, index_t max_index,
                                   index_t* Id_out, index_t* Id_in,
                                   index_t* Tag_out, index_t* Tag_in,
                                   index_t* globalDegreesOfFreedom_out, index_t* globalDegreesOfFreedom_in,
                                   dim_t numDim, double* Coordinates_out, double* Coordinates_in)
{
   dim_t i;
   register index_t k;
   const index_t range=max_index-min_index;
   const size_t numDim_size=(size_t)numDim*sizeof(double);

   #pragma omp parallel for private(i,k) schedule(static)
   for (i=0;i<n;i++) {
      k=index[i]-min_index;
      if ((k>=0) && (k <range)) {
         Id_out[k]=Id_in[i];
         Tag_out[k]=Tag_in[i];
         globalDegreesOfFreedom_out[k]=globalDegreesOfFreedom_in[i];
         memcpy(&(Coordinates_out[INDEX2(0,k,numDim)]), &(Coordinates_in[INDEX2(0,i,numDim)]), numDim_size);
      }
   }
}

void Finley_NodeFile_scatter(index_t* index, Finley_NodeFile* in, Finley_NodeFile* out)
{
   Finley_NodeFile_scatterEntries(out->numNodes, index, 0, in->numNodes,
                                  out->Id, in->Id,
                                  out->Tag, in->Tag,
                                  out->globalDegreesOfFreedom, in->globalDegreesOfFreedom,
                                  out->numDim, out->Coordinates, in->Coordinates);
}
