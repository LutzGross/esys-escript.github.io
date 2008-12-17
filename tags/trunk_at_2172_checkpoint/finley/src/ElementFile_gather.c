
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: ElementFile */

/*   gathers the ElementFile out from the  ElementFile in using index[0:out->numElements-1].  */
/*   index has to be between 0 and in->numElements-1. */
/*   a conservative assumtion on the coloring is made */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_gather(index_t* index, Finley_ElementFile* in, Finley_ElementFile* out) {
   index_t k;
   dim_t e,j;
   dim_t NN_in=in->numNodes;
   dim_t NN_out=out->numNodes;
   if (in!=NULL) {
     /*OMP */
     #pragma omp parallel for private(e,k,j) schedule(static)
     for (e=0;e<out->numElements;e++) {
        k=index[e];
        out->Id[e]=in->Id[k];
        out->Tag[e]=in->Tag[k];
        out->Owner[e]=in->Owner[k];
        out->Color[e]=in->Color[k]+out->maxColor+1;
        for(j=0;j<MIN(NN_out,NN_in);j++) out->Nodes[INDEX2(j,e,NN_out)]=in->Nodes[INDEX2(j,k,NN_in)];
     }
     out->minColor=MIN(out->minColor,in->minColor+out->maxColor+1);
     out->maxColor=MAX(out->maxColor,in->maxColor+out->maxColor+1);
   }
}
