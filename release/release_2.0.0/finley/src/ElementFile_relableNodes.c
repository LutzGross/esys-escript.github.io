
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

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_relableNodes(index_t* newNode,index_t offset,Finley_ElementFile* in) {
   dim_t i,j,NN;

   if (in!=NULL) {
     NN=in->ReferenceElement->Type->numNodes;
     #pragma omp parallel for private(j,i) schedule(static)
     for (j=0;j<in->numElements;j++) {
       for (i=0;i<NN;i++) {
         in->Nodes[INDEX2(i,j,NN)]=newNode[in->Nodes[INDEX2(i,j,NN)]-offset];
       }
     }
   }
}
