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

/*   Finley: ElementFile */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

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
