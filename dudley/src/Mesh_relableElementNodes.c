
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Dudley: Mesh */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_relableElementNodes(index_t* newNode,index_t offset,Dudley_Mesh* in) {
      Dudley_ElementFile_relableNodes(newNode,offset,in->Elements);
      Dudley_ElementFile_relableNodes(newNode,offset,in->FaceElements);
      Dudley_ElementFile_relableNodes(newNode,offset,in->Points);
}
