
/* $Id: Mesh_relableElementNodes.c 1306 2007-09-18 05:51:09Z ksteube $ */

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

/*   Finley: Mesh */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_relableElementNodes(index_t* newNode,index_t offset,Finley_Mesh* in) {
      Finley_ElementFile_relableNodes(newNode,offset,in->Elements);
      Finley_ElementFile_relableNodes(newNode,offset,in->FaceElements);
      Finley_ElementFile_relableNodes(newNode,offset,in->ContactElements);
      Finley_ElementFile_relableNodes(newNode,offset,in->Points);
}
