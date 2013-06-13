
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*   Finley: Mesh */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Finley_Mesh_relableElementNodes(index_t* newNode, index_t offset, Finley_Mesh* in)
{
    in->Elements->relabelNodes(newNode, offset);
    in->FaceElements->relabelNodes(newNode, offset);
    in->ContactElements->relabelNodes(newNode, offset);
    in->Points->relabelNodes(newNode, offset);
}
