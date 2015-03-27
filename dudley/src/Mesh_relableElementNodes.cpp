
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

/*   Dudley: Mesh */

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Dudley_Mesh_relableElementNodes(index_t * newNode, index_t offset, Dudley_Mesh * in)
{
    Dudley_ElementFile_relableNodes(newNode, offset, in->Elements);
    Dudley_ElementFile_relableNodes(newNode, offset, in->FaceElements);
    Dudley_ElementFile_relableNodes(newNode, offset, in->Points);
}
