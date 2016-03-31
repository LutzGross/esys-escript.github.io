
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************/
/*   Dudley: Mesh */
/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */
/****************************************************************************/

#include "Mesh.h"

namespace dudley {

void Dudley_Mesh_relableElementNodes(index_t* newNode, index_t offset, Dudley_Mesh* in)
{
    Dudley_ElementFile_relableNodes(newNode, offset, in->Elements);
    Dudley_ElementFile_relableNodes(newNode, offset, in->FaceElements);
    Dudley_ElementFile_relableNodes(newNode, offset, in->Points);
}

} // namespace dudley

