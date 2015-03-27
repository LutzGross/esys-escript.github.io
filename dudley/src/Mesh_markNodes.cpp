
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

/*   mark the used nodes with offset: */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Dudley_Mesh_markNodes(index_t * mask, index_t offset, Dudley_Mesh * in, bool useLinear)
{
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->Elements, useLinear);
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->FaceElements, useLinear);
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->Points, useLinear);
}

void Dudley_Mesh_markDOFsConnectedToRange(index_t * mask, index_t offset, index_t marker,
					  index_t firstDOF, index_t lastDOF, Dudley_Mesh * in, bool useLinear)
{
    index_t *dofIndex;
    if (useLinear)
    {
	dofIndex = in->Nodes->globalReducedDOFIndex;
    }
    else
    {
	dofIndex = in->Nodes->globalDegreesOfFreedom;
    }
    Dudley_ElementFile_markDOFsConnectedToRange(mask, offset, marker, firstDOF, lastDOF, dofIndex, in->Elements,
						useLinear);
    Dudley_ElementFile_markDOFsConnectedToRange(mask, offset, marker, firstDOF, lastDOF, dofIndex, in->FaceElements,
						useLinear);
    Dudley_ElementFile_markDOFsConnectedToRange(mask, offset, marker, firstDOF, lastDOF, dofIndex, in->Points,
						useLinear);
}
