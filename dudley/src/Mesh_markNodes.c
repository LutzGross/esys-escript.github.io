
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

/*   mark the used nodes with offeset: */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_markNodes(index_t * mask, index_t offset, Dudley_Mesh * in, bool_t useLinear)
{
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->Elements, useLinear);
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->FaceElements, useLinear);
    Dudley_ElementFile_markNodes(mask, offset, in->Nodes->numNodes, in->Points, useLinear);
}

void Dudley_Mesh_markDOFsConnectedToRange(index_t * mask, index_t offset, index_t marker,
					  index_t firstDOF, index_t lastDOF, Dudley_Mesh * in, bool_t useLinear)
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
