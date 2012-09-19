
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/**************************************************************/

/*   Dudley: Mesh tagmaps: provides access to the mesh tagmap */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_addTagMap(Dudley_Mesh * mesh_p, const char *name, index_t tag_key)
{
    Dudley_TagMap_insert(&(mesh_p->TagMap), name, tag_key);
}

index_t Dudley_Mesh_getTag(Dudley_Mesh * mesh_p, const char *name)
{
    return Dudley_TagMap_getTag(mesh_p->TagMap, name);
}

bool_t Dudley_Mesh_isValidTagName(Dudley_Mesh * mesh_p, const char *name)
{
    return Dudley_TagMap_isValidTagName(mesh_p->TagMap, name);
}
