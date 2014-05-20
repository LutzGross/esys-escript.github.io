
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

/*   Dudley: Mesh tagmaps: provides access to the mesh tagmap */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Dudley_Mesh_addTagMap(Dudley_Mesh * mesh_p, const char *name, index_t tag_key)
{
    Dudley_TagMap_insert(&(mesh_p->TagMap), name, tag_key);
}

index_t Dudley_Mesh_getTag(Dudley_Mesh * mesh_p, const char *name)
{
    return Dudley_TagMap_getTag(mesh_p->TagMap, name);
}

bool Dudley_Mesh_isValidTagName(Dudley_Mesh * mesh_p, const char *name)
{
    return Dudley_TagMap_isValidTagName(mesh_p->TagMap, name);
}
