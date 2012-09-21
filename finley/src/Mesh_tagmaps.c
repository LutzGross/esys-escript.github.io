
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/*   Finley: Mesh tagmaps: provides access to the mesh tag map */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Finley_Mesh_addTagMap(Finley_Mesh *mesh_p,const char* name, index_t tag_key) 
{
   Finley_TagMap_insert(&(mesh_p->TagMap),name,tag_key);
}
index_t Finley_Mesh_getTag(Finley_Mesh *mesh_p,const char* name) {
   return Finley_TagMap_getTag(mesh_p->TagMap,name);
}
bool_t Finley_Mesh_isValidTagName(Finley_Mesh *mesh_p,const char* name) {
   return Finley_TagMap_isValidTagName(mesh_p->TagMap,name);
}
