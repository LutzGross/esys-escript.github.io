
/* $Id$ */

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

/*   Finley: Mesh tagmaps: provides access to the mesh tagmap */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

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
