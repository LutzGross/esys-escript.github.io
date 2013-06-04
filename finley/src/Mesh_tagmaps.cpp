
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

/*   Finley: Mesh tagmaps: provides access to the mesh tag map */

/************************************************************************************/

#include "Mesh.h"

#include <sstream>

void Finley_Mesh_addTagMap(Finley_Mesh *mesh_p, const char* name, int tag_key) 
{
   mesh_p->tagMap[std::string(name)]=tag_key;
}

int Finley_Mesh_getTag(Finley_Mesh *mesh_p, const char* name)
{
    TagMap::iterator it = mesh_p->tagMap.find(name);
    if (it == mesh_p->tagMap.end()) {
        std::stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        const std::string errorMsg(ss.str());
        Finley_setError(VALUE_ERROR, errorMsg.c_str());
        return -1;
    }
    return it->second;
}

bool Finley_Mesh_isValidTagName(Finley_Mesh *mesh_p, const char* name)
{
   return (mesh_p->tagMap.count(std::string(name)) > 0);
}

