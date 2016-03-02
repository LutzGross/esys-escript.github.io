
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

/****************************************************************************/

/* Dudley: simple mapping from names to tag keys via a linked list */

/****************************************************************************/

#include "TagMap.h"

namespace dudley {

void Dudley_TagMap_insert(Dudley_TagMap** tag_map, const char* name, index_t tag_key)
{
    Dudley_TagMap *map = NULL;
    if (strlen(name) < 1) {
        throw escript::ValueError("TagMap_insert: empty tag name.");
    }
    if (strchr(name, 32) != NULL) {
        throw escript::ValueError("TagMap_insert: tag name may not contain a space.");
    }
    if (*tag_map == NULL) {
        map = new Dudley_TagMap;
        map->name = new char[strlen(name) + 1];
        strcpy(map->name, name);
        map->tag_key = tag_key;
        map->next = NULL;
        *tag_map = map;
    } else {
        if (strcmp((*tag_map)->name, name) == 0) {
            (*tag_map)->tag_key = tag_key;
        } else {
            Dudley_TagMap_insert(&((*tag_map)->next), name, tag_key);
        }
    }
}

index_t Dudley_TagMap_getTag(Dudley_TagMap * tag_map, const char *name)
{
    if (tag_map == NULL) {
        std::stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        const std::string errorMsg(ss.str());
        throw escript::ValueError(errorMsg);
    } else {
        if (strcmp(tag_map->name, name) == 0) {
            return tag_map->tag_key;
        } else {
            return Dudley_TagMap_getTag(tag_map->next, name);
        }
    }
}

bool Dudley_TagMap_isValidTagName(Dudley_TagMap * tag_map, const char *name)
{
    if (tag_map == NULL)
        return false;

    if (strcmp(tag_map->name, name) == 0)
        return true;

    return Dudley_TagMap_isValidTagName(tag_map->next, name);
}

/* deallocates the Dudley_TagMap in by recursive calls */

void Dudley_TagMap_free(Dudley_TagMap * in)
{
    if (in != NULL)
    {
        Dudley_TagMap_free(in->next);
        delete[] in->name;
        delete in;
    }
}

} // namespace dudley

