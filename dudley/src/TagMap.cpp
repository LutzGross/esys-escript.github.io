
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

/************************************************************************************/

/* Dudley: simple mapping from names to tag keys via a linked list */

/************************************************************************************/

#include "TagMap.h"
#include "esysUtils/mem.h"
#include "string.h"

/************************************************************************************/

void Dudley_TagMap_insert(Dudley_TagMap ** tag_map, const char *name, index_t tag_key)
{
    Dudley_TagMap *map = NULL;
    if (strlen(name) < 1)
    {
	Dudley_setError(VALUE_ERROR, "empty tag name.");
	return;
    }
    if (strchr(name, 32) != NULL)
    {				/* check for space */
	Dudley_setError(VALUE_ERROR, "tag name may not contain a space.");
	return;
    }
    if (*tag_map == NULL)
    {
	map = new Dudley_TagMap;
	if (Dudley_checkPtr(map))
	    return;
	map->name = new char[strlen(name) + 1];
	if (Dudley_checkPtr(map->name))
	{
	    delete map;
	}
	else
	{
	    strcpy(map->name, name);
	    map->tag_key = tag_key;
	    map->next = NULL;
	    *tag_map = map;
	}
    }
    else
    {
	if (strcmp((*tag_map)->name, name) == 0)
	{
	    (*tag_map)->tag_key = tag_key;
	}
	else
	{
	    Dudley_TagMap_insert(&((*tag_map)->next), name, tag_key);
	}
    }
}

index_t Dudley_TagMap_getTag(Dudley_TagMap * tag_map, const char *name)
{
    char error_msg[LenErrorMsg_MAX];
    if (tag_map == NULL)
    {
	sprintf(error_msg, "Dudley_TagMap_getTag: unknown tag name %s.", name);
	Dudley_setError(VALUE_ERROR, error_msg);
	return -1;
    }
    else
    {
	if (strcmp(tag_map->name, name) == 0)
	{
	    return tag_map->tag_key;
	}
	else
	{
	    return Dudley_TagMap_getTag(tag_map->next, name);
	}
    }
}

bool Dudley_TagMap_isValidTagName(Dudley_TagMap * tag_map, const char *name)
{
    if (tag_map == NULL)
    {
	return false;
    }
    else
    {
	if (strcmp(tag_map->name, name) == 0)
	{
	    return true;
	}
	else
	{
	    return Dudley_TagMap_isValidTagName(tag_map->next, name);
	}
    }
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
    return;
}
