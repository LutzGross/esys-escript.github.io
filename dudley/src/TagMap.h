
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

/* Dudley: simple link list to privide clear names for a tag keys */

/****************************************************************************/

#ifndef __DUDLEY_TAGMAP_H__
#define __DUDLEY_TAGMAP_H__

#include "Dudley.h"

namespace dudley {

struct Dudley_TagMap {
    char *name;
    index_t tag_key;
    struct Dudley_TagMap *next;
};

void Dudley_TagMap_insert(Dudley_TagMap **, const char *name, index_t tag_key);
index_t Dudley_TagMap_getTag(Dudley_TagMap *, const char *name);
bool Dudley_TagMap_isValidTagName(Dudley_TagMap *, const char *name);
void Dudley_TagMap_free(Dudley_TagMap *);

} // namespace dudley

#endif // __DUDLEY_TAGMAP_H__

