
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

/* Dudley: simple link list to privide clear names for a tag keys */

/************************************************************************************/

#ifndef INC_DUDLEY_TAGMAP
#define INC_DUDLEY_TAGMAP

#include "Dudley.h"

typedef struct Dudley_TagMap {
    char *name;
    index_t tag_key;
    struct Dudley_TagMap *next;
} Dudley_TagMap;

void Dudley_TagMap_insert(Dudley_TagMap **, const char *name, index_t tag_key);
index_t Dudley_TagMap_getTag(Dudley_TagMap *, const char *name);
bool_t Dudley_TagMap_isValidTagName(Dudley_TagMap *, const char *name);
void Dudley_TagMap_free(Dudley_TagMap *);
#endif				/* #ifndef INC_DUDLEY_TAGMAP */
