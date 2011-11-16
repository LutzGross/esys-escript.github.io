
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

/* Dudley: simple link list to privide clear names for a tag keys */

/**************************************************************/

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
