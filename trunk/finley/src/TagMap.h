
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

/* Finley: simple linked list to provide clear names for tag keys */

/************************************************************************************/

#ifndef INC_FINLEY_TAGMAP
#define INC_FINLEY_TAGMAP

#include "Finley.h"


typedef struct Finley_TagMap {
  char* name;
  index_t tag_key;
  struct Finley_TagMap *next;
} Finley_TagMap;

void Finley_TagMap_insert(Finley_TagMap**,const char* name, index_t tag_key);
index_t Finley_TagMap_getTag(Finley_TagMap*,const char* name);
bool_t Finley_TagMap_isValidTagName(Finley_TagMap*,const char* name);
void Finley_TagMap_free(Finley_TagMap*);
#endif /* #ifndef INC_FINLEY_TAGMAP */

