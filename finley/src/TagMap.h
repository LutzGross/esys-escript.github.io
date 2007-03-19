/*
 ************************************************************
 *          Copyright 2007 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/* Finley: simple link list to privide clear names for a tag keys */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id:*/

/**************************************************************/

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

