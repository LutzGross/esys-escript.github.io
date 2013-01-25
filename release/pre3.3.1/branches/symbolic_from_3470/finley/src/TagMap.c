
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Finley: simple mapping from names to tag keys via a linked list */

/**************************************************************/

#include "TagMap.h"
#include "esysUtils/mem.h"
#include <string.h>

/**************************************************************/

void Finley_TagMap_insert(Finley_TagMap** tag_map, 
                          const char* name,
                          index_t tag_key) 
{
  Finley_TagMap* map=NULL;
  if (strlen(name)<1) {
     Finley_setError(VALUE_ERROR,"empty tag name.");
     return;
  }
  if (strchr(name,32) != NULL) {   /* check for space */
     Finley_setError(VALUE_ERROR,"tag name may not contain a space.");
     return;
  }
  if (*tag_map == NULL) {
       map=MEMALLOC(1,Finley_TagMap);
       if (Finley_checkPtr(map)) return;
       map->name=MEMALLOC(strlen(name)+1,char);
       if (Finley_checkPtr(map->name) ) {
            MEMFREE(map);
        } else {
            strcpy(map->name,name);
            map->tag_key=tag_key;
            map->next=NULL;
            *tag_map=map;
        } 
  } else {
      if (strcmp((*tag_map)->name,name)==0) {
         (*tag_map)->tag_key=tag_key;
      } else {
         Finley_TagMap_insert(&((*tag_map)->next),name,tag_key);
      }
  }
}

index_t Finley_TagMap_getTag(Finley_TagMap* tag_map,const char* name)
{
  char error_msg[LenErrorMsg_MAX]; 
  if (tag_map == NULL) {
        sprintf(error_msg,"Finley_TagMap_getTag: unknown tag name %s.",name);
        Finley_setError(VALUE_ERROR,error_msg);
        return -1;
  } else {
      if (strcmp(tag_map->name,name)==0) {
         return tag_map->tag_key;
      } else {
         return Finley_TagMap_getTag(tag_map->next,name);
      }
  }
}
bool_t Finley_TagMap_isValidTagName(Finley_TagMap* tag_map, const char* name)
{
  if (tag_map == NULL) {
        return FALSE;
  } else {
      if (strcmp(tag_map->name,name)==0) {
         return TRUE;
      } else {
         return Finley_TagMap_isValidTagName(tag_map->next,name);
      }
  }
}
/* deallocates the Finley_TagMap in by recursive calls */

void Finley_TagMap_free(Finley_TagMap* in) {
  if (in!=NULL) {
    Finley_TagMap_free(in->next);
    MEMFREE(in->name);
    MEMFREE(in);
  }
  return;
}

