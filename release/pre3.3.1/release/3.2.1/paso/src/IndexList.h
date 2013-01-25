
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Paso: pattern                                            */

/**************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_INDEX_LIST
#define INC_PASO_INDEX_LIST

#include "Common.h"

/**************************************************************/

#define INDEXLIST_LENGTH 85

typedef struct Paso_IndexList {
  index_t index[INDEXLIST_LENGTH];
  dim_t n;
  struct Paso_IndexList *extension;
} Paso_IndexList;

typedef struct Paso_IndexListArray {
   dim_t n;
   Paso_IndexList* index_list;
} Paso_IndexListArray;

PASO_DLL_API
void Paso_IndexList_insertIndex(Paso_IndexList*, index_t);

PASO_DLL_API
void Paso_IndexList_toArray(Paso_IndexList*, index_t*, index_t, index_t, index_t);

PASO_DLL_API
dim_t Paso_IndexList_count(Paso_IndexList*,  index_t, index_t);

PASO_DLL_API
void Paso_IndexList_free(Paso_IndexList*);

PASO_DLL_API
Paso_IndexListArray* Paso_IndexListArray_alloc(const dim_t n);

PASO_DLL_API
void Paso_IndexListArray_free(Paso_IndexListArray* in);

#define Paso_IndexListArray_insertIndex(__a__,i,j) Paso_IndexList_insertIndex(&((__a__)->index_list[i]),j);

#endif /* #ifndef INC_PASO_INDEX_LIST */
