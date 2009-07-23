
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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

/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_PATTERN
#define INC_PASO_PATTERN

#include "Common.h"

/**************************************************************/

#define PATTERN_FORMAT_DEFAULT 0
#define PATTERN_FORMAT_SYM 1
#define PATTERN_FORMAT_OFFSET1 2

typedef struct Paso_Pattern {
  int type;
  dim_t numOutput;
  dim_t numInput;
  dim_t len;
  index_t* ptr;
  index_t* index;
  dim_t reference_counter;
} Paso_Pattern;

#define INDEXLIST_LENGTH 85

typedef struct Paso_IndexList {
  index_t index[INDEXLIST_LENGTH];
  dim_t n;
  struct Paso_IndexList *extension;
} Paso_IndexList;
/*  interfaces: */


PASO_DLL_API
Paso_Pattern* Paso_Pattern_alloc(int type, dim_t numOutput, dim_t numInput, index_t* ptr, index_t* index);

PASO_DLL_API

PASO_DLL_API
Paso_Pattern* Paso_Pattern_getReference(Paso_Pattern*);

PASO_DLL_API
void Paso_Pattern_free(Paso_Pattern*);

PASO_DLL_API
int Paso_comparIndex(const void *,const void *);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_unrollBlocks(Paso_Pattern*,int, dim_t,dim_t);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_getSubpattern(Paso_Pattern*,dim_t,dim_t,index_t*,index_t*);

PASO_DLL_API
bool_t Paso_Pattern_isEmpty(Paso_Pattern* in);

PASO_DLL_API
void Paso_Pattern_mis(Paso_Pattern* pattern_p, index_t* mis_marker);

PASO_DLL_API
void Paso_Pattern_reduceBandwidth(Paso_Pattern* self,index_t* oldToNew);

PASO_DLL_API
void Paso_Pattern_color(Paso_Pattern* patter, index_t* num_colors, index_t* colorOf);
Paso_Pattern* Paso_Pattern_multiply(int type, Paso_Pattern* A, Paso_Pattern* B);

PASO_DLL_API
Paso_Pattern* Paso_Pattern_binop(int type, Paso_Pattern* A, Paso_Pattern* B);


PASO_DLL_API
void Paso_IndexList_insertIndex(Paso_IndexList*, index_t);

PASO_DLL_API
void Paso_IndexList_toArray(Paso_IndexList*, index_t*, index_t, index_t, index_t);

PASO_DLL_API
dim_t Paso_IndexList_count(Paso_IndexList*,  index_t, index_t);

PASO_DLL_API
void Paso_IndexList_free(Paso_IndexList*);

PASO_DLL_API
Paso_Pattern* Paso_IndexList_createPattern(dim_t n0, dim_t n,Paso_IndexList* index_list,index_t range_min,index_t range_max, index_t index_offset);

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */
