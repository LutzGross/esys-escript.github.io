/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/* Finley: Converting an element list into a matrix shape     */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_INDEXLIST
#define INC_FINLEY_INDEXLIST

#include "Finley.h"
#include "ElementFile.h"
#include "Mesh.h"

/* structure to build system matrix */

#define INDEXLIST_LENGTH 85

typedef struct Finley_IndexList {
  index_t index[INDEXLIST_LENGTH];
  dim_t n;
  struct Finley_IndexList *extension;
} Finley_IndexList;
void Finley_IndexList_insertElements(Finley_IndexList* index_list, Finley_ElementFile* elements,
                                       bool_t reduce_row_order, index_t* row_map,
                                       bool_t reduce_col_order, index_t* col_map);
void Finley_IndexList_insertIndex(Finley_IndexList*, index_t);
void Finley_IndexList_toArray(Finley_IndexList*, index_t*, index_t, index_t);
dim_t Finley_IndexList_count(Finley_IndexList*,  index_t, index_t);
void Finley_IndexList_free(Finley_IndexList*);
Paso_Pattern* Finley_IndexList_createPattern(dim_t n,Finley_IndexList* index_list,index_t range_min,index_t range_max);
Paso_Coupler* Finley_IndexList_createCoupler(Finley_IndexList *index_list);
void Finley_IndexList_insertElementsWithRowRange(Finley_IndexList* index_list, index_t firstRow, index_t lastRow,
                                                 Finley_ElementFile* elements, index_t* row_map, index_t* col_map);



#endif /* #ifndef INC_FINLEY_INDEXLIST */
