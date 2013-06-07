
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

/* Finley: Converting an element list into a matrix shape     */

/************************************************************************************/

#ifndef INC_FINLEY_INDEXLIST
#define INC_FINLEY_INDEXLIST

#include "Finley.h"
#include "ElementFile.h"
#include "Mesh.h"

/* structure to build system matrix */

#define INDEXLIST_LENGTH 85

typedef struct Finley_IndexList {
  int index[INDEXLIST_LENGTH];
  int n;
  struct Finley_IndexList *extension;
} Finley_IndexList;

void Finley_IndexList_insertElements(Finley_IndexList* index_list, ElementFile* elements,
                                     bool_t reduce_row_order, int* row_map,
                                     bool_t reduce_col_order, int* col_map);
void Finley_IndexList_insertIndex(Finley_IndexList*, int);
void Finley_IndexList_toArray(Finley_IndexList*, int*, int, int, int);
dim_t Finley_IndexList_count(Finley_IndexList*,  int, int);
void Finley_IndexList_free(Finley_IndexList*);
Paso_Pattern* Finley_IndexList_createPattern(int n0, int n, Finley_IndexList* index_list, int range_min, int range_max, int index_offset);
void Finley_IndexList_insertElementsWithRowRange(Finley_IndexList* index_list, int firstRow, int lastRow,
                                                 ElementFile* elements, int* row_map, int* col_map);
void Finley_IndexList_insertElementsWithRowRangeNoMainDiagonal(Finley_IndexList* index_list, int firstRow, int lastRow,
                                                               ElementFile* elements, int* row_map, int* col_map);


#endif /* #ifndef INC_FINLEY_INDEXLIST */

