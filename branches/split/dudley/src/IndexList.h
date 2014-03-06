
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

/************************************************************************************/

/* Dudley: Converting an element list into a matrix shape     */

/************************************************************************************/

#ifndef INC_DUDLEY_INDEXLIST
#define INC_DUDLEY_INDEXLIST

#include "Dudley.h"
#include "ElementFile.h"
#include "Mesh.h"

/* structure to build system matrix */

#define INDEXLIST_LENGTH 85

typedef struct Dudley_IndexList {
    index_t index[INDEXLIST_LENGTH];
    dim_t n;
    struct Dudley_IndexList *extension;
} Dudley_IndexList;
void Dudley_IndexList_insertElements(Dudley_IndexList * index_list, Dudley_ElementFile * elements,
				     bool reduce_row_order, index_t * row_map,
				     bool reduce_col_order, index_t * col_map);
void Dudley_IndexList_insertIndex(Dudley_IndexList *, index_t);
void Dudley_IndexList_toArray(Dudley_IndexList *, index_t *, index_t, index_t, index_t);
dim_t Dudley_IndexList_count(Dudley_IndexList *, index_t, index_t);
void Dudley_IndexList_free(Dudley_IndexList *);
Paso_Pattern *Dudley_IndexList_createPattern(dim_t n0, dim_t n, Dudley_IndexList * index_list, index_t range_min,
					     index_t range_max, index_t index_offset);
void Dudley_IndexList_insertElementsWithRowRange(Dudley_IndexList * index_list, index_t firstRow, index_t lastRow,
						 Dudley_ElementFile * elements, index_t * row_map, index_t * col_map);
void Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(Dudley_IndexList * index_list, index_t firstRow,
							       index_t lastRow, Dudley_ElementFile * elements,
							       index_t * row_map, index_t * col_map);

#endif				/* #ifndef INC_DUDLEY_INDEXLIST */
