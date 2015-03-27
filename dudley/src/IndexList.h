
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

#include "esysUtils/IndexList.h"

using esysUtils::IndexList;

/* structure to build system matrix */

void Dudley_IndexList_insertElements(IndexList* index_list,
                     Dudley_ElementFile * elements,
                     bool reduce_row_order, index_t * row_map,
                     bool reduce_col_order, index_t * col_map);
void Dudley_IndexList_insertElementsWithRowRange(IndexList* index_list,
                     index_t firstRow, index_t lastRow,
                     Dudley_ElementFile* elements, index_t* row_map,
                     index_t* col_map);
void Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(
                     IndexList* index_list, index_t firstRow,
                     index_t lastRow, Dudley_ElementFile* elements,
                     index_t* row_map, index_t* col_map);

#endif /* #ifndef INC_DUDLEY_INDEXLIST */

