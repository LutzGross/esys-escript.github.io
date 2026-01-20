
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Finley: Converting an element list into a matrix shape

*****************************************************************************/

#ifndef __FINLEY_INDEXLIST_H__
#define __FINLEY_INDEXLIST_H__

#include "Finley.h"

#include <escript/IndexList.h>

// helpers to build system matrix

namespace finley {

using escript::IndexList;

class ElementFile;

void IndexList_insertElements(IndexList* index_list, ElementFile* elements,
                              bool reduce_row_order, const index_t* row_map,
                              bool reduce_col_order, const index_t* col_map);

void IndexList_insertElementsWithRowRangeNoMainDiagonal(
        IndexList* index_list, index_t firstRow, index_t lastRow,
        ElementFile* elements, index_t* row_map, index_t* col_map);

} // namespace finley

#endif // __FINLEY_INDEXLIST_H__

