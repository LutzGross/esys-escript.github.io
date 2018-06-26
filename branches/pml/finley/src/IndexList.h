
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
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

