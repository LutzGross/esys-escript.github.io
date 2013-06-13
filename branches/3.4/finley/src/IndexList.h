
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


/****************************************************************************

  Finley: Converting an element list into a matrix shape

*****************************************************************************/

#ifndef __FINLEY_INDEXLIST_H__
#define __FINLEY_INDEXLIST_H__

#include "Finley.h"

#include <list>

struct Paso_Pattern;

// helpers to build system matrix

namespace finley {

class ElementFile;

typedef std::list<int> IndexList;

void IndexList_insertIndex(IndexList&, int);
void IndexList_toArray(const IndexList&, int*, int, int, int);
int IndexList_count(const IndexList&,  int, int);

Paso_Pattern* IndexList_createPattern(int n0, int n,
        const IndexList* index_list, int range_min, int range_max,
        int index_offset);

void IndexList_insertElements(IndexList* index_list, ElementFile* elements,
                              bool reduce_row_order, int* row_map,
                              bool reduce_col_order, int* col_map);

void IndexList_insertElementsWithRowRangeNoMainDiagonal(IndexList* index_list,
        int firstRow, int lastRow, ElementFile* elements, int* row_map,
        int* col_map);

} // namespace finley

#endif // __FINLEY_INDEXLIST_H__

