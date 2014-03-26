
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


/****************************************************************************

  Finley: Converting an element list into a matrix shape

*****************************************************************************/

#ifndef __FINLEY_INDEXLIST_H__
#define __FINLEY_INDEXLIST_H__

#include "Finley.h"

#include <list>

namespace paso {
    struct Pattern;
}

// helpers to build system matrix

namespace finley {

class ElementFile;

typedef std::list<int> IndexList;

void IndexList_insertIndex(IndexList&, int);
void IndexList_toArray(const IndexList&, int*, int, int, int);
int IndexList_count(const IndexList&,  int, int);

paso::Pattern* IndexList_createPattern(int n0, int n,
        const IndexList* index_list, int range_min, int range_max,
        int index_offset);

void IndexList_insertElements(IndexList* index_list, ElementFile* elements,
                              bool reduce_row_order, const int* row_map,
                              bool reduce_col_order, const int* col_map);

void IndexList_insertElementsWithRowRangeNoMainDiagonal(IndexList* index_list,
        int firstRow, int lastRow, ElementFile* elements, int* row_map,
        int* col_map);

} // namespace finley

#endif // __FINLEY_INDEXLIST_H__

