
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

/****************************************************************************/

/* Dudley: Converting an element list into a matrix shape     */

/****************************************************************************/

#ifndef __DUDLEY_INDEXLIST_H__
#define __DUDLEY_INDEXLIST_H__

#include "Dudley.h"

#include <escript/IndexList.h>

namespace dudley {

using escript::IndexList;

// helpers to build system matrix

class ElementFile;

void IndexList_insertElements(IndexList* indexlist, const ElementFile* elements,
                              const index_t* map);

void IndexList_insertElementsWithRowRangeNoMainDiagonal(IndexList* index_list,
                            index_t firstRow, index_t lastRow,
                            const ElementFile* elements, const index_t* map);

} // namespace dudley

#endif // __DUDLEY_INDEXLIST_H__

