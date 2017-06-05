
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

#include "IndexList.h"
#include "ElementFile.h"

#include <escript/index.h>

namespace dudley {

/* Translate from distributed/local array indices to global indices */

/// inserts the contributions from the element matrices of elements
/// into the row index col.
void IndexList_insertElements(IndexList* index_list,
                              const ElementFile* elements, const index_t* map)
{
    // index_list is an array of linked lists. Each entry is a row (DOF) and
    // contains the indices to the non-zero columns
    if (!elements)
        return;

    const int NN = elements->numNodes;
    // number of element nodes for both column and row
    const int NN_rowcol = elements->numShapes;

    for (index_t color = elements->minColor; color <= elements->maxColor; color++) {
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            if (elements->Color[e] == color) {
                for (int kr = 0; kr < NN_rowcol; kr++) {
                    const index_t irow = map[elements->Nodes[INDEX2(kr, e, NN)]];
                    for (int kc = 0; kc < NN_rowcol; kc++) {
                        const index_t icol = map[elements->Nodes[INDEX2(kc, e, NN)]];
                        index_list[irow].insertIndex(icol);
                    }
                }
            }
        }
    }
}

void IndexList_insertElementsWithRowRangeNoMainDiagonal(IndexList* indexList,
                              index_t firstRow, index_t lastRow,
                              const ElementFile* elements, const index_t* map)
{
    if (!elements)
        return;

    const int NN = elements->numNodes;
    for (index_t color = elements->minColor; color <= elements->maxColor; color++) {
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            if (elements->Color[e] == color) {
                for (int kr = 0; kr < NN; kr++) {
                    const index_t irow = map[elements->Nodes[INDEX2(kr, e, NN)]];
                    if (firstRow <= irow && irow < lastRow) {
                        const index_t irow_loc = irow - firstRow;
                        for (int kc = 0; kc < NN; kc++) {
                            const index_t icol = map[elements->Nodes[INDEX2(kc, e, NN)]];
                            if (icol != irow)
                                indexList[irow_loc].insertIndex(icol);
                        }
                    }
                }
            }
        }
    }
}

} // namespace dudley

