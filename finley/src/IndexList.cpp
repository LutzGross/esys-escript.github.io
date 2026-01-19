
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
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

#include "IndexList.h"
#include "ElementFile.h"

#include <escript/index.h>

namespace finley {

/* Translate from distributed/local array indices to global indices */

/// inserts the contributions from the element matrices of elements
/// into the row index col.
void IndexList_insertElements(IndexList* index_list, ElementFile* elements,
                              bool reduce_row_order, const index_t* row_map,
                              bool reduce_col_order, const index_t* col_map)
{
    // index_list is an array of linked lists. Each entry is a row (DOF) and
    // contains the indices to the non-zero columns
    if (!elements)
        return;

    const int NN = elements->numNodes;
    const_ReferenceElement_ptr refElement(elements->referenceElementSet->
                                            borrowReferenceElement(false));

    int NN_row, NN_col, numSub;
    const int *row_node=NULL, *col_node=NULL;
    if (reduce_col_order) {
        numSub=1;
        col_node=refElement->Type->linearNodes;
        NN_col=refElement->LinearBasisFunctions->Type->numShapes * refElement->Type->numSides;
    } else {
        numSub=refElement->Type->numSubElements;
        col_node=refElement->Type->subElementNodes;
        NN_col=refElement->BasisFunctions->Type->numShapes * refElement->Type->numSides;
    }

    if (reduce_row_order) {
        numSub=1;
        row_node=refElement->Type->linearNodes;
        NN_row=refElement->LinearBasisFunctions->Type->numShapes * refElement->Type->numSides;
    } else {
        numSub=refElement->Type->numSubElements;
        row_node=refElement->Type->subElementNodes;
        NN_row=refElement->BasisFunctions->Type->numShapes * refElement->Type->numSides;
    }

    for (int color=elements->minColor; color<=elements->maxColor; color++) {
#pragma omp for
        for (index_t e=0; e<elements->numElements; e++) {
            if (elements->Color[e]==color) {
                for (int isub=0; isub<numSub; isub++) {
                    for (int kr=0; kr<NN_row; kr++) {
                        const index_t irow=row_map[elements->Nodes[INDEX2(row_node[INDEX2(kr,isub,NN_row)],e,NN)]];
                        for (int kc=0; kc<NN_col; kc++) {
                            const index_t icol=col_map[elements->Nodes[INDEX2(col_node[INDEX2(kc,isub,NN_col)],e,NN)]];
                            index_list[irow].insertIndex(icol);
                        }
                    }
                }
            }
        }
    }
}

void IndexList_insertElementsWithRowRangeNoMainDiagonal(
                            IndexList* index_list, index_t firstRow,
                            index_t lastRow, ElementFile* elements,
                            index_t* row_map, index_t* col_map)
{
    if (!elements)
        return;

    // this does not resolve macro elements
    const int NN = elements->numNodes;
    for (index_t color = elements->minColor; color <= elements->maxColor; color++) {
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            if (elements->Color[e] == color) {
                for (int kr = 0; kr < NN; kr++) {
                    const index_t irow = row_map[elements->Nodes[INDEX2(kr, e, NN)]];
                    if (firstRow <= irow && irow < lastRow) {
                        const index_t irow_loc = irow - firstRow;
                        for (int kc = 0; kc < NN; kc++) {
                            const index_t icol = col_map[elements->Nodes[INDEX2(kc, e, NN)]];
                            if (icol != irow)
                                index_list[irow_loc].insertIndex(icol);
                        }
                    }
                }
            }
        }
    }
}


} // namespace finley

