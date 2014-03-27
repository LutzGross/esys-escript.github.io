
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

#include "IndexList.h"
#include "ElementFile.h"

#include "paso/SystemMatrixPattern.h"

/* Translate from distributed/local array indices to global indices */

/****************************************************************************/
/* inserts the contributions from the element matrices of elements
   into the row index col. If symmetric is set, only the upper
   triangle of the matrix is stored.
*/

namespace finley {

void IndexList_insertElements(IndexListArray& index_list, ElementFile* elements,
                              bool reduce_row_order, const int* row_map,
                              bool reduce_col_order, const int* col_map)
{
    // index_list is an array of linked lists. Each entry is a row (DOF) and
    // contains the indices to the non-zero columns
    if (!elements)
        return;

    const int NN=elements->numNodes;
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
        for (int e=0; e<elements->numElements; e++) {
            if (elements->Color[e]==color) {
                for (int isub=0; isub<numSub; isub++) {
                    for (int kr=0; kr<NN_row; kr++) {
                        const int irow=row_map[elements->Nodes[INDEX2(row_node[INDEX2(kr,isub,NN_row)],e,NN)]];
                        for (int kc=0; kc<NN_col; kc++) {
                            const int icol=col_map[elements->Nodes[INDEX2(col_node[INDEX2(kc,isub,NN_col)],e,NN)]];
                            index_list[irow].insertIndex(icol);
                        }
                    }
                }
            }
        }
    }
}

void IndexList_insertElementsWithRowRangeNoMainDiagonal(
                            IndexListArray& index_list, int firstRow,
                            int lastRow, ElementFile* elements,
                            int* row_map, int* col_map)
{
    if (!elements)
        return;

    // this does not resolve macro elements
    const int NN=elements->numNodes;
    for (int color=elements->minColor; color<=elements->maxColor; color++) {
#pragma omp for
        for (int e=0; e<elements->numElements; e++) {
            if (elements->Color[e]==color) {
                for (int kr=0; kr<NN; kr++) {
                    const int irow=row_map[elements->Nodes[INDEX2(kr,e,NN)]];
                    if (firstRow<=irow && irow<lastRow) {
                        const int irow_loc=irow-firstRow;
                        for (int kc=0; kc<NN; kc++) {
                            const int icol=col_map[elements->Nodes[INDEX2(kc,e,NN)]];
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

