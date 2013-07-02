
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

void IndexList_insertElements(IndexList* index_list, ElementFile* elements,
                              bool reduce_row_order, const int* row_map,
                              bool reduce_col_order, const int* col_map)
{
    // index_list is an array of linked lists. Each entry is a row (DOF) and
    // contains the indices to the non-zero columns
    if (!elements)
        return;

    const int NN=elements->numNodes;
    ReferenceElement* refElement =
        ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, false);

    int NN_row, NN_col, *row_node=NULL, *col_node=NULL, numSub;
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
                            IndexList_insertIndex(index_list[irow], icol);
                        }
                    }
                }
            }
        }
    }
}

void IndexList_insertElementsWithRowRangeNoMainDiagonal(IndexList* index_list,
                            int firstRow, int lastRow, ElementFile* elements,
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
                                IndexList_insertIndex(index_list[irow_loc], icol);
                        }
                    }
                }
            }
        }
    }
}

/// inserts a row index into the IndexList in if it does not exist
void IndexList_insertIndex(IndexList& in, int index)
{
    if (std::find(in.begin(), in.end(), index) == in.end())
        in.push_back(index);
}

/// counts the number of row indices in the IndexList in
int IndexList_count(const IndexList& in, int range_min, int range_max)
{
    int out=0;
    for (IndexList::const_iterator it=in.begin(); it!=in.end(); it++) {
        if (*it >= range_min && range_max > *it)
            ++out;
    }
    return out;
}

/// count the number of row indices in the IndexList in
void IndexList_toArray(const IndexList& in, int* array, int range_min, int range_max, int index_offset)
{
    int idx=0;
    for (IndexList::const_iterator it=in.begin(); it!=in.end(); it++) {
        if (*it >= range_min && range_max > *it) {
            array[idx]=(*it)+index_offset;
            ++idx;
        }
    }
}

/// creates a Paso_pattern from a range of indices
Paso_Pattern* IndexList_createPattern(int n0, int n,
                                      const IndexList* index_list,
                                      int range_min, int range_max,
                                      int index_offset)
{
    int *ptr=new int[n+1-n0];
    // get the number of connections per row
#pragma omp parallel for
    for (int i=n0; i<n; ++i) {
        ptr[i-n0]=IndexList_count(index_list[i], range_min, range_max);
    }
    // accumulate ptr
    int s=0;
    for (int i=n0; i<n; ++i) {
        const int itmp=ptr[i-n0];
        ptr[i-n0]=s;
        s+=itmp;
    }
    ptr[n-n0]=s;
    // fill index
    int *index=new int[ptr[n-n0]];
#pragma omp parallel for
    for (int i=n0; i<n; ++i) {
        IndexList_toArray(index_list[i], &index[ptr[i-n0]], range_min, range_max, index_offset);
    }

    Paso_Pattern* out=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, n-n0, range_max+index_offset, ptr, index);
    if (!Finley_noError()) {
        delete[] ptr;
        delete[] index;
        Paso_Pattern_free(out);
    }
    return out;
}

} // namespace finley

