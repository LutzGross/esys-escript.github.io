
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

#include "IndexList.h"

/* Translate from distributed/local array indices to global indices */

/************************************************************************************/
/* inserts the contributions from the element matrices of elements
   into the row index col. If symmetric is set, only the upper
   triangle of the matrix is stored. */

void Dudley_IndexList_insertElements(IndexList* index_list,
                                     Dudley_ElementFile * elements,
                                     bool reduce_row_order, index_t* row_map,
                                     bool reduce_col_order, index_t* col_map)
{
    /* index_list is an array of linked lists. Each entry is a row (DOF) and contains the indices to the non-zero columns */
    index_t color;
    dim_t e, kr, kc, NN_row, NN_col, icol, irow, NN;
    if (elements != NULL)
    {
    NN = elements->numNodes;
    NN_col = (elements->numShapes);
    NN_row = (elements->numShapes);

    for (color = elements->minColor; color <= elements->maxColor; color++)
    {
#pragma omp for private(e,irow,kr,kc,icol) schedule(static)
        for (e = 0; e < elements->numElements; e++)
        {
        if (elements->Color[e] == color)
        {
            for (kr = 0; kr < NN_row; kr++)
            {
            irow = row_map[elements->Nodes[INDEX2(kr, e, NN)]];
            for (kc = 0; kc < NN_col; kc++)
            {
                icol = col_map[elements->Nodes[INDEX2(kc, e, NN)]];
                index_list[irow].insertIndex(icol);
            }
            }
        }
        }
    }
    }
}

void Dudley_IndexList_insertElementsWithRowRange(IndexList* index_list,
        index_t firstRow, index_t lastRow, Dudley_ElementFile* elements,
        index_t* row_map, index_t* col_map)
{
/* this does not resolve macro elements */
    index_t color;
    dim_t e, kr, kc, icol, irow, NN;
    if (elements != NULL)
    {
    NN = elements->numNodes;
    for (color = elements->minColor; color <= elements->maxColor; color++)
    {
#pragma omp for private(e,irow,kr,kc,icol) schedule(static)
        for (e = 0; e < elements->numElements; e++)
        {
        if (elements->Color[e] == color)
        {
            for (kr = 0; kr < NN; kr++)
            {
            irow = row_map[elements->Nodes[INDEX2(kr, e, NN)]];
            if ((firstRow <= irow) && (irow < lastRow))
            {
                irow -= firstRow;
                for (kc = 0; kc < NN; kc++)
                {
                icol = col_map[elements->Nodes[INDEX2(kc, e, NN)]];
                index_list[irow].insertIndex(icol);
                }
            }
            }
        }
        }
    }
    }
}

void Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(
        IndexList* index_list, index_t firstRow, index_t lastRow,
        Dudley_ElementFile* elements, index_t* row_map, index_t* col_map)
{
    /* this does not resolve macro elements */
    index_t color;
    dim_t e, kr, kc, icol, irow, NN, irow_loc;
    if (elements != NULL)
    {
    NN = elements->numNodes;
    for (color = elements->minColor; color <= elements->maxColor; color++)
    {
#pragma omp for private(e,irow,kr,kc,icol,irow_loc) schedule(static)
        for (e = 0; e < elements->numElements; e++)
        {
        if (elements->Color[e] == color)
        {
            for (kr = 0; kr < NN; kr++)
            {
            irow = row_map[elements->Nodes[INDEX2(kr, e, NN)]];
            if ((firstRow <= irow) && (irow < lastRow))
            {
                irow_loc = irow - firstRow;
                for (kc = 0; kc < NN; kc++)
                {
                icol = col_map[elements->Nodes[INDEX2(kc, e, NN)]];
                if (icol != irow)
                    index_list[irow_loc].insertIndex(icol);
                }
            }
            }
        }
        }
    }
    }
}

