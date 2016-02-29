
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

/****************************************************************************/

/* Dudley: Mesh */

/****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Mesh.h"
#include "IndexList.h"
#include <boost/scoped_array.hpp>

namespace dudley {

/* returns a reference to the matrix pattern                  */
paso::SystemMatrixPattern_ptr Dudley_getPattern(Dudley_Mesh* mesh, bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern_ptr out;
    /* make sure that the requested pattern is available */
    if (reduce_row_order) {
        if (reduce_col_order) {
            if (mesh->ReducedReducedPattern == NULL)
                mesh->ReducedReducedPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        } else {
            if (mesh->ReducedFullPattern == NULL)
                mesh->ReducedFullPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
    } else {
        if (reduce_col_order) {
            if (mesh->FullReducedPattern == NULL)
                mesh->FullReducedPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        } else {
            if (mesh->FullFullPattern == NULL)
                mesh->FullFullPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
    }

    if (reduce_row_order) {
        if (reduce_col_order) {
            out = mesh->ReducedReducedPattern;
        } else {
            out = mesh->ReducedFullPattern;
        }
    } else {
        if (reduce_col_order) {
            out = mesh->FullReducedPattern;
        } else {
            out = mesh->FullFullPattern;
        }
    }
    return out;
}

paso::SystemMatrixPattern_ptr Dudley_makePattern(Dudley_Mesh* mesh, bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern_ptr out;
    paso::Pattern_ptr main_pattern, col_couple_pattern, row_couple_pattern;
    paso::Connector_ptr col_connector, row_connector;
    Dudley_NodeMapping *colMap = NULL, *rowMap = NULL;
    paso::Distribution_ptr colDistribution, rowDistribution;

    if (reduce_col_order)
    {
        colMap = mesh->Nodes->reducedDegreesOfFreedomMapping;
        colDistribution = mesh->Nodes->reducedDegreesOfFreedomDistribution;
        col_connector = mesh->Nodes->reducedDegreesOfFreedomConnector;

    }
    else
    {
        colMap = mesh->Nodes->degreesOfFreedomMapping;
        colDistribution = mesh->Nodes->degreesOfFreedomDistribution;
        col_connector = mesh->Nodes->degreesOfFreedomConnector;
    }

    if (reduce_row_order)
    {
        rowMap = mesh->Nodes->reducedDegreesOfFreedomMapping;
        rowDistribution = mesh->Nodes->reducedDegreesOfFreedomDistribution;
        row_connector = mesh->Nodes->reducedDegreesOfFreedomConnector;
    }
    else
    {
        rowMap = mesh->Nodes->degreesOfFreedomMapping;
        rowDistribution = mesh->Nodes->degreesOfFreedomDistribution;
        row_connector = mesh->Nodes->degreesOfFreedomConnector;
    }

    boost::scoped_array<IndexList> index_list(new IndexList[rowMap->numTargets]);

#pragma omp parallel
    {
        /*  insert contributions from element matrices into columns index index_list: */
        Dudley_IndexList_insertElements(index_list.get(), mesh->Elements,
                                        reduce_row_order, rowMap->target, reduce_col_order, colMap->target);
        Dudley_IndexList_insertElements(index_list.get(), mesh->FaceElements,
                                        reduce_row_order, rowMap->target, reduce_col_order, colMap->target);
        Dudley_IndexList_insertElements(index_list.get(), mesh->Points,
                                        reduce_row_order, rowMap->target, reduce_col_order, colMap->target);

    }

    /* create pattern */
    main_pattern = paso::Pattern::fromIndexListArray(0,
        rowDistribution->getMyNumComponents(), index_list.get(),
        0, colDistribution->getMyNumComponents(), 0);
    col_couple_pattern = paso::Pattern::fromIndexListArray(0,
        rowDistribution->getMyNumComponents(), index_list.get(),
        colDistribution->getMyNumComponents(), colMap->numTargets,
        -colDistribution->getMyNumComponents());
    row_couple_pattern = paso::Pattern::fromIndexListArray(
        rowDistribution->getMyNumComponents(), rowMap->numTargets,
        index_list.get(), 0, colDistribution->getMyNumComponents(), 0);

    /* if everything is in order we can create the return value */
    out.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
            rowDistribution, colDistribution, main_pattern,
            col_couple_pattern, row_couple_pattern, col_connector,
            row_connector));
    return out;
}

} // namespace dudley

