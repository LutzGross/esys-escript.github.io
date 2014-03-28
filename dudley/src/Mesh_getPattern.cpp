
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

/* Dudley: Mesh */

/************************************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/************************************************************************************/

/* returns a reference to the matrix pattern                  */

paso::SystemMatrixPattern *Dudley_getPattern(Dudley_Mesh * mesh, bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern *out = NULL;
    Dudley_resetError();
    /* make sure that the requested pattern is available */
    if (reduce_row_order)
    {
        if (reduce_col_order)
        {
            if (mesh->ReducedReducedPattern == NULL)
                mesh->ReducedReducedPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
        else
        {
            if (mesh->ReducedFullPattern == NULL)
                mesh->ReducedFullPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
    }
    else
    {
        if (reduce_col_order)
        {
            if (mesh->FullReducedPattern == NULL)
                mesh->FullReducedPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
        else
        {
            if (mesh->FullFullPattern == NULL)
                mesh->FullFullPattern = Dudley_makePattern(mesh, reduce_row_order, reduce_col_order);
        }
    }
    if (Dudley_noError())
    {
        if (reduce_row_order)
        {
            if (reduce_col_order)
            {
                out = paso::SystemMatrixPattern_getReference(mesh->ReducedReducedPattern);
            }
            else
            {
                out = paso::SystemMatrixPattern_getReference(mesh->ReducedFullPattern);
            }
        }
        else
        {
            if (reduce_col_order)
            {
                out = paso::SystemMatrixPattern_getReference(mesh->FullReducedPattern);
            }
            else
            {
                out = paso::SystemMatrixPattern_getReference(mesh->FullFullPattern);
            }
        }
    }
    return out;
}

paso::SystemMatrixPattern *Dudley_makePattern(Dudley_Mesh * mesh, bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern *out = NULL;
    paso::Pattern *main_pattern = NULL, *col_couple_pattern = NULL, *row_couple_pattern = NULL;
    paso::Connector *col_connector, *row_connector;
    Dudley_NodeMapping *colMap = NULL, *rowMap = NULL;
    paso::Distribution_ptr colDistribution, rowDistribution;

    Dudley_resetError();

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

    IndexListArray index_list(rowMap->numTargets);
    {

#pragma omp parallel
        {
            /*  insert contributions from element matrices into columns index index_list: */
            Dudley_IndexList_insertElements(index_list, mesh->Elements,
                                            reduce_row_order, rowMap->target, reduce_col_order, colMap->target);
            Dudley_IndexList_insertElements(index_list, mesh->FaceElements,
                                            reduce_row_order, rowMap->target, reduce_col_order, colMap->target);
            Dudley_IndexList_insertElements(index_list, mesh->Points,
                                            reduce_row_order, rowMap->target, reduce_col_order, colMap->target);

        }

        /* create pattern */
        main_pattern = paso::Pattern_fromIndexListArray(0,
            rowDistribution->getMyNumComponents(), index_list,
            0, colDistribution->getMyNumComponents(), 0);
        col_couple_pattern = paso::Pattern_fromIndexListArray(0,
            rowDistribution->getMyNumComponents(), index_list,
            colDistribution->getMyNumComponents(), colMap->numTargets,
            -colDistribution->getMyNumComponents());
        row_couple_pattern = paso::Pattern_fromIndexListArray(
            rowDistribution->getMyNumComponents(), rowMap->numTargets,
            index_list, 0, colDistribution->getMyNumComponents(), 0);

        /* if everything is in order we can create the return value */
        if (Dudley_noError())
        {
            out = new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                    rowDistribution, colDistribution, main_pattern,
                    col_couple_pattern, row_couple_pattern, col_connector,
                    row_connector);
        }
    paso::Pattern_free(main_pattern);
    paso::Pattern_free(col_couple_pattern);
    paso::Pattern_free(row_couple_pattern);
    }
    Esys_MPIInfo_noError(mesh->MPIInfo);
    return out;
}

