
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

  Finley: Mesh

*****************************************************************************/

#include "Mesh.h"
#include "IndexList.h"

namespace finley {

/// returns a reference to the matrix pattern
Paso_SystemMatrixPattern* Mesh::getPattern(bool reduce_row_order, bool reduce_col_order)
{
    Paso_SystemMatrixPattern *out=NULL;
    resetError();
    /* make sure that the requested pattern is available */
    if (reduce_row_order) {
        if (reduce_col_order) {
            if (ReducedReducedPattern==NULL)
                ReducedReducedPattern=makePattern(reduce_row_order,reduce_col_order);
        } else {
            if (ReducedFullPattern==NULL)
                ReducedFullPattern=makePattern(reduce_row_order,reduce_col_order);
        }
    } else {
        if (reduce_col_order) {
            if (FullReducedPattern==NULL)
                FullReducedPattern=makePattern(reduce_row_order,reduce_col_order);
        } else {
            if (FullFullPattern==NULL)
                FullFullPattern=makePattern(reduce_row_order,reduce_col_order);
        }
    }
    if (noError()) {
        if (reduce_row_order) {
            if (reduce_col_order) {
                out=Paso_SystemMatrixPattern_getReference(ReducedReducedPattern);
            } else {
                out=Paso_SystemMatrixPattern_getReference(ReducedFullPattern);
            }
        } else {
            if (reduce_col_order) {
                out=Paso_SystemMatrixPattern_getReference(FullReducedPattern);
            } else {
                out=Paso_SystemMatrixPattern_getReference(FullFullPattern);
            }
        }
    }  
    return out;
}

Paso_SystemMatrixPattern* Mesh::makePattern(bool reduce_row_order, bool reduce_col_order)
{
    Paso_SystemMatrixPattern* out=NULL;
    Paso_Pattern *main_pattern = NULL, *col_couple_pattern=NULL, *row_couple_pattern=NULL;
    Paso_Connector *col_connector, *row_connector;
    Paso_Distribution *colDistribution=NULL, *rowDistribution=NULL;
  
    resetError();

    int myNumColTargets, myNumRowTargets;
    int numColTargets, numRowTargets;
    const int *colTarget, *rowTarget;

    if (reduce_col_order) {
        myNumColTargets=Nodes->getNumReducedDegreesOfFreedom();
        numColTargets=Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        colTarget=Nodes->borrowTargetReducedDegreesOfFreedom();
        colDistribution=Nodes->reducedDegreesOfFreedomDistribution;
        col_connector=Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumColTargets=Nodes->getNumDegreesOfFreedom();
        numColTargets=Nodes->degreesOfFreedomMapping.getNumTargets();
        colTarget=Nodes->borrowTargetDegreesOfFreedom();
        colDistribution=Nodes->degreesOfFreedomDistribution;
        col_connector=Nodes->degreesOfFreedomConnector;
    }

    if (reduce_row_order) {
        myNumRowTargets=Nodes->getNumReducedDegreesOfFreedom();
        numRowTargets=Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        rowTarget=Nodes->borrowTargetReducedDegreesOfFreedom();
        rowDistribution=Nodes->reducedDegreesOfFreedomDistribution;
        row_connector=Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumRowTargets=Nodes->getNumDegreesOfFreedom();
        numRowTargets=Nodes->degreesOfFreedomMapping.getNumTargets();
        rowTarget=Nodes->borrowTargetDegreesOfFreedom();
        rowDistribution=Nodes->degreesOfFreedomDistribution;
        row_connector=Nodes->degreesOfFreedomConnector;
    }
    IndexList* index_list=new IndexList[numRowTargets];
  
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index index_list:
        IndexList_insertElements(index_list, Elements, reduce_row_order,
                                 rowTarget, reduce_col_order, colTarget);
        IndexList_insertElements(index_list, FaceElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list, ContactElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list, Points, reduce_row_order,
                                 rowTarget, reduce_col_order, colTarget);
    }
 
    /* create pattern */
    main_pattern=IndexList_createPattern(0, myNumRowTargets, index_list, 0,
                                         myNumColTargets, 0);
    col_couple_pattern=IndexList_createPattern(0, myNumRowTargets, index_list,
                                               myNumColTargets, numColTargets,
                                               -myNumColTargets);
    row_couple_pattern=IndexList_createPattern(myNumRowTargets, numRowTargets,
                                            index_list, 0, myNumColTargets, 0);

    // if everything is in order we can create the return value
    if (noError()) {
        out=Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT,
                rowDistribution, colDistribution, main_pattern,
                col_couple_pattern, row_couple_pattern,
                col_connector, row_connector);
    }
    delete[] index_list;
    Paso_Pattern_free(main_pattern);
    Paso_Pattern_free(col_couple_pattern);
    Paso_Pattern_free(row_couple_pattern);
    Esys_MPIInfo_noError(MPIInfo);
    return out;
}

} // namespace finley

