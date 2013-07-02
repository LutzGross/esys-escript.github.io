
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


/// returns a reference to the matrix pattern
Paso_SystemMatrixPattern* Finley_getPattern(Finley_Mesh *mesh, bool reduce_row_order, bool reduce_col_order)
{
    Paso_SystemMatrixPattern *out=NULL;
    Finley_resetError();
    /* make sure that the requested pattern is available */
    if (reduce_row_order) {
        if (reduce_col_order) {
            if (mesh->ReducedReducedPattern==NULL)
                mesh->ReducedReducedPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
        } else {
            if (mesh->ReducedFullPattern==NULL)
                mesh->ReducedFullPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
        }
    } else {
        if (reduce_col_order) {
            if (mesh->FullReducedPattern==NULL)
                mesh->FullReducedPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
        } else {
            if (mesh->FullFullPattern==NULL)
                mesh->FullFullPattern=Finley_makePattern(mesh,reduce_row_order,reduce_col_order);
        }
    }
    if (Finley_noError()) {
        if (reduce_row_order) {
            if (reduce_col_order) {
                out=Paso_SystemMatrixPattern_getReference(mesh->ReducedReducedPattern);
            } else {
                out=Paso_SystemMatrixPattern_getReference(mesh->ReducedFullPattern);
            }
        } else {
            if (reduce_col_order) {
                out=Paso_SystemMatrixPattern_getReference(mesh->FullReducedPattern);
            } else {
                out=Paso_SystemMatrixPattern_getReference(mesh->FullFullPattern);
            }
        }
    }  
    return out;
}

Paso_SystemMatrixPattern* Finley_makePattern(Finley_Mesh *mesh, bool reduce_row_order, bool reduce_col_order)
{
    Paso_SystemMatrixPattern* out=NULL;
    Paso_Pattern *main_pattern = NULL, *col_couple_pattern=NULL, *row_couple_pattern=NULL;
    Paso_Connector *col_connector, *row_connector;
    Paso_Distribution *colDistribution=NULL, *rowDistribution=NULL;
  
    Finley_resetError();

    int myNumColTargets, myNumRowTargets;
    int numColTargets, numRowTargets;
    const int *colTarget, *rowTarget;

    if (reduce_col_order) {
        myNumColTargets=mesh->Nodes->getNumReducedDegreesOfFreedom();
        numColTargets=mesh->Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        colTarget=mesh->Nodes->borrowTargetReducedDegreesOfFreedom();
        colDistribution=mesh->Nodes->reducedDegreesOfFreedomDistribution;
        col_connector=mesh->Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumColTargets=mesh->Nodes->getNumDegreesOfFreedom();
        numColTargets=mesh->Nodes->degreesOfFreedomMapping.getNumTargets();
        colTarget=mesh->Nodes->borrowTargetDegreesOfFreedom();
        colDistribution=mesh->Nodes->degreesOfFreedomDistribution;
        col_connector=mesh->Nodes->degreesOfFreedomConnector;
    }

    if (reduce_row_order) {
        myNumRowTargets=mesh->Nodes->getNumReducedDegreesOfFreedom();
        numRowTargets=mesh->Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        rowTarget=mesh->Nodes->borrowTargetReducedDegreesOfFreedom();
        rowDistribution=mesh->Nodes->reducedDegreesOfFreedomDistribution;
        row_connector=mesh->Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumRowTargets=mesh->Nodes->getNumDegreesOfFreedom();
        numRowTargets=mesh->Nodes->degreesOfFreedomMapping.getNumTargets();
        rowTarget=mesh->Nodes->borrowTargetDegreesOfFreedom();
        rowDistribution=mesh->Nodes->degreesOfFreedomDistribution;
        row_connector=mesh->Nodes->degreesOfFreedomConnector;
    }
    IndexList* index_list=new IndexList[numRowTargets];
  
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index index_list:
        IndexList_insertElements(index_list, mesh->Elements, reduce_row_order,
                                 rowTarget, reduce_col_order, colTarget);
        IndexList_insertElements(index_list, mesh->FaceElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list, mesh->ContactElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list, mesh->Points, reduce_row_order,
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

    /* if everything is in order we can create the return value */
    if (Finley_noError()) {
        out=Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT,
                rowDistribution, colDistribution, main_pattern,
                col_couple_pattern, row_couple_pattern,
                col_connector, row_connector);
    }
    delete[] index_list;
    Paso_Pattern_free(main_pattern);
    Paso_Pattern_free(col_couple_pattern);
    Paso_Pattern_free(row_couple_pattern);
    Esys_MPIInfo_noError(mesh->MPIInfo);
    return out;
}

