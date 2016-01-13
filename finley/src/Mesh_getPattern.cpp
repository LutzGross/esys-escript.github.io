
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


/****************************************************************************

  Finley: Mesh

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "Mesh.h"
#include "IndexList.h"
#include <boost/scoped_array.hpp>

namespace finley {

/// returns a reference to the matrix pattern
paso::SystemMatrixPattern_ptr Mesh::getPattern(bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern_ptr out;
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
                out = ReducedReducedPattern;
            } else {
                out = ReducedFullPattern;
            }
        } else {
            if (reduce_col_order) {
                out = FullReducedPattern;
            } else {
                out = FullFullPattern;
            }
        }
    }  
    return out;
}

paso::SystemMatrixPattern_ptr Mesh::makePattern(bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern_ptr out;
    paso::Connector_ptr col_connector, row_connector;
    paso::Distribution_ptr colDistribution, rowDistribution;
  
    resetError();

    int myNumColTargets, myNumRowTargets;
    int numColTargets, numRowTargets;
    const index_t *colTarget, *rowTarget;

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
    boost::scoped_array<IndexList> index_list(new IndexList[numRowTargets]);
  
#pragma omp parallel
    {
        // insert contributions from element matrices into columns in indexlist:
        IndexList_insertElements(index_list.get(), Elements, reduce_row_order,
                                 rowTarget, reduce_col_order, colTarget);
        IndexList_insertElements(index_list.get(), FaceElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list.get(), ContactElements,
                                 reduce_row_order, rowTarget, reduce_col_order,
                                 colTarget);
        IndexList_insertElements(index_list.get(), Points, reduce_row_order,
                                 rowTarget, reduce_col_order, colTarget);
    }
 
    /* create pattern */
    paso::Pattern_ptr main_pattern, col_couple_pattern, row_couple_pattern;
    main_pattern=paso::Pattern::fromIndexListArray(
            0, myNumRowTargets, index_list.get(), 0, myNumColTargets, 0);
    col_couple_pattern=paso::Pattern::fromIndexListArray(
            0, myNumRowTargets, index_list.get(), myNumColTargets,
            numColTargets, -myNumColTargets);
    row_couple_pattern=paso::Pattern::fromIndexListArray(
            myNumRowTargets, numRowTargets, index_list.get(), 0, myNumColTargets, 0);

    // if everything is in order we can create the return value
    if (noError()) {
        out.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                rowDistribution, colDistribution, main_pattern,
                col_couple_pattern, row_couple_pattern,
                col_connector, row_connector));
    }
    Esys_MPIInfo_noError(MPIInfo);
    return out;
}

} // namespace finley

