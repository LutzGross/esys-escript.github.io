
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "Mesh.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

namespace finley {

#ifdef ESYS_HAVE_PASO
paso::SystemMatrixPattern_ptr Mesh::getPasoPattern(bool reduce_row_order, bool reduce_col_order)
{
    paso::SystemMatrixPattern_ptr out;
    // make sure that the requested pattern is available
    if (reduce_row_order) {
        if (reduce_col_order) {
            if (!ReducedReducedPattern)
                ReducedReducedPattern = makePasoPattern(reduce_row_order, reduce_col_order);
        } else {
            if (!ReducedFullPattern)
                ReducedFullPattern = makePasoPattern(reduce_row_order, reduce_col_order);
        }
    } else {
        if (reduce_col_order) {
            if (!FullReducedPattern)
                FullReducedPattern = makePasoPattern(reduce_row_order, reduce_col_order);
        } else {
            if (!FullFullPattern)
                FullFullPattern = makePasoPattern(reduce_row_order, reduce_col_order);
        }
    }
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
    return out;
}

paso::SystemMatrixPattern_ptr Mesh::makePasoPattern(bool reduce_row_order, bool reduce_col_order) const
{
    paso::Connector_ptr col_connector, row_connector;
    escript::Distribution_ptr colDistribution, rowDistribution;
  
    dim_t myNumColTargets, myNumRowTargets;
    dim_t numColTargets, numRowTargets;
    const index_t *colTarget, *rowTarget;

    if (reduce_col_order) {
        myNumColTargets = Nodes->getNumReducedDegreesOfFreedom();
        numColTargets = Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        colTarget = Nodes->borrowTargetReducedDegreesOfFreedom();
        colDistribution = Nodes->reducedDegreesOfFreedomDistribution;
        col_connector = Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumColTargets = Nodes->getNumDegreesOfFreedom();
        numColTargets = Nodes->degreesOfFreedomMapping.getNumTargets();
        colTarget = Nodes->borrowTargetDegreesOfFreedom();
        colDistribution = Nodes->degreesOfFreedomDistribution;
        col_connector = Nodes->degreesOfFreedomConnector;
    }

    if (reduce_row_order) {
        myNumRowTargets = Nodes->getNumReducedDegreesOfFreedom();
        numRowTargets = Nodes->reducedDegreesOfFreedomMapping.getNumTargets();
        rowTarget = Nodes->borrowTargetReducedDegreesOfFreedom();
        rowDistribution = Nodes->reducedDegreesOfFreedomDistribution;
        row_connector = Nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumRowTargets = Nodes->getNumDegreesOfFreedom();
        numRowTargets = Nodes->degreesOfFreedomMapping.getNumTargets();
        rowTarget = Nodes->borrowTargetDegreesOfFreedom();
        rowDistribution = Nodes->degreesOfFreedomDistribution;
        row_connector = Nodes->degreesOfFreedomConnector;
    }
    boost::scoped_array<IndexList> index_list(new IndexList[numRowTargets]);
  
#pragma omp parallel
    {
        // insert contributions from element matrices into columns in indexlist
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

    // create pattern
    paso::Pattern_ptr mainPattern(paso::Pattern::fromIndexListArray(0,
              myNumRowTargets, index_list.get(), 0, myNumColTargets, 0));
    paso::Pattern_ptr colCouplePattern(paso::Pattern::fromIndexListArray(0,
              myNumRowTargets, index_list.get(), myNumColTargets,
              numColTargets, -myNumColTargets));
    paso::Pattern_ptr rowCouplePattern(paso::Pattern::fromIndexListArray(
              myNumRowTargets, numRowTargets, index_list.get(), 0,
              myNumColTargets, 0));

    paso::SystemMatrixPattern_ptr out(new paso::SystemMatrixPattern(
                MATRIX_FORMAT_DEFAULT, rowDistribution, colDistribution,
                mainPattern, colCouplePattern, rowCouplePattern,
                col_connector, row_connector));
    return out;
}
#endif // ESYS_HAVE_PASO

} // namespace finley

