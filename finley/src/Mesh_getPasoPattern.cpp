
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

#ifdef ESYS_HAVE_PASO
#include "FinleyDomain.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

namespace finley {

paso::SystemMatrixPattern_ptr FinleyDomain::getPasoPattern(
                              bool reducedRowOrder, bool reducedColOrder) const
{
    paso::SystemMatrixPattern_ptr out;
    // make sure that the requested pattern is available
    if (reducedRowOrder) {
        if (reducedColOrder) {
            if (!ReducedReducedPattern)
                ReducedReducedPattern = makePasoPattern(reducedRowOrder, reducedColOrder);
        } else {
            if (!ReducedFullPattern)
                ReducedFullPattern = makePasoPattern(reducedRowOrder, reducedColOrder);
        }
    } else {
        if (reducedColOrder) {
            if (!FullReducedPattern)
                FullReducedPattern = makePasoPattern(reducedRowOrder, reducedColOrder);
        } else {
            if (!FullFullPattern)
                FullFullPattern = makePasoPattern(reducedRowOrder, reducedColOrder);
        }
    }
    if (reducedRowOrder) {
        if (reducedColOrder) {
            out = ReducedReducedPattern;
        } else {
            out = ReducedFullPattern;
        }
    } else {
        if (reducedColOrder) {
            out = FullReducedPattern;
        } else {
            out = FullFullPattern;
        }
    }  
    return out;
}

paso::SystemMatrixPattern_ptr FinleyDomain::makePasoPattern(
                              bool reducedRowOrder, bool reducedColOrder) const
{
    paso::Connector_ptr col_connector, row_connector;
    escript::Distribution_ptr colDistribution, rowDistribution;
  
    dim_t myNumColTargets, myNumRowTargets;
    dim_t numColTargets, numRowTargets;
    const index_t *colTarget, *rowTarget;

    if (reducedColOrder) {
        myNumColTargets = m_nodes->getNumReducedDegreesOfFreedom();
        numColTargets = m_nodes->getNumReducedDegreesOfFreedomTargets();
        colTarget = m_nodes->borrowTargetReducedDegreesOfFreedom();
        colDistribution = m_nodes->reducedDegreesOfFreedomDistribution;
        col_connector = m_nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumColTargets = m_nodes->getNumDegreesOfFreedom();
        numColTargets = m_nodes->getNumDegreesOfFreedomTargets();
        colTarget = m_nodes->borrowTargetDegreesOfFreedom();
        colDistribution = m_nodes->degreesOfFreedomDistribution;
        col_connector = m_nodes->degreesOfFreedomConnector;
    }

    if (reducedRowOrder) {
        myNumRowTargets = m_nodes->getNumReducedDegreesOfFreedom();
        numRowTargets = m_nodes->getNumReducedDegreesOfFreedomTargets();
        rowTarget = m_nodes->borrowTargetReducedDegreesOfFreedom();
        rowDistribution = m_nodes->reducedDegreesOfFreedomDistribution;
        row_connector = m_nodes->reducedDegreesOfFreedomConnector;
    } else {
        myNumRowTargets = m_nodes->getNumDegreesOfFreedom();
        numRowTargets = m_nodes->getNumDegreesOfFreedomTargets();
        rowTarget = m_nodes->borrowTargetDegreesOfFreedom();
        rowDistribution = m_nodes->degreesOfFreedomDistribution;
        row_connector = m_nodes->degreesOfFreedomConnector;
    }
    boost::scoped_array<IndexList> index_list(new IndexList[numRowTargets]);
  
#pragma omp parallel
    {
        // insert contributions from element matrices into columns in indexlist
        IndexList_insertElements(index_list.get(), m_elements, reducedRowOrder,
                                 rowTarget, reducedColOrder, colTarget);
        IndexList_insertElements(index_list.get(), m_faceElements,
                                 reducedRowOrder, rowTarget, reducedColOrder,
                                 colTarget);
        IndexList_insertElements(index_list.get(), m_contactElements,
                                 reducedRowOrder, rowTarget, reducedColOrder,
                                 colTarget);
        IndexList_insertElements(index_list.get(), m_points, reducedRowOrder,
                                 rowTarget, reducedColOrder, colTarget);
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

} // namespace finley

#endif // ESYS_HAVE_PASO

