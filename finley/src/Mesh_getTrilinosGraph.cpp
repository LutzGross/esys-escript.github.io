
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifdef ESYS_HAVE_TRILINOS
#include "FinleyDomain.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

using namespace esys_trilinos;

namespace finley {

const_TrilinosGraph_ptr FinleyDomain::getTrilinosGraph(bool reducedOrder) const
{
    const_TrilinosGraph_ptr out;
    // make sure that the requested graph is available
    if (reducedOrder) {
        if (m_reducedGraph.is_null())
            m_reducedGraph.reset(createTrilinosGraph(reducedOrder));
        out = m_reducedGraph;
    } else {
        if (m_fullGraph.is_null())
            m_fullGraph.reset(createTrilinosGraph(reducedOrder));
        out = m_fullGraph;
    }
    return out;
}

GraphType* FinleyDomain::createTrilinosGraph(bool reducedOrder) const
{
    index_t myNumTargets;
    index_t numTargets;
    const index_t* target;
    const_TrilinosMap_ptr rowMap;
    const_TrilinosMap_ptr colMap;
    if (reducedOrder) {
        myNumTargets = m_nodes->getNumReducedDegreesOfFreedom();
        numTargets = m_nodes->getNumReducedDegreesOfFreedomTargets();
        target = m_nodes->borrowTargetReducedDegreesOfFreedom();
        rowMap = m_nodes->trilinosReducedRowMap;
        colMap = m_nodes->trilinosReducedColMap;
    } else {
        myNumTargets = m_nodes->getNumDegreesOfFreedom();
        numTargets = m_nodes->getNumDegreesOfFreedomTargets();
        target = m_nodes->borrowTargetDegreesOfFreedom();
        rowMap = m_nodes->trilinosRowMap;
        colMap = m_nodes->trilinosColMap;
    }

    boost::scoped_array<IndexList> indexList(new IndexList[numTargets]);

#pragma omp parallel
    {
        // insert contributions from element matrices into columns in
        // index list
        IndexList_insertElements(indexList.get(), m_elements, reducedOrder,
                                 target, reducedOrder, target);
        IndexList_insertElements(indexList.get(), m_faceElements,
                                 reducedOrder, target, reducedOrder, target);
        IndexList_insertElements(indexList.get(), m_contactElements,
                                 reducedOrder, target, reducedOrder, target);
        IndexList_insertElements(indexList.get(), m_points, reducedOrder,
                                 target, reducedOrder, target);
    }

    Teuchos::ArrayRCP<size_t> rowPtr(myNumTargets + 1);
    for (size_t i = 0; i < myNumTargets; i++) {
        rowPtr[i+1] = rowPtr[i] + indexList[i].count(0, numTargets);
    }

    Teuchos::ArrayRCP<LO> colInd(rowPtr[myNumTargets]);

#pragma omp parallel for
    for (index_t i = 0; i < myNumTargets; i++) {
        indexList[i].toArray(&colInd[rowPtr[i]], 0, numTargets, 0);
        std::sort(&colInd[rowPtr[i]], &colInd[rowPtr[i+1]]);
    }

    GraphType* graph = new GraphType(rowMap, colMap, rowPtr, colInd);
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("Optimize Storage", true);
    graph->fillComplete(rowMap, rowMap, params);
    return graph;
}

} // namespace finley

#endif // ESYS_HAVE_TRILINOS

