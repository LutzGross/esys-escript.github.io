
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#ifdef ESYS_HAVE_TRILINOS

#include "DudleyDomain.h"
#include "IndexList.h"

using namespace esys_trilinos;

namespace dudley {

void NodeFile::createTrilinosGraph(const IndexList* indexList)
{
    const index_t* gNI = borrowGlobalNodesIndex();
    const index_t* dofMap = borrowDegreesOfFreedomTarget();

    const index_t myNumTargets = getNumDegreesOfFreedom();
    const index_t numTargets = getNumDegreesOfFreedomTargets();
    const index_t* target = borrowTargetDegreesOfFreedom();
    IndexVector myRows(myNumTargets);
    IndexVector columns(numTargets);

#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t i = 0; i < myRows.size(); i++) {
            myRows[i] = gNI[dofMap[i]];
        }
#pragma omp for
        for (size_t i = 0; i < columns.size(); i++) {
            columns[target[i]] = gNI[i];
        }
    } // end parallel section

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

    TrilinosMap_ptr rowMap(new MapType(getGlobalNumNodes(), myRows,
                0, TeuchosCommFromEsysComm(MPIInfo->comm)));

    TrilinosMap_ptr colMap(new MapType(getGlobalNumNodes(), columns,
                0, TeuchosCommFromEsysComm(MPIInfo->comm)));

    GraphType* graph = new GraphType(rowMap, colMap, rowPtr, colInd);
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("Optimize Storage", true);
    graph->fillComplete(rowMap, rowMap, params);
    m_graph.reset(graph);
}

} // namespace dudley

#endif // ESYS_HAVE_TRILINOS

