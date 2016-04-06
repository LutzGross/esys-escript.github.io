
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


/****************************************************************************

  Finley: Mesh

*****************************************************************************/

#ifdef ESYS_HAVE_TRILINOS

#include "Mesh.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

using namespace esys_trilinos;

namespace finley {

typedef std::vector<index_t> IndexVector;

esys_trilinos::const_TrilinosGraph_ptr Mesh::createTrilinosGraph() const
{
    const index_t numColTargets = Nodes->degreesOfFreedomMapping.getNumTargets();
    const index_t* colTarget = Nodes->borrowTargetDegreesOfFreedom();
    const index_t* gNI = Nodes->borrowGlobalNodesIndex();
    const IndexVector& dofMap = Nodes->degreesOfFreedomMapping.map;

    const index_t myNumRowTargets = Nodes->getNumDegreesOfFreedom();
    const index_t numRowTargets = Nodes->degreesOfFreedomMapping.getNumTargets();
    const index_t* rowTarget = Nodes->borrowTargetDegreesOfFreedom();
    boost::scoped_array<IndexList> index_list(new IndexList[numRowTargets]);
    IndexVector myRows(myNumRowTargets);
    IndexVector columns(numColTargets);

#pragma omp parallel
    {
        // insert contributions from element matrices into columns in indexlist:
        IndexList_insertElements(index_list.get(), Elements, false,
                                 rowTarget, false, colTarget);
        IndexList_insertElements(index_list.get(), FaceElements,
                                 false, rowTarget, false, colTarget);
        IndexList_insertElements(index_list.get(), ContactElements,
                                 false, rowTarget, false, colTarget);
        IndexList_insertElements(index_list.get(), Points, false,
                                 rowTarget, false, colTarget);

#pragma omp for nowait
        for (size_t i=0; i<myRows.size(); i++) {
            myRows[i] = gNI[dofMap[i]];
        }
#pragma omp for
        for (size_t i=0; i<columns.size(); i++) {
            columns[colTarget[i]] = gNI[i];
        }
    } // end parallel section

    Teuchos::ArrayRCP<size_t> rowPtr(myNumRowTargets+1);
    for (size_t i=0; i < myNumRowTargets; i++) {
        rowPtr[i+1] = rowPtr[i] + index_list[i].count(0, numColTargets);
    }

    Teuchos::ArrayRCP<LO> colInd(rowPtr[myNumRowTargets]);

#pragma omp parallel for
    for (index_t i=0; i < myNumRowTargets; i++) {
        index_list[i].toArray(&colInd[rowPtr[i]], 0, numColTargets, 0);
        std::sort(&colInd[rowPtr[i]], &colInd[rowPtr[i+1]]);
    }

    TrilinosMap_ptr rowMap(new MapType(Nodes->getGlobalNumNodes(), myRows,
                0, TeuchosCommFromEsysComm(MPIInfo->comm)));

    TrilinosMap_ptr colMap(new MapType(Nodes->getGlobalNumNodes(), columns,
                0, TeuchosCommFromEsysComm(MPIInfo->comm)));

    TrilinosGraph_ptr graph(new GraphType(rowMap, colMap, rowPtr, colInd));
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("Optimize Storage", true);
    graph->fillComplete(rowMap, rowMap, params);
    return graph;
}

} // namespace finley

#endif // ESYS_HAVE_TRILINOS

