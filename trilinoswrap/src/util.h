
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESYS_TRILINOSWRAP_UTIL_H__
#define __ESYS_TRILINOSWRAP_UTIL_H__

#include <escript/EsysException.h>
#include <escript/SolverOptions.h>

#include <Teuchos_ParameterList.hpp>

#include <boost/python.hpp>
#include <boost/python/dict.hpp>

#include <string>

//#ifdef ESYS_TRILINOS_14
#include "Tpetra_CrsMatrix.hpp"
//#endif

namespace esys_trilinos {

namespace util {

template<typename T>
void extractParamIfSet(const std::string& name,
                       const boost::python::dict& pyDict,
                       Teuchos::ParameterList& params)
{
    if (pyDict.has_key(name)) {
        boost::python::object bpo = pyDict.get(name);
        if (boost::python::extract<T>(bpo).check()) {
            T val = boost::python::extract<T>(bpo);
            params.set(name, val);
        } else {
            throw escript::ValueError("Wrong type for option " + name);
        }
    }
}

inline
TrilinosGraph_ptr unrollCrsGraph(const_TrilinosGraph_ptr graph, int blockSize)
{
    // template type double is arbitrary
    MapType cpm = BlockVectorType<double>::makePointMap(*graph->getColMap(), blockSize);
    MapType rpm = BlockVectorType<double>::makePointMap(*graph->getRowMap(), blockSize);
    TrilinosMap_ptr colPointMap(new MapType(cpm));
    TrilinosMap_ptr rowPointMap(new MapType(rpm));

    const LO numMatrixRows = graph->getRowMap()->getLocalNumElements();
    const LO numUnrolledRows = rpm.getLocalNumElements();
    Teuchos::ArrayRCP<size_t> rowPtr(numUnrolledRows + 1);
    Teuchos::ArrayRCP<GO> colInd(graph->getLocalNumEntries() * blockSize * blockSize);
    for (LO row = 0; row < numMatrixRows; row++) {
        size_t numColumns = graph->getNumEntriesInLocalRow(row);
        Tpetra::CrsMatrix<>::nonconst_local_inds_host_view_type indices("indices",numColumns);
        graph->getLocalRowCopy(row, indices, numColumns);
        for (int b = 0; b < blockSize; b++) {
            for (size_t c = 0; c < numColumns; c++) {
                for (int cb = 0; cb < blockSize; cb++) {
                    colInd[rowPtr[row * blockSize + b] + c * blockSize + cb] =
                        indices[c] * blockSize + cb;
                }
            }
            rowPtr[row * blockSize + b + 1] = rowPtr[row * blockSize + b]
                    + numColumns * blockSize;
        }
    }

    GraphType* unrolledGraph = new GraphType(rowPointMap, colPointMap, rowPtr, colInd);

    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("Optimize Storage", true);
    unrolledGraph->fillComplete(rowPointMap, rowPointMap, params);
    return Teuchos::rcp(unrolledGraph);
}

} // namespace util
} // namespace esys_trilinos

#endif // __ESYS_TRILINOSWRAP_UTIL_H__

