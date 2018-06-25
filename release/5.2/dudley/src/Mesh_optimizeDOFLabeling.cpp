
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "DudleyDomain.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

namespace dudley {

/// optimizes the labeling of the DOFs on each processor
void DudleyDomain::optimizeDOFLabeling(const IndexVector& distribution)
{
    // this method relies on Pattern::reduceBandwidth so requires PASO
    // at the moment
#ifdef ESYS_HAVE_PASO
    const int myRank = m_mpiInfo->rank;
    const int mpiSize = m_mpiInfo->size;
    const index_t myFirstVertex = distribution[myRank];
    const index_t myLastVertex = distribution[myRank+1];
    const dim_t myNumVertices = myLastVertex-myFirstVertex;
    dim_t len = 0;
    for (int p=0; p<mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);

    boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
    boost::scoped_array<index_t> newGlobalDOFID(new index_t[len]);

    // create the adjacency structure xadj and adjncy
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
            myFirstVertex, myLastVertex, m_elements,
            m_nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
            myFirstVertex, myLastVertex, m_faceElements,
            m_nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
            myFirstVertex, myLastVertex, m_points,
            m_nodes->globalDegreesOfFreedom);
    }
    // create the local matrix pattern
    paso::Pattern_ptr pattern = paso::Pattern::fromIndexListArray(0,
            myNumVertices, index_list.get(), myFirstVertex, myLastVertex,
            -myFirstVertex);

    pattern->reduceBandwidth(&newGlobalDOFID[0]);

    // shift new labeling to create a global id
#pragma omp parallel for
    for (index_t i = 0; i < myNumVertices; ++i)
        newGlobalDOFID[i] += myFirstVertex;

    // distribute new labeling to other processors
#ifdef ESYS_MPI
    const int dest = m_mpiInfo->mod_rank(myRank + 1);
    const int source = m_mpiInfo->mod_rank(myRank - 1);
#endif
    int current_rank = myRank;
    for (int p = 0; p < mpiSize; ++p) {
        const index_t firstVertex = distribution[current_rank];
        const index_t lastVertex = distribution[current_rank + 1];
#pragma omp parallel for
        for (index_t i = 0; i < m_nodes->getNumNodes(); ++i) {
            const index_t k = m_nodes->globalDegreesOfFreedom[i];
            if (firstVertex <= k && k < lastVertex) {
                m_nodes->globalDegreesOfFreedom[i] =
                                            newGlobalDOFID[k-firstVertex];
            }
        }

        if (p < mpiSize - 1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_DIM_T, dest,
                                 m_mpiInfo->counter(), source,
                                 m_mpiInfo->counter(), m_mpiInfo->comm, &status);
            m_mpiInfo->incCounter();
#endif
            current_rank = m_mpiInfo->mod_rank(current_rank - 1);
        }
    }
#if 0
    for (index_t i = 0; i < m_nodes->getNumNodes(); ++i)
        std::cout << m_nodes->globalDegreesOfFreedom[i] << " ";
    std::cout << std::endl;
#endif
#endif // ESYS_HAVE_PASO
}

} // namespace dudley

