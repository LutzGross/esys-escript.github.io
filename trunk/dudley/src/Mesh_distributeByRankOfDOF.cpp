
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

#include "DudleyDomain.h"

namespace dudley {

/// redistributes the Nodes and Elements including overlap
/// according to the DOF distribution. It will create an element colouring
/// but will not create any mappings.
void DudleyDomain::distributeByRankOfDOF(const std::vector<index_t>& dofDistribution)
{
    int* mpiRankOfDOF = new int[m_nodes->getNumNodes()];
    m_nodes->assignMPIRankToDOFs(mpiRankOfDOF, dofDistribution);

    // first, the elements are redistributed according to mpiRankOfDOF
    // at the input the Node tables refer to a the local labeling of the nodes
    // while at the output they refer to the global labeling which is rectified
    // in the next step
    m_elements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);
    m_faceElements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);
    m_points->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);

    // this will replace the node file!
    resolveNodeIds();

    // create a local labeling of the DOFs
    const std::pair<index_t,index_t> dofRange(m_nodes->getDOFRange());
    const dim_t len = dofRange.second - dofRange.first + 1;
    // local mask for used nodes
    index_t* localDOF_mask = new index_t[len];
    index_t* localDOF_map = new index_t[m_nodes->getNumNodes()];

#pragma omp parallel for
    for (index_t n = 0; n < len; n++)
        localDOF_mask[n] = -1;

#pragma omp parallel for
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++)
        localDOF_map[n] = -1;

#pragma omp parallel for
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
#ifdef BOUNDS_CHECK
        ESYS_ASSERT(m_nodes->globalDegreesOfFreedom[n] - dofRange.first < len, "BOUNDS_CHECK");
        ESYS_ASSERT(m_nodes->globalDegreesOfFreedom[n] - dofRange.first >= 0, "BOUNDS_CHECK");
#endif
        localDOF_mask[m_nodes->globalDegreesOfFreedom[n] - dofRange.first] = n;
    }

    dim_t numDOFs = 0;
    for (index_t n = 0; n < len; n++) {
        const index_t k = localDOF_mask[n];
        if (k >= 0) {
            localDOF_mask[n] = numDOFs;
            numDOFs++;
        }
    }
#pragma omp parallel for
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
        localDOF_map[n] = localDOF_mask[
                            m_nodes->globalDegreesOfFreedom[n] - dofRange.first];
    }
    // create element coloring
    createColoring(localDOF_map);

    delete[] localDOF_mask;
    delete[] localDOF_map;
    delete[] mpiRankOfDOF;
}

} // namespace dudley

