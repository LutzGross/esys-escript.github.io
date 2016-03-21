
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

#include "Mesh.h"

namespace dudley {

/// redistributes the Nodes and Elements including overlap
/// according to the DOF distribution. It will create an element colouring
/// but will not create any mappings.
void Mesh::distributeByRankOfDOF(const std::vector<index_t>& dofDistribution)
{
    const dim_t numNodes = Nodes->getNumNodes();
    int* mpiRankOfDOF = new int[numNodes];
    Nodes->assignMPIRankToDOFs(mpiRankOfDOF, dofDistribution);

    // first, the elements are redistributed according to mpiRankOfDOF
    // at the input the Node tables refer to a the local labeling of the nodes
    // while at the output they refer to the global labeling which is rectified
    // in the next step
    Elements->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);
    FaceElements->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);
    Points->distributeByRankOfDOF(mpiRankOfDOF, Nodes->Id);

    resolveNodeIds();

    // create a local labeling of the DOFs
    const std::pair<index_t,index_t> dofRange(Nodes->getDOFRange());
    const dim_t len = dofRange.second - dofRange.first + 1;
    // local mask for used nodes
    index_t* localDOF_mask = new index_t[len];
    index_t* localDOF_map = new index_t[numNodes];

#pragma omp parallel for
    for (index_t n = 0; n < len; n++)
        localDOF_mask[n] = -1;

#pragma omp parallel for
    for (index_t n = 0; n < numNodes; n++)
        localDOF_map[n] = -1;

#pragma omp parallel for
    for (index_t n = 0; n < numNodes; n++) {
#ifdef BOUNDS_CHECK
        if ((Nodes->globalDegreesOfFreedom[n] - dofRange.first) >= len
            || (Nodes->globalDegreesOfFreedom[n] - dofRange.first) < 0) {
            printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
            exit(1);
        }
#endif
        localDOF_mask[Nodes->globalDegreesOfFreedom[n] - dofRange.first] = n;
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
    for (index_t n = 0; n < numNodes; n++) {
        localDOF_map[n] = localDOF_mask[
                            Nodes->globalDegreesOfFreedom[n] - dofRange.first];
    }
    // create element coloring
    createColoring(localDOF_map);

    delete[] localDOF_mask;
    delete[] localDOF_map;
    delete[] mpiRankOfDOF;
}

} // namespace dudley

