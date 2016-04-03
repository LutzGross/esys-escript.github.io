
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

namespace dudley {

/// prepares the mesh for further use
void Mesh::prepare(bool optimize)
{
    // first step is to distribute the elements according to a global
    // distribution of DOF
    std::vector<index_t> distribution(MPIInfo->size+1);

    // first we create dense labeling for the DOFs
    dim_t newGlobalNumDOFs = Nodes->createDenseDOFLabeling();

    // create a distribution of the global DOFs and determine the MPI rank
    // controlling the DOFs on this processor
    MPIInfo->setDistribution(0, newGlobalNumDOFs - 1, &distribution[0]);

    // now the mesh is re-distributed according to the distribution vector
    // this will redistribute the Nodes and Elements including overlap and
    // will create an element colouring but will not create any mappings
    // (see later in this function)
    distributeByRankOfDOF(distribution);

    // at this stage we are able to start an optimization of the DOF
    // distribution using ParaMetis. On return distribution is altered and
    // new DOF IDs have been assigned
    if (optimize && MPIInfo->size > 1) {
        optimizeDOFDistribution(distribution);
        distributeByRankOfDOF(distribution);
    }
    // the local labelling of the degrees of freedom is optimized
    if (optimize) {
        optimizeDOFLabeling(distribution);
    }

    // rearrange elements with the aim of bringing elements closer to memory
    // locations of the nodes (distributed shared memory!):
    optimizeElementOrdering();

    // create the global indices
    std::vector<index_t> node_distribution(MPIInfo->size + 1);

    Nodes->createDenseNodeLabeling(node_distribution, distribution);
    // create the missing mappings
    Nodes->createNodeMappings(distribution, node_distribution);

    updateTagList();
}

/// tries to reduce the number of colours for all element files
void Mesh::createColoring(const index_t* node_localDOF_map)
{
    Elements->createColoring(Nodes->getNumNodes(), node_localDOF_map);
    FaceElements->createColoring(Nodes->getNumNodes(), node_localDOF_map);
    Points->createColoring(Nodes->getNumNodes(), node_localDOF_map);
}

/// redistributes elements to minimize communication during assemblage
void Mesh::optimizeElementOrdering()
{
    Elements->optimizeOrdering();
    FaceElements->optimizeOrdering();
    Points->optimizeOrdering();
}

/// regenerates list of tags in use for node file and element files
void Mesh::updateTagList()
{
    Nodes->updateTagList();
    Elements->updateTagList();
    FaceElements->updateTagList();
    Points->updateTagList();
}

} // namespace dudley

