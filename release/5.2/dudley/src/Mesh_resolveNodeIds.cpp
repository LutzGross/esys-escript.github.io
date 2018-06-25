
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
#include "Util.h"

namespace dudley {

void DudleyDomain::resolveNodeIds()
{
    // find the minimum and maximum id used by elements
    index_t min_id = escript::DataTypes::index_t_max();
    index_t max_id = -escript::DataTypes::index_t_max();
    std::pair<index_t,index_t> range(m_elements->getNodeRange());
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
    range = m_faceElements->getNodeRange();
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
    range = m_points->getNodeRange();
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
#ifdef Dudley_TRACE
    index_t global_min_id, global_max_id;
#ifdef ESYS_MPI
    index_t id_range[2], global_id_range[2];
    id_range[0] = -min_id;
    id_range[1] = max_id;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_DIM_T, MPI_MAX, m_mpiInfo->comm);
    global_min_id = -global_id_range[0];
    global_max_id = global_id_range[1];
#else
    global_min_id = min_id;
    global_max_id = max_id;
#endif
    printf("Node id range used by elements is %d:%d\n", global_min_id, global_max_id);
#endif
    if (min_id > max_id) {
        max_id = -1;
        min_id = 0;
    }

    // allocate mappings for new local node labeling to global node labeling
    // (newLocalToGlobalNodeLabels) and global node labeling to the new local
    // node labeling (globalToNewLocalNodeLabels[i-min_id] is the new local id
    // of global node i)
    index_t len = (max_id >= min_id) ? max_id - min_id + 1 : 0;

    // mark the nodes referred by elements in usedMask
    std::vector<short> usedMask(len, -1);
    markNodes(usedMask, min_id);

    // create a local labeling newLocalToGlobalNodeLabels of the local nodes
    // by packing the mask usedMask
    std::vector<index_t> newLocalToGlobalNodeLabels =  util::packMask(usedMask);
    const dim_t newNumNodes = newLocalToGlobalNodeLabels.size();

    usedMask.clear();

    // invert the new labeling and shift the index newLocalToGlobalNodeLabels
    // to global node IDs
    index_t* globalToNewLocalNodeLabels = new index_t[len];

#pragma omp parallel for
    for (index_t n = 0; n < newNumNodes; n++) {
#ifdef BOUNDS_CHECK
        ESYS_ASSERT(newLocalToGlobalNodeLabels[n] < len, "BOUNDS_CHECK");
        ESYS_ASSERT(newLocalToGlobalNodeLabels[n] >= 0, "BOUNDS_CHECK");
#endif
        globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]] = n;
        newLocalToGlobalNodeLabels[n] += min_id;
    }
    // create a new node file
    NodeFile* newNodeFile = new NodeFile(getDim(), m_mpiInfo);
    newNodeFile->allocTable(newNumNodes);
    if (len)
        newNodeFile->gather_global(&newLocalToGlobalNodeLabels[0], m_nodes);
    else
        newNodeFile->gather_global(NULL, m_nodes);

    delete m_nodes;
    m_nodes = newNodeFile;
    // relabel nodes of the elements
    relabelElementNodes(globalToNewLocalNodeLabels, min_id);
    delete[] globalToNewLocalNodeLabels;
}

} // namespace dudley

