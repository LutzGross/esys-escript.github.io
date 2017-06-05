
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

#include "NodeFile.h"

namespace dudley {

dim_t NodeFile::createDenseDOFLabeling()
{
    const index_t UNSET_ID = -1, SET_ID = 1;

    // get the global range of DOF IDs
    const std::pair<index_t,index_t> idRange(getGlobalDOFRange());

    // distribute the range of DOF IDs
    std::vector<index_t> distribution(MPIInfo->size + 1);
    dim_t bufferLen = MPIInfo->setDistribution(idRange.first, idRange.second,
                                              &distribution[0]);

    index_t* DOF_buffer = new index_t[bufferLen];
    // fill buffer by the UNSET_ID marker to check if nodes are defined
#pragma omp parallel for
    for (index_t n = 0; n < bufferLen; n++)
        DOF_buffer[n] = UNSET_ID;

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
#ifdef ESYS_MPI
        if (p > 0) { // the initial send can be skipped
            MPI_Sendrecv_replace(DOF_buffer, bufferLen, MPI_DIM_T, dest,
                                 MPIInfo->counter(), source, MPIInfo->counter(),
                                 MPIInfo->comm, &status);
            MPIInfo->incCounter();
        }
#endif
        buffer_rank = MPIInfo->mod_rank(buffer_rank - 1);
        const index_t dof0 = distribution[buffer_rank];
        const index_t dof1 = distribution[buffer_rank + 1];
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            const index_t k = globalDegreesOfFreedom[n];
            if (dof0 <= k && k < dof1) {
                DOF_buffer[k - dof0] = SET_ID;
            }
        }
    }
    // count the entries in the buffer
    const dim_t myDOFs = distribution[MPIInfo->rank + 1] - distribution[MPIInfo->rank];
    // TODO: OMP parallel
    dim_t myNewDOFs = 0;
    for (index_t n = 0; n < myDOFs; ++n) {
        if (DOF_buffer[n] == SET_ID) {
            DOF_buffer[n] = myNewDOFs;
            myNewDOFs++;
        }
    }

    std::vector<index_t> loc_offsets(MPIInfo->size);
    std::vector<index_t> offsets(MPIInfo->size);
    dim_t new_numGlobalDOFs;
    bool* set_new_DOF = new bool[numNodes];

#ifdef ESYS_MPI
    new_numGlobalDOFs = 0;
    loc_offsets[MPIInfo->rank] = myNewDOFs;
    MPI_Allreduce(&loc_offsets[0], &offsets[0], MPIInfo->size, MPI_DIM_T,
                  MPI_SUM, MPIInfo->comm);
    for (int n = 0; n < MPIInfo->size; ++n) {
        loc_offsets[n] = new_numGlobalDOFs;
        new_numGlobalDOFs += offsets[n];
    }
#else
    new_numGlobalDOFs = myNewDOFs;
#endif

#pragma omp parallel
    {
#pragma omp for
        for (index_t n = 0; n < myDOFs; ++n)
            DOF_buffer[n] += loc_offsets[MPIInfo->rank];
#pragma omp for
        for (index_t n = 0; n < numNodes; ++n)
            set_new_DOF[n] = true;
    }

    // now entries are collected from the buffer again by sending them around
    // in a circle
#ifdef ESYS_MPI
    dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
        const index_t dof0 = distribution[buffer_rank];
        const index_t dof1 = distribution[buffer_rank + 1];
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            const index_t k = globalDegreesOfFreedom[n];
            if (set_new_DOF[n] && dof0 <= k && k < dof1) {
                globalDegreesOfFreedom[n] = DOF_buffer[k - dof0];
                set_new_DOF[n] = false;
            }
        }
#ifdef ESYS_MPI
        if (p < MPIInfo->size - 1) { // the last send can be skipped
            MPI_Sendrecv_replace(DOF_buffer, bufferLen, MPI_DIM_T, dest,
                                 MPIInfo->counter(), source, MPIInfo->counter(),
                                 MPIInfo->comm, &status);
            MPIInfo->incCounter();
        }
#endif
        buffer_rank = MPIInfo->mod_rank(buffer_rank - 1);
    }
    delete[] DOF_buffer;
    delete[] set_new_DOF;
    return new_numGlobalDOFs;
}

dim_t NodeFile::createDenseNodeLabeling(std::vector<index_t>& nodeDistribution,
                                   const std::vector<index_t>& dofDistribution)
{
    const index_t UNSET_ID = -1, SET_ID = 1;
    const index_t myFirstDOF = dofDistribution[MPIInfo->rank];
    const index_t myLastDOF = dofDistribution[MPIInfo->rank + 1];

    // find the range of node IDs controlled by me
    index_t min_id = escript::DataTypes::index_t_max();
    index_t max_id = escript::DataTypes::index_t_min();
#pragma omp parallel
    {
        index_t loc_min_id = min_id;
        index_t loc_max_id = max_id;
#pragma omp for
        for (index_t n = 0; n < numNodes; n++) {
            const index_t dof = globalDegreesOfFreedom[n];
            if (myFirstDOF <= dof && dof < myLastDOF) {
                loc_min_id = std::min(loc_min_id, Id[n]);
                loc_max_id = std::max(loc_max_id, Id[n]);
            }
        }
#pragma omp critical
        {
            min_id = std::min(loc_min_id, min_id);
            max_id = std::max(loc_max_id, max_id);
        }
    }
    dim_t myBufferLen = (max_id >= min_id ? max_id - min_id + 1 : 0);
    dim_t bufferLen;

#ifdef ESYS_MPI
    MPI_Allreduce(&myBufferLen, &bufferLen, 1, MPI_DIM_T, MPI_MAX,
                  MPIInfo->comm);
#else
    bufferLen = myBufferLen;
#endif

    const dim_t headerLen = 2;

    index_t* Node_buffer = new index_t[bufferLen + headerLen];
    // mark and count the nodes in use
#pragma omp parallel
    {
#pragma omp for
        for (index_t n = 0; n < bufferLen + headerLen; n++)
            Node_buffer[n] = UNSET_ID;
#pragma omp for
        for (index_t n = 0; n < numNodes; n++) {
            globalNodesIndex[n] = -1;
            const index_t dof = globalDegreesOfFreedom[n];
            if (myFirstDOF <= dof && dof < myLastDOF)
                Node_buffer[Id[n] - min_id + headerLen] = SET_ID;
        }
    }
    dim_t myNewNumNodes = 0;
    for (index_t n = 0; n < myBufferLen; n++) {
        if (Node_buffer[headerLen + n] == SET_ID) {
            Node_buffer[headerLen + n] = myNewNumNodes;
            myNewNumNodes++;
        }
    }
    // make the local number of nodes globally available
#ifdef ESYS_MPI
    MPI_Allgather(&myNewNumNodes, 1, MPI_DIM_T, &nodeDistribution[0], 1,
                  MPI_DIM_T, MPIInfo->comm);
#else
    nodeDistribution[0] = myNewNumNodes;
#endif

    dim_t globalNumNodes = 0;
    for (int p = 0; p < MPIInfo->size; ++p) {
        const dim_t itmp = nodeDistribution[p];
        nodeDistribution[p] = globalNumNodes;
        globalNumNodes += itmp;
    }
    nodeDistribution[MPIInfo->size] = globalNumNodes;

    // offset node buffer
#pragma omp parallel for
    for (index_t n = 0; n < myBufferLen; n++)
        Node_buffer[n + headerLen] += nodeDistribution[MPIInfo->rank];

    // now we send this buffer around to assign global node index
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    Node_buffer[0] = min_id;
    Node_buffer[1] = max_id;
    int buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
        const index_t nodeID0 = Node_buffer[0];
        const index_t nodeID1 = Node_buffer[1];
        const index_t dof0 = dofDistribution[buffer_rank];
        const index_t dof1 = dofDistribution[buffer_rank + 1];
        if (nodeID0 <= nodeID1) {
#pragma omp parallel for
            for (index_t n = 0; n < numNodes; n++) {
                const index_t dof = globalDegreesOfFreedom[n];
                const index_t id = Id[n] - nodeID0;
                if (dof0 <= dof && dof < dof1 && id >= 0 &&
                        id <= nodeID1 - nodeID0)
                    globalNodesIndex[n] = Node_buffer[id + headerLen];
            }
        }
#ifdef ESYS_MPI
        if (p < MPIInfo->size - 1) { // the last send can be skipped
            MPI_Status status;
            MPI_Sendrecv_replace(Node_buffer, bufferLen + headerLen, MPI_DIM_T,
                                 dest, MPIInfo->counter(), source,
                                 MPIInfo->counter(), MPIInfo->comm, &status);
            MPIInfo->incCounter();
        }
#endif
        buffer_rank = MPIInfo->mod_rank(buffer_rank - 1);
    }
    delete[] Node_buffer;
    return globalNumNodes;
}

} // namespace dudley

