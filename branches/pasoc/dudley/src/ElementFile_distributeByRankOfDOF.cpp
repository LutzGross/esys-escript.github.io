
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

#include "ElementFile.h"

#include <escript/index.h>

namespace dudley {

void ElementFile::distributeByRankOfDOF(const int* mpiRankOfDOF,
                                        const index_t* nodesId)
{
    const int size = MPIInfo->size;
    if (size > 1) {
#ifdef ESYS_MPI
        const int myRank = MPIInfo->rank;
        int numRequests = 0;
        std::vector<MPI_Request> mpi_requests(8 * size);
        std::vector<MPI_Status> mpi_stati(8 * size);

        // count the number elements that have to be sent to each processor
        // (send_count) and define a new element owner as the processor with
        // the largest number of DOFs and the smallest id
        std::vector<dim_t> send_count(size);
        std::vector<dim_t> recv_count(size);
        int* newOwner = new int[numElements];
#pragma omp parallel
        {
            std::vector<dim_t> loc_proc_mask(size);
            std::vector<dim_t> loc_send_count(size);
#pragma omp for
            for (index_t e = 0; e < numElements; e++) {
                if (Owner[e] == myRank) {
                    newOwner[e] = myRank;
                    loc_proc_mask.assign(size, 0);
                    for (int j = 0; j < numNodes; j++) {
                        const int p = mpiRankOfDOF[Nodes[INDEX2(j, e, numNodes)]];
                        loc_proc_mask[p]++;
                    }
                    dim_t loc_proc_mask_max = 0;
                    for (int p = 0; p < size; ++p) {
                        if (loc_proc_mask[p] > 0)
                            loc_send_count[p]++;
                        if (loc_proc_mask[p] > loc_proc_mask_max) {
                            newOwner[e] = p;
                            loc_proc_mask_max = loc_proc_mask[p];
                        }
                    }
                } else {
                    newOwner[e] = -1;
                }
            }
#pragma omp critical
            {
                for (int p = 0; p < size; ++p)
                    send_count[p] += loc_send_count[p];
            }
        } // end parallel section
        MPI_Alltoall(&send_count[0], 1, MPI_DIM_T, &recv_count[0], 1,
                     MPI_DIM_T, MPIInfo->comm);
        // get the new number of elements for this processor
        dim_t newNumElements = 0;
        dim_t numElementsInBuffer = 0;
        for (int p = 0; p < size; ++p) {
            newNumElements += recv_count[p];
            numElementsInBuffer += send_count[p];
        }

        std::vector<index_t> Id_buffer(numElementsInBuffer);
        std::vector<int> Tag_buffer(numElementsInBuffer);
        std::vector<int> Owner_buffer(numElementsInBuffer);
        std::vector<index_t> Nodes_buffer(numElementsInBuffer * numNodes);
        std::vector<index_t> send_offset(size);
        std::vector<index_t> recv_offset(size);
        std::vector<unsigned char> proc_mask(size);

        // calculate the offsets for the processor buffers
        for (int p = 0; p < size - 1; ++p) {
            recv_offset[p + 1] = recv_offset[p] + recv_count[p];
            send_offset[p + 1] = send_offset[p] + send_count[p];
        }

        send_count.assign(size, 0);
        // copy element into buffers. proc_mask makes sure that an element is
        // copied once only for each processor
        for (index_t e = 0; e < numElements; e++) {
            if (Owner[e] == myRank) {
                proc_mask.assign(size, 1);
                for (int j = 0; j < numNodes; j++) {
                    const int p = mpiRankOfDOF[Nodes[INDEX2(j, e, numNodes)]];
                    if (proc_mask[p]) {
                        const index_t k = send_offset[p] + send_count[p];
                        Id_buffer[k] = Id[e];
                        Tag_buffer[k] = Tag[e];
                        Owner_buffer[k] = newOwner[e];
                        for (int i = 0; i < numNodes; i++)
                            Nodes_buffer[INDEX2(i, k, numNodes)] =
                                         nodesId[Nodes[INDEX2(i, e, numNodes)]];
                        send_count[p]++;
                        proc_mask[p] = 0;
                    }
                }
            }
        }
        // allocate new tables
        allocTable(newNumElements);

        // start to receive new elements
        for (int p = 0; p < size; ++p) {
            if (recv_count[p] > 0) {
                MPI_Irecv(&Id[recv_offset[p]], recv_count[p], MPI_DIM_T, p,
                          MPIInfo->counter() + myRank, MPIInfo->comm,
                          &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Tag[recv_offset[p]], recv_count[p], MPI_INT, p,
                          MPIInfo->counter() + size + myRank, MPIInfo->comm,
                          &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Owner[recv_offset[p]], recv_count[p], MPI_INT, p,
                          MPIInfo->counter() + 2 * size + myRank,
                          MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&Nodes[recv_offset[p] * numNodes],
                          recv_count[p] * numNodes, MPI_DIM_T, p,
                          MPIInfo->counter() + 3 * size + myRank,
                          MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        // now the buffers can be sent away
        for (int p = 0; p < size; ++p) {
            if (send_count[p] > 0) {
                MPI_Issend(&Id_buffer[send_offset[p]], send_count[p],
                           MPI_DIM_T, p, MPIInfo->counter() + p,
                           MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Tag_buffer[send_offset[p]], send_count[p],
                           MPI_INT, p, MPIInfo->counter() + size + p,
                           MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Owner_buffer[send_offset[p]], send_count[p],
                           MPI_INT, p, MPIInfo->counter() + 2 * size + p,
                           MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&Nodes_buffer[send_offset[p] * numNodes],
                           send_count[p] * numNodes, MPI_DIM_T, p,
                           MPIInfo->counter() + 3 * size + p,
                           MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        MPIInfo->incCounter(4 * size);
        // wait for the requests to be finalized
        MPI_Waitall(numRequests, &mpi_requests[0], &mpi_stati[0]);
        delete[] newOwner;
#endif
    } else { // single rank
#pragma omp parallel for
        for (index_t e = 0; e < numElements; e++) {
            Owner[e] = 0;
            for (int i = 0; i < numNodes; i++)
                Nodes[INDEX2(i, e, numNodes)] =
                                     nodesId[Nodes[INDEX2(i, e, numNodes)]];
        }
    }
}

} // namespace dudley

