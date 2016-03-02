
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

/****************************************************************************

 Dudley: ElementFile: this will redistribute the Elements including overlap by

*****************************************************************************/

#include "ElementFile.h"

namespace dudley {

void Dudley_ElementFile_distributeByRankOfDOF(Dudley_ElementFile* self, int* mpiRankOfDOF, index_t* Id)
{
    if (self == NULL)
        return;
    dim_t e, i;
    int myRank = self->MPIInfo->rank;
    dim_t NN = self->numNodes;
    dim_t size = self->MPIInfo->size;
    if (size > 1)
    {
#ifdef ESYS_MPI
        int p, *Owner_buffer = NULL, loc_proc_mask_max;
        dim_t j, *send_count = NULL, *recv_count = NULL, *newOwner = NULL;
        dim_t *loc_proc_mask = NULL, *loc_send_count = NULL;
        dim_t newNumElements, numElementsInBuffer;
        index_t *send_offset = NULL, *recv_offset = NULL, *Id_buffer = NULL, *Tag_buffer = NULL, *Nodes_buffer = NULL, k;
        bool *proc_mask = NULL;
        size_t size_size = size * sizeof(dim_t);
        dim_t numRequests = 0;
        MPI_Request *mpi_requests = NULL;
        MPI_Status *mpi_stati = NULL;
        mpi_requests = new  MPI_Request[8 * size];
        mpi_stati = new  MPI_Status[8 * size];
        // count the number elements that have to be sent to each processor
        // (send_count) and define a new element owner as the processor with
        // the largest number of DOFs and the smallest id
        send_count = new dim_t[size];
        recv_count = new dim_t[size];
        newOwner = new int[self->numElements];
        memset(send_count, 0, size_size);
#pragma omp parallel private(p,loc_proc_mask,loc_send_count)
        {
            loc_proc_mask = new dim_t[size];
            loc_send_count = new dim_t[size];
            memset(loc_send_count, 0, size_size);
#pragma omp for private(e,j,loc_proc_mask_max) schedule(static)
            for (e = 0; e < self->numElements; e++)
            {
                if (self->Owner[e] == myRank)
                {
                    newOwner[e] = myRank;
                    memset(loc_proc_mask, 0, size_size);
                    for (j = 0; j < NN; j++)
                    {
                        p = mpiRankOfDOF[self->Nodes[INDEX2(j, e, NN)]];
                        loc_proc_mask[p]++;
                    }
                    loc_proc_mask_max = 0;
                    for (p = 0; p < size; ++p)
                    {
                        if (loc_proc_mask[p] > 0)
                            loc_send_count[p]++;
                        if (loc_proc_mask[p] > loc_proc_mask_max)
                        {
                            newOwner[e] = p;
                            loc_proc_mask_max = loc_proc_mask[p];
                        }
                    }
                }
                else
                {
                    newOwner[e] = -1;
                }
            }
#pragma omp critical
            {
                for (p = 0; p < size; ++p)
                    send_count[p] += loc_send_count[p];
            }
            delete[] loc_proc_mask;
            delete[] loc_send_count;
        }
        MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, self->MPIInfo->comm);
        /* get the new number of elements for this processor */
        newNumElements = 0;
        for (p = 0; p < size; ++p)
            newNumElements += recv_count[p];

        /* get the new number of elements for this processor */
        numElementsInBuffer = 0;
        for (p = 0; p < size; ++p)
            numElementsInBuffer += send_count[p];
        /* allocate buffers */
        Id_buffer = new  index_t[numElementsInBuffer];
        Tag_buffer = new  index_t[numElementsInBuffer];
        Owner_buffer = new  int[numElementsInBuffer];
        Nodes_buffer = new  index_t[numElementsInBuffer * NN];
        send_offset = new  index_t[size];
        recv_offset = new  index_t[size];
        proc_mask = new  bool[size];

        /* calculate the offsets for the processor buffers */
        recv_offset[0] = 0;
        for (p = 0; p < size - 1; ++p)
            recv_offset[p + 1] = recv_offset[p] + recv_count[p];
        send_offset[0] = 0;
        for (p = 0; p < size - 1; ++p)
            send_offset[p + 1] = send_offset[p] + send_count[p];

        memset(send_count, 0, size_size);
        /* copy element into buffers. proc_mask makes sure that an 
         * element is copied once only for each processor */
        for (e = 0; e < self->numElements; e++)
        {
            if (self->Owner[e] == myRank)
            {
                memset(proc_mask, true, size*sizeof(bool));
                for (j = 0; j < NN; j++)
                {
                    p = mpiRankOfDOF[self->Nodes[INDEX2(j, e, NN)]];
                    if (proc_mask[p])
                    {
                        k = send_offset[p] + send_count[p];
                        Id_buffer[k] = self->Id[e];
                        Tag_buffer[k] = self->Tag[e];
                        Owner_buffer[k] = newOwner[e];
                        for (i = 0; i < NN; i++)
                            Nodes_buffer[INDEX2(i, k, NN)] = Id[self->Nodes[INDEX2(i, e, NN)]];
                        send_count[p]++;
                        proc_mask[p] = false;
                    }
                }
            }
        }
        /* allocate new tables */
        Dudley_ElementFile_allocTable(self, newNumElements);

        /* start to receive new elements */
        for (p = 0; p < size; ++p)
        {
            if (recv_count[p] > 0)
            {
                MPI_Irecv(&(self->Id[recv_offset[p]]), recv_count[p],
                          MPI_INT, p, self->MPIInfo->counter() + myRank,
                          self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&(self->Tag[recv_offset[p]]), recv_count[p],
                          MPI_INT, p, self->MPIInfo->counter() + size + myRank,
                          self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&(self->Owner[recv_offset[p]]), recv_count[p],
                          MPI_INT, p, self->MPIInfo->counter() + 2 * size + myRank,
                          self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Irecv(&(self->Nodes[recv_offset[p] * NN]), recv_count[p] * NN,
                          MPI_INT, p, self->MPIInfo->counter() + 3 * size + myRank,
                          self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        /* now the buffers can be send away */
        for (p = 0; p < size; ++p)
        {
            if (send_count[p] > 0)
            {
                MPI_Issend(&(Id_buffer[send_offset[p]]), send_count[p],
                           MPI_INT, p, self->MPIInfo->counter() + p,
                           self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&(Tag_buffer[send_offset[p]]), send_count[p],
                           MPI_INT, p, self->MPIInfo->counter() + size + p,
                           self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&(Owner_buffer[send_offset[p]]), send_count[p],
                           MPI_INT, p, self->MPIInfo->counter() + 2 * size + p,
                           self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
                MPI_Issend(&(Nodes_buffer[send_offset[p] * NN]), send_count[p] * NN,
                           MPI_INT, p, self->MPIInfo->counter() + 3 * size + p,
                           self->MPIInfo->comm, &mpi_requests[numRequests]);
                numRequests++;
            }
        }
        /* wait for the requests to be finalized */
        self->MPIInfo->incCounter(4 * size);
        MPI_Waitall(numRequests, mpi_requests, mpi_stati);
        /* clear buffer */
        delete[] Id_buffer;
        delete[] Tag_buffer;
        delete[] Owner_buffer;
        delete[] Nodes_buffer;
        delete[] send_offset;
        delete[] recv_offset;
        delete[] proc_mask;
        delete[] mpi_requests;
        delete[] mpi_stati;
        delete[] send_count;
        delete[] recv_count;
        delete[] newOwner;
#endif
    } else { // single rank
#pragma omp for private(e,i) schedule(static)
        for (e = 0; e < self->numElements; e++)
        {
            self->Owner[e] = myRank;
            for (i = 0; i < NN; i++)
                self->Nodes[INDEX2(i, e, NN)] = Id[self->Nodes[INDEX2(i, e, NN)]];
        }
    }
}

} // namespace dudley

