
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

#include "Coupler.h"

#include <cstring> // memcpy

namespace paso {

/****************************************************************************
 *
 * allocates a Coupler
 *
 ****************************************************************************/

Coupler::Coupler(const_Connector_ptr conn, dim_t blockSize,
                 escript::JMPI mpiInfo) :
    connector(conn),
    block_size(blockSize),
    in_use(false),
    data(NULL),
    send_buffer(NULL),
    recv_buffer(NULL),
    mpi_requests(NULL),
    mpi_stati(NULL),
    mpi_info(mpiInfo)
{
#ifdef ESYS_MPI
    mpi_requests = new MPI_Request[conn->send->neighbour.size() +
                                   conn->recv->neighbour.size()];
    mpi_stati = new MPI_Status[conn->send->neighbour.size() +
                               conn->recv->neighbour.size()];
    if (mpi_info->size > 1) {
        send_buffer = new double[conn->send->numSharedComponents * block_size];
        recv_buffer = new double[conn->recv->numSharedComponents * block_size];
    }
#endif
}

Coupler::~Coupler()
{
#ifdef ESYS_MPI
    delete[] send_buffer;
    delete[] recv_buffer;
    delete[] mpi_requests;
    delete[] mpi_stati;
#endif
}

void Coupler::startCollect(const double* in)
{
    data = const_cast<double*>(in);
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        if (in_use) {
            throw PasoException("Coupler::startCollect: Coupler in use.");
        }
        // start receiving input
        for (dim_t i=0; i < connector->recv->neighbour.size(); ++i) {
            MPI_Irecv(&recv_buffer[connector->recv->offsetInShared[i]*block_size],
                    (connector->recv->offsetInShared[i+1]-connector->recv->offsetInShared[i])*block_size,
                    MPI_DOUBLE, connector->recv->neighbour[i],
                    mpi_info->counter()+connector->recv->neighbour[i],
                    mpi_info->comm, &mpi_requests[i]);
        }
        // collect values into buffer
        const int numSharedSend = connector->send->numSharedComponents;
        if (block_size > 1) {
            const size_t block_size_size=block_size*sizeof(double);
#pragma omp parallel for
            for (dim_t i=0; i < numSharedSend; ++i) {
                memcpy(&(send_buffer[(block_size)*i]),
                       &(in[block_size*connector->send->shared[i]]),
                       block_size_size);
            }
        } else {
#pragma omp parallel for
            for (dim_t i=0; i < numSharedSend; ++i) {
                send_buffer[i]=in[connector->send->shared[i]];
            }
        }
        // send buffer out
        for (dim_t i=0; i < connector->send->neighbour.size(); ++i) {
            MPI_Issend(&send_buffer[connector->send->offsetInShared[i]*block_size],
                    (connector->send->offsetInShared[i+1] - connector->send->offsetInShared[i])*block_size,
                    MPI_DOUBLE, connector->send->neighbour[i],
                    mpi_info->counter()+mpi_info->rank, mpi_info->comm,
                    &mpi_requests[i+connector->recv->neighbour.size()]);
        }
        mpi_info->incCounter(mpi_info->size);
        in_use = true;
    }
#endif
}

double* Coupler::finishCollect()
{
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        if (!in_use) {
            throw PasoException("Coupler::finishCollect: Communication has not been initiated.");
        }
        // wait for receive
        MPI_Waitall(connector->recv->neighbour.size() +
                    connector->send->neighbour.size(), mpi_requests, mpi_stati);
        in_use = false;
    }
#endif
    return recv_buffer;
}

void Coupler::copyAll(Coupler_ptr target) const
{
    const dim_t overlap = getNumOverlapValues();
    const dim_t localSize = getLocalLength()*block_size;
#pragma omp parallel
    {
#pragma omp for
        for (dim_t i=0; i < overlap; ++i) {
            target->recv_buffer[i] = recv_buffer[i];
        }
#pragma omp for
        for (dim_t i=0; i < localSize; ++i) {
            target->data[i] = data[i];
        }
    }
}

void Coupler::fillOverlap(dim_t n, double* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n= n - overlap_n;
    const dim_t offset = block_size * my_n;

    startCollect(x);
    double* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n * block_size; ++i) {
        x[offset+i] = remote_values[i];
    }
}

/* adjusts max values across shared values x */
void Coupler::max(dim_t n, double* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n = n - overlap_n;

    startCollect(x);
    double* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n; ++i)
        x[my_n+i] = std::max(x[my_n+i], remote_values[i]);
}

} // namespace paso

