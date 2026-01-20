
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
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

template<typename Scalar>
Coupler<Scalar>::Coupler(const_Connector_ptr conn, dim_t blockSize,
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
        send_buffer = new Scalar[conn->send->numSharedComponents * block_size];
        recv_buffer = new Scalar[conn->recv->numSharedComponents * block_size];
    }
#endif
}

template<typename Scalar>
Coupler<Scalar>::~Coupler()
{
#ifdef ESYS_MPI
    delete[] send_buffer;
    delete[] recv_buffer;
    delete[] mpi_requests;
    delete[] mpi_stati;
#endif
}

template<typename Scalar>
void Coupler<Scalar>::startCollect(const Scalar* in)
{
    data = const_cast<Scalar*>(in);
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        if (in_use) {
            throw PasoException("Coupler::startCollect: Coupler in use.");
        }
        MPI_Datatype mpiType = (sizeof(Scalar) == sizeof(double) ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX);
        // start receiving input
        for (dim_t i=0; i < connector->recv->neighbour.size(); ++i) {
            MPI_Irecv(&recv_buffer[connector->recv->offsetInShared[i]*block_size],
                    (connector->recv->offsetInShared[i+1]-connector->recv->offsetInShared[i])*block_size,
                    mpiType, connector->recv->neighbour[i],
                    mpi_info->counter()+connector->recv->neighbour[i],
                    mpi_info->comm, &mpi_requests[i]);
        }
        // collect values into buffer
        const int numSharedSend = connector->send->numSharedComponents;
        if (block_size > 1) {
            const size_t block_size_size=block_size*sizeof(Scalar);
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
                    mpiType, connector->send->neighbour[i],
                    mpi_info->counter()+mpi_info->rank, mpi_info->comm,
                    &mpi_requests[i+connector->recv->neighbour.size()]);
        }
        mpi_info->incCounter(mpi_info->size);
        in_use = true;
    }
#endif
}

template<typename Scalar>
Scalar* Coupler<Scalar>::finishCollect()
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

template<typename Scalar>
void Coupler<Scalar>::copyAll(Coupler_ptr<Scalar> target) const
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

template<typename Scalar>
void Coupler<Scalar>::fillOverlap(dim_t n, Scalar* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n= n - overlap_n;
    const dim_t offset = block_size * my_n;

    startCollect(x);
    Scalar* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n * block_size; ++i) {
        x[offset+i] = remote_values[i];
    }
}

/* adjusts max values across shared values x */
template<>
void Coupler<real_t>::max(dim_t n, real_t* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n = n - overlap_n;

    startCollect(x);
    real_t* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n; ++i)
        x[my_n+i] = std::max(x[my_n+i], remote_values[i]);
}

template<>
void Coupler<cplx_t>::max(dim_t n, cplx_t* x)
{
    throw PasoException("Coupler::max: invalid call for complex data"); 
}

// instantiate
template struct Coupler<real_t>;
template struct Coupler<cplx_t>;

} // namespace paso

