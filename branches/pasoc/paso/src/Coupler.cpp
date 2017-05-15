
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

    
template class Coupler<escript::DataTypes::cplx_t>;
template class Coupler<escript::DataTypes::real_t>;
    
/****************************************************************************
 *
 * allocates a Coupler
 *
 ****************************************************************************/

template<class T>
Coupler<T>::Coupler(const_Connector_ptr conn, dim_t blockSize,
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
        send_buffer = new T[conn->send->numSharedComponents * block_size];
        recv_buffer = new T[conn->recv->numSharedComponents * block_size];
    }
#endif
}

template<class T>
Coupler<T>::~Coupler()
{
#ifdef ESYS_MPI
    delete[] send_buffer;
    delete[] recv_buffer;
    delete[] mpi_requests;
    delete[] mpi_stati;
#endif
}

template<class T>
void Coupler<T>::startCollect(const T* in)
{
    data = const_cast<T*>(in);
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        if (in_use) {
            throw PasoException("Coupler::startCollect: Coupler in use.");
        }
        // start receiving input
        for (dim_t i=0; i < connector->recv->neighbour.size(); ++i) {
            MPI_Irecv(&recv_buffer[connector->recv->offsetInShared[i]*block_size],
                    (connector->recv->offsetInShared[i+1]-connector->recv->offsetInShared[i])*block_size,
                    escript::getMPIType(T(0)), connector->recv->neighbour[i],
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
                    escript::getMPIType(T(0)), connector->send->neighbour[i],
                    mpi_info->counter()+mpi_info->rank, mpi_info->comm,
                    &mpi_requests[i+connector->recv->neighbour.size()]);
        }
        mpi_info->incCounter(mpi_info->size);
        in_use = true;
    }
#endif
}

template<class T>
T* Coupler<T>::finishCollect()
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

template<class T>
void Coupler<T>::copyAll(Coupler_ptr<T> target) const
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

template<class T>
void Coupler<T>::fillOverlap(dim_t n, T* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n= n - overlap_n;
    const dim_t offset = block_size * my_n;

    startCollect(x);
    T* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n * block_size; ++i) {
        x[offset+i] = remote_values[i];
    }
}

/* adjusts max values across shared values x */
template<class T>
void Coupler<T>::max(dim_t n, T* x)
{
    const dim_t overlap_n = getNumOverlapValues();
    const dim_t my_n = n - overlap_n;

    startCollect(x);
    T* remote_values = finishCollect();

#pragma omp parallel for
    for (dim_t i=0; i < overlap_n; ++i)
        x[my_n+i] = std::max(x[my_n+i], remote_values[i]);
}

/* adjusts max values across shared values x */
template<>
void Coupler<escript::DataTypes::cplx_t>::max(dim_t n, escript::DataTypes::cplx_t* x)
{
    throw PasoException("Coupler::max: max operation not defined for complex values.");
}

} // namespace paso

