
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#include "Coupler.h"
#include "esysUtils/error.h"

namespace paso {

/****************************************************************************
 *
 * allocates a Coupler
 *
 ****************************************************************************/

Coupler::Coupler(Connector_ptr conn, dim_t blockSize) :
    connector(conn),
    block_size(blockSize),
    in_use(false),
    data(NULL),
    send_buffer(NULL),
    recv_buffer(NULL),
    mpi_requests(NULL),
    mpi_stati(NULL)
{
    Esys_resetError();
    mpi_info = Esys_MPIInfo_getReference(conn->mpi_info);
#ifdef ESYS_MPI
    mpi_requests = new MPI_Request[conn->send->numNeighbors +
                                   conn->recv->numNeighbors];
    mpi_stati = new MPI_Status[conn->send->numNeighbors +
                               conn->recv->numNeighbors];
#endif
    if (mpi_info->size > 1) {
        send_buffer=new double[conn->send->numSharedComponents * block_size];
        recv_buffer=new double[conn->recv->numSharedComponents * block_size];
    }
}

Coupler::~Coupler()
{
    delete[] send_buffer;
    delete[] recv_buffer;
#ifdef ESYS_MPI
    delete[] mpi_requests;
    delete[] mpi_stati;
#endif                
    Esys_MPIInfo_free(mpi_info);
}

void Coupler::startCollect(const double* in)
{
    data = const_cast<double*>(in);
    if (mpi_info->size > 1) {
        if (in_use) {
            Esys_setError(SYSTEM_ERROR,"Coupler::startCollect: Coupler in use.");
        }
        // start receiving input
        for (dim_t i=0; i < connector->recv->numNeighbors; ++i) {
#ifdef ESYS_MPI
            MPI_Irecv(&recv_buffer[connector->recv->offsetInShared[i]*block_size],
                    (connector->recv->offsetInShared[i+1]-connector->recv->offsetInShared[i])*block_size,
                    MPI_DOUBLE, connector->recv->neighbor[i], 
                    mpi_info->msg_tag_counter+connector->recv->neighbor[i],
                    mpi_info->comm, &mpi_requests[i]);
#endif
        }
        // collect values into buffer
        if (block_size > 1) {
            const size_t block_size_size=block_size*sizeof(double);
#pragma omp parallel for
            for (dim_t i=0; i < connector->send->numSharedComponents; ++i) {
                memcpy(&(send_buffer[(block_size)*i]),
                       &(in[block_size*connector->send->shared[i]]),
                       block_size_size);
            }
        } else {
#pragma omp parallel for
            for (dim_t i=0; i < connector->send->numSharedComponents; ++i) {
                send_buffer[i]=in[connector->send->shared[i]];
            }
        }
        // send buffer out
        for (dim_t i=0; i < connector->send->numNeighbors; ++i) {
#ifdef ESYS_MPI
            MPI_Issend(&send_buffer[connector->send->offsetInShared[i]*block_size],
                    (connector->send->offsetInShared[i+1] - connector->send->offsetInShared[i])*block_size,
                    MPI_DOUBLE, connector->send->neighbor[i], 
                    mpi_info->msg_tag_counter+mpi_info->rank, mpi_info->comm,
                    &mpi_requests[i+connector->recv->numNeighbors]);
#endif 
        }
        ESYS_MPI_INC_COUNTER(*mpi_info, mpi_info->size)
        in_use = true;
    }
}

double* Coupler::finishCollect()
{
    if (mpi_info->size > 1) {
        if (!in_use) {
            Esys_setError(SYSTEM_ERROR, "Coupler::finishCollect: Communication has not been initiated.");
            return NULL;
        }
        // wait for receive
#ifdef ESYS_MPI
        MPI_Waitall(connector->recv->numNeighbors+connector->send->numNeighbors,
                    mpi_requests, mpi_stati);
#endif
        in_use = false;
    }

    return recv_buffer;
}

void Coupler::copyAll(Coupler_ptr target) const
{
#pragma omp parallel 
    {
#pragma omp for
        for (dim_t i=0; i < getNumOverlapValues(); ++i) {
            target->recv_buffer[i] = recv_buffer[i];
        }
#pragma omp for
        for (dim_t i=0; i < getLocalLength()*block_size; ++i) {
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
        x[my_n+i] = MAX(x[my_n+i], remote_values[i]);
}

} // namespace paso

