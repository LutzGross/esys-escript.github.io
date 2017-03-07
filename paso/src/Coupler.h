
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


/****************************************************************************/

/*   Paso: coupler                                            */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_COUPLER_H__
#define __PASO_COUPLER_H__

#include "Paso.h"
#include "PasoException.h"
#include "SharedComponents.h"

namespace paso {

struct Connector;
typedef boost::shared_ptr<Connector> Connector_ptr;
typedef boost::shared_ptr<const Connector> const_Connector_ptr;

struct Coupler;
typedef boost::shared_ptr<Coupler> Coupler_ptr;
typedef boost::shared_ptr<const Coupler> const_Coupler_ptr;

PASO_DLL_API
struct Connector
{
    SharedComponents_ptr send;
    SharedComponents_ptr recv;

    Connector(SharedComponents_ptr s, SharedComponents_ptr r)
    {
        if (s->local_length != r->local_length) {
            throw PasoException("Connector: local length of send and recv "
                                "SharedComponents must match.");
        }
        send = s;
        recv = r;
    }

    /// creates a copy
    inline Connector_ptr copy() const { return unroll(1); }

    inline Connector_ptr unroll(index_t block_size) const
    {
        SharedComponents_ptr new_send_shcomp, new_recv_shcomp;
        Connector_ptr out;
        if (block_size > 1) {
            new_send_shcomp.reset(new SharedComponents(send->local_length,
                        send->neighbour, send->shared, send->offsetInShared,
                        block_size, 0));

            new_recv_shcomp.reset(new SharedComponents(recv->local_length,
                        recv->neighbour, recv->shared, recv->offsetInShared,
                        block_size, 0));
        } else {
            new_send_shcomp = send;
            new_recv_shcomp = recv;
        }
        out.reset(new Connector(new_send_shcomp, new_recv_shcomp));
        return out;
    }

    //inline debug() const
    //{
    //    for (int i=0; i<recv->neighbour.size(); ++i)
    //        printf("Coupler: %d receive %d data at %d from %d\n",
    //            s->mpi_info->rank,recv->offsetInShared[i+1]-recv->offsetInShared[i],
    //            recv->offsetInShared[i],recv->neighbour[i]);
    //    for (int i=0; i<send->neighbour.size(); ++i)
    //        printf("Coupler: %d send %d data at %d to %d\n",
    //            s->mpi_info->rank,send->offsetInShared[i+1]-send->offsetInShared[i],
    //            send->offsetInShared[i],send->neighbour[i]);
    //}
};


PASO_DLL_API
struct Coupler
{
    Coupler(const_Connector_ptr, dim_t blockSize, escript::JMPI mpiInfo);
    ~Coupler();

    void startCollect(const double* in);
    double* finishCollect();
    void copyAll(Coupler_ptr target) const;
    void fillOverlap(dim_t n, double* x);
    void max(dim_t n, double* x);

    inline const double* borrowLocalData() const { return data; }

    inline const double* borrowRemoteData() const { return recv_buffer; }

    inline dim_t getNumSharedComponents() const
    {
        return connector->send->numSharedComponents;
    }

    inline dim_t getNumOverlapComponents() const
    {
        return connector->recv->numSharedComponents;
    }

    inline dim_t getNumSharedValues() const
    {
        return getNumSharedComponents() * block_size;
    }

    inline dim_t getNumOverlapValues() const
    {
        return getNumOverlapComponents() * block_size;
    }

    inline dim_t getLocalLength() const
    {
        return connector->send->local_length;
    }

    const_Connector_ptr connector;
    dim_t block_size;
    bool in_use;

    // unmanaged pointer to data to be sent
    double* data;
    double* send_buffer;
    double* recv_buffer;
    MPI_Request* mpi_requests;
    MPI_Status* mpi_stati;
    escript::JMPI mpi_info;
};


} // namespace paso

#endif // __PASO_COUPLER_H__

