
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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

template<typename Scalar> struct Coupler;
template<typename T> using Coupler_ptr = boost::shared_ptr<Coupler<T> >;
template<typename T> using const_Coupler_ptr = boost::shared_ptr<const Coupler<T> >;

struct PASO_DLL_API Connector
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


template<typename Scalar>
struct PASO_DLL_API Coupler
{
    Coupler(const_Connector_ptr, dim_t blockSize, escript::JMPI mpiInfo);
    ~Coupler();

    void startCollect(const Scalar* in);
    Scalar* finishCollect();
    void copyAll(Coupler_ptr<Scalar> target) const;
    void fillOverlap(dim_t n, Scalar* x);
    void max(dim_t n, Scalar* x);

    inline const Scalar* borrowLocalData() const { return data; }

    inline const Scalar* borrowRemoteData() const { return recv_buffer; }

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
    Scalar* data;
    Scalar* send_buffer;
    Scalar* recv_buffer;
    MPI_Request* mpi_requests;
    MPI_Status* mpi_stati;
    escript::JMPI mpi_info;
};


} // namespace paso

#endif // __PASO_COUPLER_H__

