
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

#include "EsysMPI.h"
#include "EsysException.h"

#include <vector>

namespace escript
{
  
using DataTypes::dim_t;
using DataTypes::index_t;

JMPI makeInfo(MPI_Comm comm, bool owncom)
{
    if (NoCOMM_WORLD::active() && comm==MPI_COMM_WORLD)
        throw EsysException("Attempt to use the MPI_COMM_WORLD "
                            "communicator when it is blocked.");
    JMPI_* p = new JMPI_(comm, owncom);
    return JMPI(p);
}

JMPI_::JMPI_(MPI_Comm mpicomm, bool owncom)
        : comm(mpicomm), ownscomm(owncom), msg_tag_counter(0)
{
#ifdef ESYS_MPI
    if (mpicomm != MPI_COMM_NULL) {
        if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS ||
                MPI_Comm_size(comm, &size) != MPI_SUCCESS) {
            throw EsysException("JMPI::JMPI: error finding comm rank/size" );
        }
    } else {
        rank = 0;
        size = 0;
    }
#else
    rank = 0;
    size = 1;        
#endif        
}

JMPI_::~JMPI_()
{
#ifdef ESYS_MPI
    if (ownscomm && comm != MPI_COMM_NULL)
        MPI_Comm_free(&comm);
#endif
}

dim_t JMPI_::setDistribution(index_t min_id, index_t max_id,
                             index_t* distribution)
{
    const dim_t N = max_id-min_id+1;
    if (N > 0) {
        const dim_t local_N = N/size;
        const dim_t rest = N-local_N*size;
        for (int p=0; p<size; ++p) {
            if (p < rest) {
                distribution[p]=min_id+(local_N+1)*p;
            } else {
                distribution[p]=min_id+rest+local_N*p;
            }
        }
        distribution[size]=max_id+1;
        if (rest==0) {
            return local_N;
        } else {
            return local_N+1;
        }
    } else {
        for (int p=0; p<size+1; ++p)
            distribution[p]=min_id;
        return 0;
    }
}

void JMPI_::split(dim_t N, dim_t* local_N, index_t* offset) 
{
    *local_N = N/size;
    dim_t rest = N-(*local_N)*size;
    if (rank < rest) {
        (*local_N)++;
        *offset = (*local_N)*rank;
    } else {
        *offset = (*local_N)*rank + rest;
    }
}

// Throw all values in and get the maximum --- used for error checking.
// This used to be implemented as a simple AllReduce.
// However, if there are other (overlapping) communicators in the system,
// they don't react well to getting unexpected/untagged messages.
// To avoid this, we do individual sends to the root which sends the
// result back.
bool checkResult(int res, int& mres, const JMPI& info)
{
    if (info->size==1) {
        mres = res;
        return true;
    }
#ifdef ESYS_MPI
    const int leader = 0;
    const int BIGTAG = getSubWorldTag();
    if (info->rank != leader) {  
        if (MPI_Send(&res, 1, MPI_INT, leader, BIGTAG, info->comm) != MPI_SUCCESS)
            return false;
        MPI_Status status;
        if (MPI_Recv(&mres, 1, MPI_INT, leader, BIGTAG, info->comm, &status) != MPI_SUCCESS)
            return false;
    } else {
        std::vector<MPI_Status> status(info->size - 1);
        MPI_Request* reqs = new MPI_Request[info->size-1];
        int* eres = new int[info->size-1];
        for (int i=0; i<info->size-1; ++i) {
            MPI_Irecv(eres+i, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);
        }
        if (MPI_Waitall(info->size-1, reqs, &status[0]) != MPI_SUCCESS) {
            delete[] reqs;
            delete[] eres;
            return false;
        }
        // now we have them all, find the max
        mres = res;
        for (int i=0; i<info->size-1; ++i) {
            if (mres < eres[i])
                mres = eres[i];
        }
        delete[] eres;
        // now we know what the result should be, send it to the others
        for (int i=0; i<info->size-1; ++i)
            MPI_Isend(&mres, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);

        if (MPI_Waitall(info->size-1, reqs, &status[0]) != MPI_SUCCESS) {
            delete[] reqs;
            return false;
        }
        delete[] reqs;
    }
#endif
    return true;
}


// ensure that the any ranks with an empty src argument end up with the string
// from one of the other ranks.
// without mpi, it makes dest point at a copy of src.
// Expected use case for this code is to ship error messages between ranks.
// As such, it is not written to be speedy
bool shipString(const char* src, char** dest, MPI_Comm& comm)
{
#ifdef ESYS_MPI  
    int rank=0;
    if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS)
        return false; // we have no reason to believe MPI works anymore
    
    int slen = strlen(src);
    // everybody needs to tell everyone if they have a string
    // send your rank if you have a non-empty string else
    // send -1
    int in = (slen ? rank : -1);
    int out;
    if (MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MAX, comm)!=MPI_SUCCESS)
        return false;

    // should not be called under these conditions, but noone had a string
    if (out==-1) {
        *dest = new char[1];
        *dest[0] = '\0';
        return true;
    }
    // since we will be using broadcast, we need to tell everyone how big the
    // string is going to be with an additional bcast
    
    if (MPI_Bcast(&slen, 1, MPI_INT, out, comm) != MPI_SUCCESS)
        return false;

    // now broadcast that string to everyone
    if (rank==out) {
        // I could const_cast src but instead I'll make a copy
        *dest = new char[slen+1];
        strcpy(*dest, src);

        // this guy should just send the string
        if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm) != MPI_SUCCESS)
            return false;

        return true;
    } else {
        *dest = new char[slen+1];
        if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm)!=MPI_SUCCESS)
            return false;

        return true;
    }
#else
    *dest = new char[strlen(src)+1];
    strcpy(*dest, src);
    return true;
#endif
}

namespace 
{
    // true if a split world call is currently running and MPI_COMM_WORLD
    // should not be allowed by default
    bool nocommworldplease=false;
}

NoCOMM_WORLD::NoCOMM_WORLD()
{
    if (nocommworldplease)
        throw EsysException("NoCOMM_WORLD does not nest.");

    nocommworldplease=true;
}

NoCOMM_WORLD::~NoCOMM_WORLD()
{
    nocommworldplease=false;
}  

bool NoCOMM_WORLD::active()
{
    return nocommworldplease;
}

} // namespace escript

