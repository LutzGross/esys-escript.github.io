
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

#include "Esys_MPI.h"
#include "error.h"
#include "index.h"
#include "EsysException.h"

#include <vector>


namespace esysUtils
{
  
JMPI makeInfo(MPI_Comm comm, bool owncom)
{
    if (esysUtils::NoCOMM_WORLD::active() && comm==MPI_COMM_WORLD)
    {
        throw EsysException("Attempt to use the MPI_COMM_WORLD communicator when it is blocked.");
    }
    JMPI_* p=new JMPI_(comm, owncom);
    return JMPI(p);
}


JMPI_::JMPI_(MPI_Comm mpicomm, bool owncom)
        : comm(mpicomm), ownscomm(owncom)
{
    msg_tag_counter = 0;
#ifdef ESYS_MPI
    if (mpicomm!=MPI_COMM_NULL)
    {
        if (MPI_Comm_rank(comm, &rank)!=MPI_SUCCESS || MPI_Comm_size(comm, &size)!=MPI_SUCCESS)
        {
            throw EsysException("JMPI::JMPI: error finding comm rank/size" );
        }
    }
    else
    {
        rank=0;
        size=0;
    }
#else
    rank=0;
    size=1;        
#endif        
}

JMPI_::~JMPI_()
{
#ifdef ESYS_MPI
    if (ownscomm && (comm!=MPI_COMM_NULL))
    {
        MPI_Comm_free(&comm);
    }
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

dim_t Esys_MPIInfo_setDistribution(JMPI& mpi_info, index_t min_id,
                                   index_t max_id,index_t* distribution)
{
    const int s = mpi_info->size;
    const dim_t N = max_id-min_id+1;
    if (N > 0) {
        const dim_t local_N = N/s;
        const dim_t rest = N-local_N*s;
        for (int p=0; p<s; ++p) {
            if (p < rest) {
                distribution[p]=min_id+(local_N+1)*p;
            } else {
                distribution[p]=min_id+rest+local_N*p;
            }
        }
        distribution[s]=max_id+1;
        if (rest==0) {
            return local_N;
        } else {
            return local_N+1;
        }
    } else {
        for (int p=0; p<s+1; ++p)
            distribution[p]=min_id;
        return 0;
    }
}


// N = #CPUs, k is a CPU number but out of range or even negative.
// Return a CPU number in 0...n-1.
index_t mod_rank(index_t n, index_t k) 
{
    index_t out=0;
    if (n > 1) {
        const index_t q = k/n;
        if (k > 0) {
           out=k-n*q;
        } else if (k < 0) {
           out=k-n*(q-1);
        }
    }
    return out;
}


/* checks that there is no error across all processes in a communicator */
/* NOTE: does not guarantee consistency of error string on each process */
bool Esys_MPIInfo_noError(const JMPI& mpi_info)
{
    int errorLocal = Esys_noError() ? 0 : 1;
    int errorGlobal = errorLocal;

#ifdef ESYS_MPI
    if (!checkResult(errorLocal, errorGlobal, mpi_info))
    {
        return false;
    }
    if (errorLocal==0 && errorGlobal==1) 
    {
        Esys_setError(ESYS_MPI_ERROR, "Esys_MPIInfo_noError(): there was an error on another MPI process" );
    }
#endif

    return (errorGlobal==0);
}

// Throw all values in and get the maximum --- used for error checking.
// This used to be implemented as a simple AllReduce.
// However, if there are other (overlapping) communicators in the system, they don't
// react well to getting unexpected/untagged messages.
// To avoid this, we do individual sends to the root which sends the result back.
bool checkResult(int res, int& mres, const esysUtils::JMPI& info)
{
    if (info->size==1)
    {
        mres=res;
        return true;
    }
#ifdef ESYS_MPI
    const int leader=0;
    const int BIGTAG=esysUtils::getSubWorldTag();
    if (info->rank!=leader)
    {  
        if (MPI_Send(&res, 1, MPI_INT, leader, BIGTAG, info->comm)!=MPI_SUCCESS)
            return false;
        MPI_Status status;
        if (MPI_Recv(&mres, 1, MPI_INT, leader, BIGTAG, info->comm, &status)!=MPI_SUCCESS)
            return false;
    }
    else
    {
        std::vector<MPI_Status> status(info->size - 1);
        MPI_Request* reqs=new MPI_Request[info->size-1];
        int* eres=new int[info->size-1];
        for (int i=0;i<info->size-1;++i)
        {
            MPI_Irecv(eres+i, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);          
        }  
        if (MPI_Waitall(info->size-1, reqs, &status[0])!=MPI_SUCCESS)
        {
            delete[] reqs;
            delete[] eres;
            return false;
        }
        // now we have them all, find the max
        mres=res;
        for (int i=0;i<info->size-1;++i)
        {
            if (mres<eres[i])
            {
                mres=eres[i];
            }
        }
        delete[] eres;
        // now we know what the result should be
        // send it to the others
        for (int i=0;i<info->size-1;++i)
        {
            MPI_Isend(&mres, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);          
        }
        if (MPI_Waitall(info->size-1, reqs, &status[0])!=MPI_SUCCESS)
        {
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
    Esys_MPI_rank rank=0;
    if (MPI_Comm_rank( comm, &rank )!=MPI_SUCCESS)
    {
        return false;        // we have no reason to believe MPI works anymore
    }
    
    int slen=strlen(src);
    // everybody needs to tell everyone if they have a string
    // send your rank if you have a non-empty string else
    // send -1
    int in=(slen?rank:-1);
    int out;
    if (MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MAX, comm)!=MPI_SUCCESS)
    {
        return false;
    }
    if (out==-1)                // should not be called under these conditions, but noone had a string
    {
        *dest=new char[1];
        *dest[0]='\0';
        return true;
    }
    // since we will be using broadcast, we need to tell everyone how big the string is going to be
    // with an additional bcast
    
    if (MPI_Bcast(&slen, 1, MPI_INT, out, comm)!=MPI_SUCCESS)
    {
        return false;
    }
    // now broadcast that string to everyone
    if (rank==out)
    {
        // I could const _ cast src but instead I'll make a copy
        
        *dest=new char[slen+1];
        strcpy(*dest, src);
        
        // this guy should just send the string
        if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm)!=MPI_SUCCESS)
        {
            return false;
        }
        return true;
    }
    else
    {
        *dest=new char[slen+1];
        if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm)!=MPI_SUCCESS)
        {
            return false;
        }
        return true;
    }
#else
    *dest=new char[strlen(src)+1];
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
    {
        throw EsysException("NoCOMM_WORLD does not nest.");
    }
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

/**************************************************
                 WRAPPERS 
**************************************************/

#ifndef _OPENMP 
int serial_get_max_threads(void) {
   return 1;
}
int serial_get_thread_num(void) {
   return 0;
}
#endif

} // namespace esysUtils


