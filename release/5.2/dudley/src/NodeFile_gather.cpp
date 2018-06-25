
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "NodeFile.h"

#include <escript/index.h>

using escript::DataTypes::real_t;

namespace dudley {

// helper function
static void gatherEntries(dim_t n, const index_t* index,
                          index_t min_index, index_t max_index,
                          index_t* Id_out, const index_t* Id_in,
                          int* Tag_out, const int* Tag_in,
                          index_t* globalDegreesOfFreedom_out,
                          const index_t* globalDegreesOfFreedom_in,
                          int numDim, real_t* Coordinates_out,
                          const real_t* Coordinates_in)
{
    const dim_t range = max_index - min_index;
    const size_t numDim_size = numDim * sizeof(real_t);
#pragma omp parallel for
    for (index_t i = 0; i < n; i++) {
        const index_t k = index[i] - min_index;
        if (k >= 0 && k < range) {
            Id_out[i] = Id_in[k];
            Tag_out[i] = Tag_in[k];
            globalDegreesOfFreedom_out[i] = globalDegreesOfFreedom_in[k];
            memcpy(&Coordinates_out[INDEX2(0, i, numDim)],
                   &Coordinates_in[INDEX2(0, k, numDim)], numDim_size);
        }
    }
}

// helper function
static void scatterEntries(dim_t n, const index_t* index,
                           index_t min_index, index_t max_index,
                           index_t* Id_out, const index_t* Id_in,
                           int* Tag_out, const int* Tag_in,
                           index_t* globalDegreesOfFreedom_out,
                           const index_t* globalDegreesOfFreedom_in,
                           int numDim, real_t* Coordinates_out,
                           const real_t* Coordinates_in)
{
    const dim_t range = max_index - min_index;
    const size_t numDim_size = numDim * sizeof(real_t);

#pragma omp parallel for
    for (index_t i = 0; i < n; i++) {
        const index_t k = index[i] - min_index;
        if (k >= 0 && k < range) {
            Id_out[k] = Id_in[i];
            Tag_out[k] = Tag_in[i];
            globalDegreesOfFreedom_out[k] = globalDegreesOfFreedom_in[i];
            memcpy(&Coordinates_out[INDEX2(0, k, numDim)],
                   &Coordinates_in[INDEX2(0, i, numDim)], numDim_size);
        }
    }
}

void NodeFile::gather(const index_t* index, const NodeFile* in)
{
    const std::pair<index_t,index_t> idRange(in->getGlobalIdRange());
    gatherEntries(numNodes, index, idRange.first, idRange.second, Id, in->Id,
             Tag, in->Tag, globalDegreesOfFreedom, in->globalDegreesOfFreedom,
             numDim, Coordinates, in->Coordinates);
}

void NodeFile::gather_global(const index_t* index, const NodeFile* in)
{
    // get the global range of node IDs
    const std::pair<index_t,index_t> idRange(in->getGlobalIdRange());
    const index_t UNDEFINED = idRange.first - 1;
    std::vector<index_t> distribution(in->MPIInfo->size + 1);

    // distribute the range of node IDs
    dim_t buffer_len = MPIInfo->setDistribution(idRange.first, idRange.second,
                                                &distribution[0]);

    // allocate buffers
    index_t* Id_buffer = new index_t[buffer_len];
    int* Tag_buffer = new int[buffer_len];
    index_t* globalDegreesOfFreedom_buffer = new index_t[buffer_len];
    real_t* Coordinates_buffer = new real_t[buffer_len * numDim];

    // fill Id_buffer by the UNDEFINED marker to check if nodes are
    // defined
#pragma omp parallel for
    for (index_t n = 0; n < buffer_len; n++)
        Id_buffer[n] = UNDEFINED;

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
#ifdef ESYS_MPI
        if (p > 0) { // the initial send can be skipped
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_DIM_T, dest,
                        MPIInfo->counter(), source, MPIInfo->counter(),
                        MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                        MPIInfo->counter() + 1, source,
                        MPIInfo->counter() + 1, MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                        MPI_DIM_T, dest, MPIInfo->counter() + 2, source,
                        MPIInfo->counter() + 2, MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len * numDim,
                        MPI_DOUBLE, dest, MPIInfo->counter() + 3, source,
                        MPIInfo->counter() + 3, MPIInfo->comm, &status);
            MPIInfo->incCounter(4);
        }
#endif
        buffer_rank = MPIInfo->mod_rank(buffer_rank - 1);
        scatterEntries(in->numNodes, in->Id, distribution[buffer_rank],
                       distribution[buffer_rank + 1], Id_buffer, in->Id,
                       Tag_buffer, in->Tag, globalDegreesOfFreedom_buffer,
                       in->globalDegreesOfFreedom, numDim, Coordinates_buffer,
                       in->Coordinates);
    }
    // now entries are collected from the buffer again by sending the entries
    // around in a circle
#ifdef ESYS_MPI
    dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
        gatherEntries(numNodes, index, distribution[buffer_rank],
                      distribution[buffer_rank + 1], Id, Id_buffer,
                      Tag, Tag_buffer, globalDegreesOfFreedom,
                      globalDegreesOfFreedom_buffer, numDim,
                      Coordinates, Coordinates_buffer);
#ifdef ESYS_MPI
        if (p < MPIInfo->size - 1) { // the last send can be skipped
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_DIM_T, dest,
                              MPIInfo->counter(), source,
                              MPIInfo->counter(), MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                              MPIInfo->counter() + 1, source,
                              MPIInfo->counter() + 1, MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                              MPI_DIM_T, dest, MPIInfo->counter() + 2, source,
                              MPIInfo->counter() + 2, MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len * numDim,
                              MPI_DOUBLE, dest, MPIInfo->counter() + 3, source,
                              MPIInfo->counter() + 3, MPIInfo->comm, &status);
            MPIInfo->incCounter(4);
        }
#endif
        buffer_rank = MPIInfo->mod_rank(buffer_rank - 1);
    }
    delete[] Id_buffer;
    delete[] Tag_buffer;
    delete[] globalDegreesOfFreedom_buffer;
    delete[] Coordinates_buffer;
#if DOASSERT
    // check if all nodes are set
    index_t err = -1;
#pragma omp parallel for
    for (index_t n = 0; n < numNodes; ++n) {
        if (Id[n] == UNDEFINED) {
#pragma omp critical
            err = n;
        }
    }
    if (err >= 0) {
        std::stringstream ss;
        ss << "NodeFile::gather_global: Node id " << Id[err]
            << " at position " << err << " is referenced but not defined.";
        throw escript::AssertException(ss.str());
    }
#endif // DOASSERT
}

} // namespace dudley

