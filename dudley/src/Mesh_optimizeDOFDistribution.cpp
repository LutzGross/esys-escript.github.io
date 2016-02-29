
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

/****************************************************************************/

/*   Dudley: Mesh: optimizes the distribution of DOFs across processors */
/*   using ParMETIS. On return a new distribution is given and the globalDOF are relabelled */
/*   accordingly but the mesh has not been redistributed yet */

/****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Mesh.h"
#include "IndexList.h"

#ifdef USE_PARMETIS
#include "parmetis.h"
#ifndef REALTYPEWIDTH
typedef float real_t;
#endif
#endif

#include <boost/scoped_array.hpp>

namespace dudley {

/*****************************************************************************
   Check whether there is any node which has no vertex. In case 
   such node exists, we don't use parmetis since parmetis requires
   that every node has at least 1 vertex (at line 129 of file
   "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
   any node has no vertex).
 *****************************************************************************/
#ifdef USE_PARMETIS
int Check_Inputs_For_Parmetis(dim_t mpiSize, dim_t rank, dim_t * distribution, MPI_Comm * comm)
{
    dim_t i, len;
    int ret_val = 1;

    if (rank == 0) {
        for (i = 0; i < mpiSize; i++) {
            len = distribution[i + 1] - distribution[i];
            if (len == 0) {
                ret_val = 0;
                break;
            }
        }
    }
    MPI_Bcast(&ret_val, 1, MPI_INTEGER, 0, *comm);
    if (ret_val == 0)
        printf("INFO: Parmetis is not used since some nodes have no vertex!\n");
    return ret_val;
}
#endif

/*****************************************************************************/

void Dudley_Mesh_optimizeDOFDistribution(Dudley_Mesh* in, dim_t* distribution)
{
    if (in == NULL || in->Nodes == NULL)
        return;

    dim_t i, k;
    int rank;
    int c;

    const int myRank = in->MPIInfo->rank;
    dim_t mpiSize = in->MPIInfo->size;

    // first step is to distribute the elements according to a global X of DOF

    index_t myFirstVertex = distribution[myRank];
    index_t myLastVertex = distribution[myRank + 1];
    dim_t myNumVertices = myLastVertex - myFirstVertex;
    dim_t globalNumVertices = distribution[mpiSize];
    dim_t len = 0;
    for (dim_t p = 0; p < mpiSize; ++p)
        len = std::max(len, distribution[p + 1] - distribution[p]);

    index_t* partition = new index_t[len];
    dim_t* partition_count = new dim_t[mpiSize + 1];
    dim_t* new_distribution = new dim_t[mpiSize + 1];
    index_t* newGlobalDOFID = new index_t[len];
    bool* setNewDOFId = new bool[in->Nodes->numNodes];
    dim_t* recvbuf = new dim_t[mpiSize * mpiSize];

#ifdef USE_PARMETIS
    dim_t dim = in->Nodes->numDim;
    real_t* xyz = new real_t[myNumVertices * dim];

    /* set the coordinates: */
    /* it is assumed that at least one node on this processor provides a coordinate */
#pragma omp parallel for private(i,k)
    for (i = 0; i < in->Nodes->numNodes; ++i)
    {
        k = in->Nodes->globalDegreesOfFreedom[i] - myFirstVertex;
        if ((k >= 0) && (k < myNumVertices))
        {
            for (dim_t j = 0; j < dim; ++j)
                xyz[k * dim + j] = (real_t)(in->Nodes->Coordinates[INDEX2(j, i, dim)]);
        }
    }
#endif // USE_PARMETIS

    boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
    /* ksteube CSR of DOF IDs */
    /* create the adjacency structure xadj and adjncy */
#pragma omp parallel
    {
        /* ksteube build CSR format */
        /*  insert contributions from element matrices into columns index index_list: */
        Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
        myFirstVertex, myLastVertex, in->Elements,
        in->Nodes->globalDegreesOfFreedom,
        in->Nodes->globalDegreesOfFreedom);
        Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
        myFirstVertex, myLastVertex, in->FaceElements,
        in->Nodes->globalDegreesOfFreedom,
        in->Nodes->globalDegreesOfFreedom);
        Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
        myFirstVertex, myLastVertex, in->Points,
        in->Nodes->globalDegreesOfFreedom,
        in->Nodes->globalDegreesOfFreedom);
    }

    /* create the local matrix pattern */
    paso::Pattern_ptr pattern = paso::Pattern::fromIndexListArray(0,
            myNumVertices, index_list.get(), 0, globalNumVertices, 0);

#ifdef USE_PARMETIS
    if (mpiSize > 1 && Check_Inputs_For_Parmetis(mpiSize, myRank, distribution, &(in->MPIInfo->comm)) > 0)
    {
        int i;
        int wgtflag = 0;
        int numflag = 0;    /* pattern->ptr is C style: starting from 0 instead of 1 */
        int ncon = 1;
        int edgecut;
        int options[3];
        real_t *tpwgts = new real_t[ncon * mpiSize];
        real_t *ubvec = new real_t[ncon];
        for (i = 0; i < ncon * mpiSize; i++)
            tpwgts[i] = 1.0 / (real_t)mpiSize;
        for (i = 0; i < ncon; i++)
            ubvec[i] = 1.05;
        options[0] = 1;
        options[1] = 15;
        options[2] = 0;
        ParMETIS_V3_PartGeomKway(distribution, pattern->ptr,
                pattern->index, NULL, NULL, &wgtflag, &numflag, &dim,
                xyz, &ncon, &mpiSize, tpwgts, ubvec, options, &edgecut,
                partition, /* new CPU ownership of elements */
                &in->MPIInfo->comm);
        //printf("ParMETIS number of edges cut by partitioning per processor: %d\n", edgecut/std::max(in->MPIInfo->size,1));
        delete[] xyz;
        delete[] ubvec;
        delete[] tpwgts;
    }
    else
    {
        for (i = 0; i < myNumVertices; ++i)
            partition[i] = 0;       /* CPU 0 owns it */
    }
#else
    for (i = 0; i < myNumVertices; ++i)
        partition[i] = myRank;      /* CPU 0 owns it */
#endif // USE_PARMETIS

    // create a new distribution and labelling of the DOF
    const size_t mpiSize_size = mpiSize * sizeof(dim_t);
    memset(new_distribution, 0, mpiSize_size);
#pragma omp parallel
    {
        dim_t* loc_partition_count = new dim_t[mpiSize];
        memset(loc_partition_count, 0, mpiSize_size);
#pragma omp for private(i)
        for (i = 0; i < myNumVertices; ++i)
            loc_partition_count[partition[i]]++;
#pragma omp critical
        {
            for (i = 0; i < mpiSize; ++i)
                new_distribution[i] += loc_partition_count[i];
        }
        delete[] loc_partition_count;
    }
#ifdef ESYS_MPI
    // recvbuf will be the concatenation of each CPU's contribution to
    // new_distribution
    MPI_Allgather(new_distribution, mpiSize, MPI_INT, recvbuf, mpiSize, MPI_INT, in->MPIInfo->comm);
#else
    for (i = 0; i < mpiSize; ++i)
        recvbuf[i] = new_distribution[i];
#endif
    new_distribution[0] = 0;
    for (rank = 0; rank < mpiSize; rank++)
    {
        c = 0;
        for (i = 0; i < myRank; ++i)
            c += recvbuf[rank + mpiSize * i];
        for (i = 0; i < myNumVertices; ++i)
        {
            if (rank == partition[i])
            {
                newGlobalDOFID[i] = new_distribution[rank] + c;
                c++;
            }
        }
        for (i = myRank + 1; i < mpiSize; ++i)
            c += recvbuf[rank + mpiSize * i];
        new_distribution[rank + 1] = new_distribution[rank] + c;
    }
    delete[] recvbuf;

    // now the overlap needs to be created by sending the partition around
#ifdef ESYS_MPI
    int dest = in->MPIInfo->mod_rank(myRank + 1);
    int source = in->MPIInfo->mod_rank(myRank - 1);
#endif
    int current_rank = myRank;
#pragma omp parallel for private(i)
    for (i = 0; i < in->Nodes->numNodes; ++i)
        setNewDOFId[i] = true;

    for (dim_t p = 0; p < mpiSize; ++p)
    {
        index_t firstVertex = distribution[current_rank];
        index_t lastVertex = distribution[current_rank + 1];
#pragma omp parallel for private(i,k)
        for (i = 0; i < in->Nodes->numNodes; ++i)
        {
            k = in->Nodes->globalDegreesOfFreedom[i];
            if (setNewDOFId[i] && (firstVertex <= k) && (k < lastVertex))
            {
                in->Nodes->globalDegreesOfFreedom[i] = newGlobalDOFID[k - firstVertex];
                setNewDOFId[i] = false;
            }
        }

        if (p < mpiSize - 1)
        {               /* the final send can be skipped */
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(newGlobalDOFID, len, MPI_INT,
                                 dest, in->MPIInfo->counter(),
                                 source, in->MPIInfo->counter(), in->MPIInfo->comm, &status);
            in->MPIInfo->incCounter();
#endif
            current_rank = in->MPIInfo->mod_rank(current_rank - 1);
        }
    }
    for (i = 0; i < mpiSize + 1; ++i)
        distribution[i] = new_distribution[i];

    delete[] newGlobalDOFID;
    delete[] setNewDOFId;
    delete[] new_distribution;
    delete[] partition_count;
    delete[] partition;
}

} // namespace dudley

