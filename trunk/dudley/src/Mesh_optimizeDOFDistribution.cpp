
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "Mesh.h"
#include "IndexList.h"

#include <escript/index.h>

#ifdef ESYS_HAVE_PARMETIS
#include <parmetis.h>
#ifndef REALTYPEWIDTH
typedef float real_t;
#endif
#endif

#include <boost/scoped_array.hpp>

namespace dudley {

#ifdef ESYS_HAVE_PARMETIS
// Checks whether there is any rank which has no vertex. In case 
// such a rank exists, we don't use parmetis since parmetis requires
// that every rank has at least 1 vertex (at line 129 of file
// "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
// any rank has no vertex).
static bool allRanksHaveNodes(escript::JMPI mpiInfo,
                              const IndexVector& distribution)
{
    int ret = 1;

    if (mpiInfo->rank == 0) {
        for (int i = 0; i < mpiInfo->size; i++) {
            if (distribution[i + 1] == distribution[i]) {
                ret = 0;
                break;
            }
        }
        if (ret == 0)
            std::cerr << "INFO: ParMetis is not used since at least one rank "
                         "has no vertex." << std::endl;
    }
    MPI_Bcast(&ret, 1, MPI_INTEGER, 0, mpiInfo->comm);
    return ret==1;
}
#endif

/// optimizes the distribution of DOFs across processors using ParMETIS.
/// On return a new distribution is given and the globalDOF are relabeled
/// accordingly but the mesh has not been redistributed yet
void Mesh::optimizeDOFDistribution(std::vector<index_t>& distribution)
{
    int mpiSize = MPIInfo->size;
    const int myRank = MPIInfo->rank;
    const index_t myFirstVertex = distribution[myRank];
    const index_t myLastVertex = distribution[myRank + 1];
    const dim_t myNumVertices = myLastVertex - myFirstVertex;
    const dim_t numNodes = Nodes->getNumNodes();

    // first step is to distribute the elements according to a global X of DOF
    dim_t len = 0;
    for (int p = 0; p < mpiSize; ++p)
        len = std::max(len, distribution[p + 1] - distribution[p]);

    index_t* partition = new index_t[len];

#ifdef ESYS_HAVE_PARMETIS
    if (mpiSize > 1 && allRanksHaveNodes(MPIInfo, distribution)) {
        boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
        int dim = Nodes->numDim;

        // create the adjacency structure xadj and adjncy
#pragma omp parallel
        {
            // insert contributions from element matrices into columns index
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, Elements,
                    Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, FaceElements,
                    Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, Points,
                    Nodes->globalDegreesOfFreedom);
        }

        // set the coordinates
        real_t* xyz = new real_t[myNumVertices * dim];
#pragma omp parallel for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = Nodes->globalDegreesOfFreedom[i] - myFirstVertex;
            if (k >= 0 && k < myNumVertices) {
                for (int j = 0; j < dim; ++j)
                    xyz[k * dim + j] = (real_t)(Nodes->Coordinates[INDEX2(j, i, dim)]);
            }
        }

        // create the local CSR matrix pattern
        const dim_t globalNumVertices = distribution[mpiSize];
        index_t* ptr = new index_t[myNumVertices + 1];
#pragma omp parallel for
        for (index_t i = 0; i < myNumVertices; ++i) {
            ptr[i] = index_list[i].count(0, globalNumVertices);
        }
        // accumulate ptr
        dim_t s = 0;
        for (index_t i = 0; i < myNumVertices; ++i) {
            const index_t itmp = ptr[i];
            ptr[i] = s;
            s += itmp;
        }
        ptr[myNumVertices] = s;

        // create index
        index_t* index = new index_t[s];
#pragma omp parallel for
        for (index_t i = 0; i < myNumVertices; ++i) {
            index_list[i].toArray(&index[ptr[i]], 0, globalNumVertices, 0);
        }

        index_t wgtflag = 0;
        index_t numflag = 0;
        index_t ncon = 1;
        index_t edgecut;
        index_t impiSize = mpiSize;
        index_t options[3] = { 1, 0, 0 };
        std::vector<real_t> tpwgts(ncon * mpiSize, 1.f / mpiSize);
        std::vector<real_t> ubvec(ncon, 1.05f);
        ParMETIS_V3_PartGeomKway(&distribution[0], ptr, index, NULL, NULL,
                                 &wgtflag, &numflag, &dim, xyz, &ncon,
                                 &impiSize, &tpwgts[0], &ubvec[0], options,
                                 &edgecut, partition, &MPIInfo->comm);
        delete[] xyz;
        delete[] index;
        delete[] ptr;
    } else {
        for (index_t i = 0; i < myNumVertices; ++i)
            partition[i] = 0; // CPU 0 owns all
    }
#else
#pragma omp parallel for
    for (index_t i = 0; i < myNumVertices; ++i)
        partition[i] = myRank;
#endif // ESYS_HAVE_PARMETIS

    // create a new distribution and labeling of the DOF
    IndexVector new_distribution(mpiSize + 1);
#pragma omp parallel
    {
        IndexVector loc_partition_count(mpiSize);
#pragma omp for
        for (index_t i = 0; i < myNumVertices; ++i)
            loc_partition_count[partition[i]]++;
#pragma omp critical
        {
            for (int i = 0; i < mpiSize; ++i)
                new_distribution[i] += loc_partition_count[i];
        }
    }

    IndexVector recvbuf(mpiSize * mpiSize);
#ifdef ESYS_MPI
    // recvbuf will be the concatenation of each CPU's contribution to
    // new_distribution
    MPI_Allgather(&new_distribution[0], mpiSize, MPI_DIM_T, &recvbuf[0],
                  mpiSize, MPI_DIM_T, MPIInfo->comm);
#else
    for (int i = 0; i < mpiSize; ++i)
        recvbuf[i] = new_distribution[i];
#endif
    new_distribution[0] = 0;
    index_t* newGlobalDOFID = new index_t[len];
    for (int rank = 0; rank < mpiSize; rank++) {
        index_t c = 0;
        for (int i = 0; i < myRank; ++i)
            c += recvbuf[rank + mpiSize * i];
        for (index_t i = 0; i < myNumVertices; ++i) {
            if (rank == partition[i]) {
                newGlobalDOFID[i] = new_distribution[rank] + c;
                c++;
            }
        }
        for (int i = myRank + 1; i < mpiSize; ++i)
            c += recvbuf[rank + mpiSize * i];
        new_distribution[rank + 1] = new_distribution[rank] + c;
    }

    // now the overlap needs to be created by sending the partition around
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(myRank + 1);
    int source = MPIInfo->mod_rank(myRank - 1);
#endif
    int current_rank = myRank;
    bool* setNewDOFId = new bool[numNodes];
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i)
        setNewDOFId[i] = true;

    for (int p = 0; p < mpiSize; ++p) {
        const index_t firstVertex = distribution[current_rank];
        const index_t lastVertex = distribution[current_rank + 1];
#pragma omp parallel for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = Nodes->globalDegreesOfFreedom[i];
            if (setNewDOFId[i] && firstVertex <= k && k < lastVertex) {
                Nodes->globalDegreesOfFreedom[i] = newGlobalDOFID[k - firstVertex];
                setNewDOFId[i] = false;
            }
        }

        if (p < mpiSize - 1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(newGlobalDOFID, len, MPI_DIM_T,
                                 dest, MPIInfo->counter(),
                                 source, MPIInfo->counter(),
                                 MPIInfo->comm, &status);
            MPIInfo->incCounter();
#endif
            current_rank = MPIInfo->mod_rank(current_rank - 1);
        }
    }
    for (int i = 0; i < mpiSize + 1; ++i)
        distribution[i] = new_distribution[i];

    delete[] newGlobalDOFID;
    delete[] setNewDOFId;
    delete[] partition;
}

} // namespace dudley

