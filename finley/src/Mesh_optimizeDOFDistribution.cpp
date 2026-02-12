
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

#include "FinleyDomain.h"
#include "IndexList.h"

#include <escript/index.h>

#ifdef ESYS_HAVE_PARMETIS
#include <parmetis.h>
#ifndef REALTYPEWIDTH
typedef float real_t;
#endif
#endif

#ifdef ESYS_HAVE_SCOTCH
#include <ptscotch.h>
#ifndef SCOTCH_NUMIDX
typedef SCOTCH_Num scotch_idx_t;
#else
typedef SCOTCH_NUMIDX scotch_idx_t;
#endif
// Define real_t for coordinate data if not already defined
#ifndef REALTYPEWIDTH
typedef float real_t;
#endif
#endif

#include <iostream>
#include <boost/scoped_array.hpp>

namespace finley {

#if defined(ESYS_HAVE_PARMETIS) || defined(ESYS_HAVE_SCOTCH)
// Checks whether there is any rank which has no vertex. In case
// such a rank exists, we don't use graph partitioning since both
// ParMETIS and PT-Scotch require that every rank has at least 1 vertex.
static bool allRanksHaveNodes(escript::JMPI mpiInfo,
                              const IndexVector& distribution,
                              const char* libname)
{
    int ret = 1;

    if (mpiInfo->rank == 0) {
        for (int i = 0; i < mpiInfo->size; i++) {
            if (distribution[i + 1] == distribution[i]) {
                ret = 0;
                break;
            }
        }
        if (ret == 0) {
            std::cerr << "INFO: " << libname << " is not used since at least one rank "
                         "has no vertex." << std::endl;
        }
    }
    MPI_Bcast(&ret, 1, MPI_INTEGER, 0, mpiInfo->comm);
    return ret==1;
}
#endif

/// Optimizes the distribution of DOFs across processors using graph partitioning.
/// Uses ParMETIS if available, otherwise PT-Scotch if available, otherwise no redistribution.
/// On return a new distribution is given and the globalDOF are relabeled
/// accordingly but the mesh has not been redistributed yet.
void FinleyDomain::optimizeDOFDistribution(IndexVector& distribution)
{
    int mpiSize = m_mpiInfo->size;
    const int myRank = m_mpiInfo->rank;
    const index_t myFirstVertex = distribution[myRank];
    const index_t myLastVertex = distribution[myRank + 1];
    const dim_t myNumVertices = myLastVertex - myFirstVertex;
    const dim_t numNodes = m_nodes->getNumNodes();

    // first step is to distribute the elements according to a global X of DOF
    dim_t len = 0;
    for (int p = 0; p < mpiSize; ++p)
        len = std::max(len, distribution[p + 1] - distribution[p]);

    index_t* partition = new index_t[len];

#if defined(ESYS_HAVE_PARMETIS) || defined(ESYS_HAVE_SCOTCH)
    bool useGraphPartitioner = false;
#ifdef ESYS_HAVE_PARMETIS
    useGraphPartitioner = (mpiSize > 1 && allRanksHaveNodes(m_mpiInfo, distribution, "ParMETIS"));
#elif defined(ESYS_HAVE_SCOTCH)
    useGraphPartitioner = (mpiSize > 1 && allRanksHaveNodes(m_mpiInfo, distribution, "PT-Scotch"));
#endif

    if (useGraphPartitioner) {
        boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
        int dim = m_nodes->numDim;

        // create the adjacency structure xadj and adjncy
#pragma omp parallel
        {
            // insert contributions from element matrices into columns index
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, m_elements,
                    m_nodes->globalDegreesOfFreedom, m_nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, m_faceElements,
                    m_nodes->globalDegreesOfFreedom, m_nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, m_contactElements,
                    m_nodes->globalDegreesOfFreedom, m_nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, m_points,
                    m_nodes->globalDegreesOfFreedom, m_nodes->globalDegreesOfFreedom);
        }

        // set the coordinates
        real_t* xyz = new real_t[myNumVertices * dim];
#pragma omp parallel for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = m_nodes->globalDegreesOfFreedom[i] - myFirstVertex;
            if (k >= 0 && k < myNumVertices) {
                for (int j = 0; j < dim; ++j)
                    xyz[k * dim + j] = static_cast<real_t>(m_nodes->Coordinates[INDEX2(j, i, dim)]);
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

        // Partition using available library
#ifdef ESYS_HAVE_PARMETIS
        // Use ParMETIS geometric partitioning
        index_t wgtflag = 0;
        index_t numflag = 0;
        index_t ncon = 1;
        index_t edgecut;
        index_t impiSize = mpiSize;
        index_t idim = dim;
        // options[0]=1 -> non-default values, evaluate rest of options
        // options[1]=0 -> debug level (no output)
        // options[2] -> random seed
        index_t options[3] = { 1, 0, 0 };
        std::vector<real_t> tpwgts(ncon * mpiSize, 1.f / mpiSize);
        std::vector<real_t> ubvec(ncon, 1.05f);
        ParMETIS_V3_PartGeomKway(&distribution[0], ptr, index, NULL, NULL,
                                 &wgtflag, &numflag, &idim, xyz, &ncon,
                                 &impiSize, &tpwgts[0], &ubvec[0], options,
                                 &edgecut, partition, &m_mpiInfo->comm);
#elif defined(ESYS_HAVE_SCOTCH)
        // Use PT-Scotch distributed graph partitioning
        SCOTCH_Dgraph grafdat;
        SCOTCH_dgraphInit(&grafdat, m_mpiInfo->comm);

        // Convert index_t arrays to SCOTCH_Num arrays
        std::vector<SCOTCH_Num> scotch_vertloctab(myNumVertices + 1);
        std::vector<SCOTCH_Num> scotch_edgeloctab(s);
        std::vector<SCOTCH_Num> scotch_parttab(myNumVertices);

        for (index_t i = 0; i <= myNumVertices; ++i)
            scotch_vertloctab[i] = static_cast<SCOTCH_Num>(ptr[i]);
        for (index_t i = 0; i < s; ++i)
            scotch_edgeloctab[i] = static_cast<SCOTCH_Num>(index[i]);

        // Build distributed graph
        SCOTCH_Num baseval = 0;  // 0-based indexing
        SCOTCH_dgraphBuild(&grafdat,
                          baseval,                                    // baseval
                          static_cast<SCOTCH_Num>(myNumVertices),    // vertlocnbr
                          static_cast<SCOTCH_Num>(myNumVertices),    // vertlocmax
                          &scotch_vertloctab[0],                     // vertloctab
                          NULL,                                       // vendloctab (use vertloctab[i+1])
                          NULL,                                       // veloloctab (no vertex weights)
                          NULL,                                       // vlblloctab (no vertex labels)
                          static_cast<SCOTCH_Num>(s),                // edgelocnbr
                          static_cast<SCOTCH_Num>(s),                // edgelocsiz
                          &scotch_edgeloctab[0],                     // edgeloctab
                          NULL,                                       // edgegsttab (let Scotch compute)
                          NULL);                                      // edloloctab (no edge weights)

        // Check graph consistency
        SCOTCH_dgraphCheck(&grafdat);

        // Create partitioning strategy (use default for now)
        SCOTCH_Strat stradat;
        SCOTCH_stratInit(&stradat);

        // Perform partitioning
        SCOTCH_dgraphPart(&grafdat, static_cast<SCOTCH_Num>(mpiSize),
                         &stradat, &scotch_parttab[0]);

        // Copy partition results back to index_t array
        for (index_t i = 0; i < myNumVertices; ++i)
            partition[i] = static_cast<index_t>(scotch_parttab[i]);

        // Clean up Scotch structures
        SCOTCH_stratExit(&stradat);
        SCOTCH_dgraphExit(&grafdat);
#endif
        delete[] xyz;
        delete[] index;
        delete[] ptr;
    } else {
        // Graph partitioning not available or not applicable
#pragma omp parallel for
        for (index_t i = 0; i < myNumVertices; ++i)
            partition[i] = myRank; // CPU myRank owns all
    }
#else
    // No graph partitioning library available
#pragma omp parallel for
    for (index_t i = 0; i < myNumVertices; ++i)
        partition[i] = myRank; // CPU myRank owns all
#endif // ESYS_HAVE_PARMETIS || ESYS_HAVE_SCOTCH

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
                  mpiSize, MPI_DIM_T, m_mpiInfo->comm);
#else
    for (int i = 0; i < mpiSize; ++i)
        recvbuf[i] = new_distribution[i];
#endif
    new_distribution[0] = 0;
    std::vector<index_t> newGlobalDOFID(len);
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
    int dest = m_mpiInfo->mod_rank(myRank + 1);
    int source = m_mpiInfo->mod_rank(myRank - 1);
#endif
    int current_rank = myRank;
    std::vector<short> setNewDOFId(numNodes, 1);

    for (int p = 0; p < mpiSize; ++p) {
        const index_t firstVertex = distribution[current_rank];
        const index_t lastVertex = distribution[current_rank + 1];
#pragma omp parallel for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = m_nodes->globalDegreesOfFreedom[i];
            if (setNewDOFId[i] && firstVertex <= k && k < lastVertex) {
                m_nodes->globalDegreesOfFreedom[i] = newGlobalDOFID[k - firstVertex];
                setNewDOFId[i] = 0;
            }
        }

        if (p < mpiSize - 1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_DIM_T,
                                 dest, m_mpiInfo->counter(),
                                 source, m_mpiInfo->counter(),
                                 m_mpiInfo->comm, &status);
            m_mpiInfo->incCounter();
#endif
            current_rank = m_mpiInfo->mod_rank(current_rank - 1);
        }
    }
    for (int i = 0; i < mpiSize + 1; ++i)
        distribution[i] = new_distribution[i];

    delete[] partition;
}

} // namespace finley

