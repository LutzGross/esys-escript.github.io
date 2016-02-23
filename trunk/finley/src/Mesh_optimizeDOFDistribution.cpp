
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


/****************************************************************************

  Finley: Mesh: optimizes the distribution of DOFs across processors
  using ParMETIS. On return a new distribution is given and the globalDOF
  are relabeled accordingly but the mesh is not redistributed yet.

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "Mesh.h"
#include "IndexList.h"
#ifdef USE_PARMETIS
#include <parmetis.h>
#ifndef REALTYPEWIDTH
typedef float real_t;
#endif
#endif
#include <boost/scoped_array.hpp>

namespace finley {

#ifdef USE_PARMETIS
// Checks whether there is any rank which has no vertex. In case 
// such a rank exists, we don't use parmetis since parmetis requires
// that every rank has at least 1 vertex (at line 129 of file
// "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
// any rank has no vertex).
bool allRanksHaveNodes(esysUtils::JMPI& mpiInfo, const std::vector<index_t>& distribution)
{
    int ret = 1;

    if (mpiInfo->rank == 0) {
        for (int i=0; i<mpiInfo->size; i++) {
            if (distribution[i+1] == distribution[i]) {
                ret = 0;
                break;
            }
        }
        if (ret == 0) {
            std::cout << "Mesh::optimizeDOFDistribution: "
                << "Parmetis is not used since at least one rank has no vertex!"
                << std::endl;
        }
    }
    MPI_Bcast(&ret, 1, MPI_INTEGER, 0, mpiInfo->comm);
    return (ret==1);
}
#endif


/****************************************************************************/

void Mesh::optimizeDOFDistribution(std::vector<index_t>& distribution)
{
    // these two are not const because of parmetis call
    int mpiSize=MPIInfo->size;
    const int myRank=MPIInfo->rank;
    const index_t myFirstVertex=distribution[myRank];
    const index_t myLastVertex=distribution[myRank+1];
    const dim_t myNumVertices=myLastVertex-myFirstVertex;

    // first step is to distribute the elements according to a global X of DOF
    // len is used for the sending around of partition later on
    index_t len=0;
    for (int p=0; p<mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);
    std::vector<index_t> partition(len);

#ifdef USE_PARMETIS
    if (mpiSize>1 && allRanksHaveNodes(MPIInfo, distribution)) {
        boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
        int dim=Nodes->numDim;
        // create the adjacency structure xadj and adjncy
#pragma omp parallel
        {
            // insert contributions from element matrices into columns index
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, Elements,
                    Nodes->globalDegreesOfFreedom, Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, FaceElements,
                    Nodes->globalDegreesOfFreedom, Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, ContactElements,
                    Nodes->globalDegreesOfFreedom, Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                    myFirstVertex, myLastVertex, Points,
                    Nodes->globalDegreesOfFreedom, Nodes->globalDegreesOfFreedom);
        }
       
        // create the local matrix pattern
        const dim_t globalNumVertices=distribution[mpiSize];
        paso::Pattern_ptr pattern(paso::Pattern::fromIndexListArray(0,
                myNumVertices, index_list.get(), 0, globalNumVertices, 0));
        // set the coordinates
        std::vector<real_t> xyz(myNumVertices*dim);
#pragma omp parallel for
        for (index_t i=0; i<Nodes->numNodes; ++i) {
            const index_t k=Nodes->globalDegreesOfFreedom[i]-myFirstVertex;
            if (k>=0 && k<myNumVertices) {
                for (int j=0; j<dim; ++j)
                    xyz[k*dim+j]=static_cast<real_t>(Nodes->Coordinates[INDEX2(j,i,dim)]); 
            }
        }

        index_t wgtflag = 0;
        index_t numflag = 0;
        index_t ncon = 1;
        index_t edgecut;
        index_t impiSize = mpiSize;
        index_t idim = dim;
        // options[0]=1 -> non-default values, evaluate rest of options
        // options[1]=15 -> DBG_TIME | DBG_INFO | DBG_PROGRESS | DBG_REFINEINFO
        // options[2] -> random seed
        index_t options[3] = { 1, 15, 0 };
        std::vector<real_t> tpwgts(ncon*mpiSize, 1.f/mpiSize);
        std::vector<real_t> ubvec(ncon, 1.05f);
        ParMETIS_V3_PartGeomKway(&distribution[0], pattern->ptr, pattern->index,
                              NULL, NULL, &wgtflag, &numflag, &idim, &xyz[0],
                              &ncon, &impiSize, &tpwgts[0], &ubvec[0], options,
                              &edgecut, &partition[0], &MPIInfo->comm);
    } else {
        for (index_t i=0; i<myNumVertices; ++i)
            partition[i]=0; // CPU 0 owns all
    }
#else
    for (index_t i=0; i<myNumVertices; ++i)
        partition[i]=myRank; // CPU myRank owns all
#endif

    // create a new distribution and labeling of the DOF
    std::vector<index_t> new_distribution(mpiSize+1, 0);
#pragma omp parallel
    {
        std::vector<int> loc_partition_count(mpiSize, 0);
#pragma omp for
        for (index_t i=0; i<myNumVertices; ++i)
            loc_partition_count[partition[i]]++;
#pragma omp critical
        {
            for (int i=0; i<mpiSize; ++i)
                new_distribution[i]+=loc_partition_count[i];
        }
    }
    index_t *recvbuf=new index_t[mpiSize*mpiSize];
#ifdef ESYS_MPI
    // recvbuf will be the concatenation of each CPU's contribution to
    // new_distribution
    MPI_Allgather(&new_distribution[0], mpiSize, MPI_DIM_T, recvbuf, mpiSize,
                  MPI_INT, MPIInfo->comm);
#else
    for (int i=0; i<mpiSize; ++i)
        recvbuf[i]=new_distribution[i];
#endif
    new_distribution[0]=0;
    std::vector<index_t> newGlobalDOFID(len);
    for (int rank=0; rank<mpiSize; rank++) {
        index_t c=0;
        for (int i=0; i<myRank; ++i)
            c+=recvbuf[rank+mpiSize*i];
        for (index_t i=0; i<myNumVertices; ++i) {
            if (rank==partition[i]) {
                newGlobalDOFID[i]=new_distribution[rank]+c;
                c++;
            }
        }
        for (int i=myRank+1; i<mpiSize; ++i)
            c+=recvbuf[rank+mpiSize*i];
        new_distribution[rank+1]=new_distribution[rank]+c;
    }
    delete[] recvbuf;

    // now the overlap needs to be created by sending the partition around
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(myRank + 1);
    int source = MPIInfo->mod_rank(myRank - 1);
#endif
    int current_rank=myRank;
    std::vector<short> setNewDOFId(Nodes->numNodes, 1);

    for (int p=0; p<mpiSize; ++p) {
        const index_t firstVertex=distribution[current_rank];
        const index_t lastVertex=distribution[current_rank+1];
#pragma omp parallel for
        for (index_t i=0; i<Nodes->numNodes; ++i) {
            const index_t k=Nodes->globalDegreesOfFreedom[i];
            if (setNewDOFId[i] && (firstVertex<=k) && (k<lastVertex)) {
                Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                setNewDOFId[i]=0;
            }
        }

        if (p<mpiSize-1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_DIM_T,
                               dest, MPIInfo->counter(),
                               source, MPIInfo->counter(),
                               MPIInfo->comm, &status);
            MPIInfo->incCounter();
#endif
            current_rank=MPIInfo->mod_rank(current_rank-1);
        }
    }
    for (int i=0; i<mpiSize+1; ++i)
        distribution[i]=new_distribution[i];
}

} // namespace finley

