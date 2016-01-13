
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


/****************************************************************************

  Finley: Mesh: optimizes the distribution of DOFs across processors
  using ParMETIS. On return a new distribution is given and the globalDOF
  are relabeled accordingly but the mesh is not redistributed yet.

*****************************************************************************/

#include "Mesh.h"
#include "IndexList.h"
#ifdef USE_PARMETIS
#include <parmetis.h>
#endif
#include <boost/scoped_array.hpp>

namespace finley {

#ifdef USE_PARMETIS
// Checks whether there is any rank which has no vertex. In case 
// such a rank exists, we don't use parmetis since parmetis requires
// that every rank has at least 1 vertex (at line 129 of file
// "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
// any rank has no vertex).
bool allRanksHaveNodes(esysUtils::JMPI& mpiInfo, const std::vector<int>& distribution)
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

void Mesh::optimizeDOFDistribution(std::vector<int>& distribution)
{
    // these two are not const because of parmetis call
    int mpiSize=MPIInfo->size;
    const int myRank=MPIInfo->rank;
    const int myFirstVertex=distribution[myRank];
    const int myLastVertex=distribution[myRank+1];
    const int myNumVertices=myLastVertex-myFirstVertex;

    // first step is to distribute the elements according to a global X of DOF
    // len is used for the sending around of partition later on
    int len=0;
    for (int p=0; p<mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);
    std::vector<int> partition(len);

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
        const int globalNumVertices=distribution[mpiSize];
        paso::Pattern_ptr pattern(paso::Pattern::fromIndexListArray(0,
                myNumVertices, index_list.get(), 0, globalNumVertices, 0));
        // set the coordinates
        std::vector<float> xyz(myNumVertices*dim);
#pragma omp parallel for
        for (int i=0; i<Nodes->numNodes; ++i) {
            const int k=Nodes->globalDegreesOfFreedom[i]-myFirstVertex;
            if (k>=0 && k<myNumVertices) {
                for (int j=0; j<dim; ++j)
                    xyz[k*dim+j]=static_cast<float>(Nodes->Coordinates[INDEX2(j,i,dim)]); 
            }
        }

        int wgtflag = 0;
        int numflag = 0;
        int ncon = 1;
        int edgecut;
        int options[2] = { 3, 15 };
        std::vector<float> tpwgts(ncon*mpiSize, 1.f/mpiSize);
        std::vector<float> ubvec(ncon, 1.05f);
        ParMETIS_V3_PartGeomKway(&distribution[0], pattern->ptr, pattern->index,
                              NULL, NULL, &wgtflag, &numflag, &dim, &xyz[0],
                              &ncon, &mpiSize, &tpwgts[0], &ubvec[0], options,
                              &edgecut, &partition[0], &MPIInfo->comm);
    } else {
        for (int i=0; i<myNumVertices; ++i)
            partition[i]=0; // CPU 0 owns all
    }
#else
    for (int i=0; i<myNumVertices; ++i)
        partition[i]=myRank; // CPU myRank owns all
#endif

    // create a new distribution and labeling of the DOF
    std::vector<int> new_distribution(mpiSize+1, 0);
#pragma omp parallel
    {
        std::vector<int> loc_partition_count(mpiSize, 0);
#pragma omp for
        for (int i=0; i<myNumVertices; ++i)
            loc_partition_count[partition[i]]++;
#pragma omp critical
        {
            for (int i=0; i<mpiSize; ++i)
                new_distribution[i]+=loc_partition_count[i];
        }
    }
    int *recvbuf=new int[mpiSize*mpiSize];
#ifdef ESYS_MPI
    // recvbuf will be the concatenation of each CPU's contribution to
    // new_distribution
    MPI_Allgather(&new_distribution[0], mpiSize, MPI_INT, recvbuf, mpiSize,
                  MPI_INT, MPIInfo->comm);
#else
    for (int i=0; i<mpiSize; ++i)
        recvbuf[i]=new_distribution[i];
#endif
    new_distribution[0]=0;
    std::vector<int> newGlobalDOFID(len);
    for (int rank=0; rank<mpiSize; rank++) {
        int c=0;
        for (int i=0; i<myRank; ++i)
            c+=recvbuf[rank+mpiSize*i];
        for (int i=0; i<myNumVertices; ++i) {
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
    int dest=esysUtils::mod_rank(mpiSize, myRank + 1);
    int source=esysUtils::mod_rank(mpiSize, myRank - 1);
#endif
    int current_rank=myRank;
    std::vector<short> setNewDOFId(Nodes->numNodes, 1);

    for (int p=0; p<mpiSize; ++p) {
        const int firstVertex=distribution[current_rank];
        const int lastVertex=distribution[current_rank+1];
#pragma omp parallel for
        for (int i=0; i<Nodes->numNodes; ++i) {
            const int k=Nodes->globalDegreesOfFreedom[i];
            if (setNewDOFId[i] && (firstVertex<=k) && (k<lastVertex)) {
                Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                setNewDOFId[i]=0;
            }
        }

        if (p<mpiSize-1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_INT,
                               dest, MPIInfo->msg_tag_counter,
                               source, MPIInfo->msg_tag_counter,
                               MPIInfo->comm, &status);
#endif
            MPIInfo->msg_tag_counter++;
            current_rank=esysUtils::mod_rank(mpiSize, current_rank-1);
        }
    }
    for (int i=0; i<mpiSize+1; ++i)
        distribution[i]=new_distribution[i];
}

} // namespace finley

