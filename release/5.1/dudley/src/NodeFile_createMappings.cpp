
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

/*   Dudley: NodeFile : creates the mappings using the indexReducedNodes */
/*                 no distribution is happening                          */

/****************************************************************************/

#include "NodeFile.h"

namespace dudley {

void NodeFile::createDOFMappingAndCoupling()
{
    const index_t myFirstDOF = dofDistribution->getFirstComponent();
    const index_t myLastDOF = dofDistribution->getLastComponent();
    const int mpiSize = MPIInfo->size;
    const int myRank = MPIInfo->rank;

    index_t min_DOF, max_DOF;
    std::pair<index_t,index_t> DOF_range(util::getFlaggedMinMaxInt(
                                            numNodes, globalDegreesOfFreedom, -1));

    if (DOF_range.second < DOF_range.first) {
        min_DOF = myFirstDOF;
        max_DOF = myLastDOF - 1;
    } else {
        min_DOF = DOF_range.first;
        max_DOF = DOF_range.second;
    }

    int p_min = mpiSize;
    int p_max = -1;
    if (max_DOF >= min_DOF) {
        for (int p = 0; p < mpiSize; ++p) {
            if (dofDistribution->first_component[p] <= min_DOF)
                p_min = p;
            if (dofDistribution->first_component[p] <= max_DOF)
                p_max = p;
        }
    }

    std::stringstream ss;
    if (myFirstDOF<myLastDOF && !(min_DOF <= myFirstDOF && myLastDOF-1 <= max_DOF)) {
        ss << "createDOFMappingAndCoupling: Local elements do not span local "
              "degrees of freedom. min_DOF=" << min_DOF << ", myFirstDOF="
           << myFirstDOF << ", myLastDOF-1=" << myLastDOF-1
           << ", max_DOF=" << max_DOF << " on rank=" << MPIInfo->rank;
    }
    const std::string msg(ss.str());
    int error = msg.length();
    int gerror = error;
    escript::checkResult(error, gerror, MPIInfo);
    if (gerror > 0) {
        char* gmsg;
        escript::shipString(msg.c_str(), &gmsg, MPIInfo->comm);
        throw DudleyException(gmsg);
    }

    const index_t UNUSED = -1;
    const dim_t len_loc_dof = max_DOF - min_DOF + 1;
    index_t* shared = new index_t[numNodes * (p_max - p_min + 1)];
    index_t* locDOFMask = new index_t[len_loc_dof];
    index_t* nodeMask = new index_t[numNodes];
#ifdef BOUNDS_CHECK
    ESYS_ASSERT(myLastDOF-min_DOF <= len_loc_dof, "BOUNDS_CHECK");
#endif

#pragma omp parallel
    {
#pragma omp for
        for (index_t i = 0; i < len_loc_dof; ++i)
            locDOFMask[i] = UNUSED;
#pragma omp for
        for (index_t i = 0; i < numNodes; ++i)
            nodeMask[i] = UNUSED;
#pragma omp for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = globalDegreesOfFreedom[i];
            if (k > -1) {
#ifdef BOUNDS_CHECK
                ESYS_ASSERT(k-min_DOF < len_loc_dof, "BOUNDS_CHECK");
#endif
                locDOFMask[k - min_DOF] = UNUSED - 1;
            }
        }
#pragma omp for
        for (index_t i = myFirstDOF - min_DOF; i < myLastDOF - min_DOF; ++i) {
            locDOFMask[i] = i - myFirstDOF + min_DOF;
        }
    }

    index_t* wanted_DOFs = new index_t[numNodes];
    std::vector<index_t> rcv_len(mpiSize);
    std::vector<index_t> snd_len(mpiSize);
    std::vector<int> neighbour;
    std::vector<index_t> offsetInShared;
    dim_t n = 0;
    dim_t lastn = n;

    for (int p = p_min; p <= p_max; ++p) {
        if (p != myRank) {
            const index_t firstDOF = std::max(min_DOF, dofDistribution->first_component[p]);
            const index_t lastDOF = std::min(max_DOF + 1, dofDistribution->first_component[p + 1]);
#ifdef BOUNDS_CHECK
            ESYS_ASSERT(lastDOF-min_DOF <= len_loc_dof, "BOUNDS_CHECK");
#endif
            for (index_t i = firstDOF - min_DOF; i < lastDOF - min_DOF; ++i) {
                if (locDOFMask[i] == UNUSED - 1) {
                    locDOFMask[i] = myLastDOF - myFirstDOF + n;
                    wanted_DOFs[n] = i + min_DOF;
                    ++n;
                }
            }
            if (n > lastn) {
                rcv_len[p] = n - lastn;
                neighbour.push_back(p);
                offsetInShared.push_back(lastn);
                lastn = n;
            }
        } // if p!=myRank
    } // for p

    offsetInShared.push_back(lastn);

    // assign new DOF labels to nodes
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i) {
        const index_t k = globalDegreesOfFreedom[i];
        if (k > -1)
            nodeMask[i] = locDOFMask[k - min_DOF];
    }

    degreesOfFreedomMapping.assign(nodeMask, numNodes, UNUSED);

    // define how to get DOF values for controlled but other processors
#ifdef BOUNDS_CHECK
    ESYS_ASSERT(numNodes == 0 || offsetInShared.back() < numNodes * (p_max - p_min + 1), "BOUNDS_CHECK");
#endif
#pragma omp parallel for
    for (index_t i = 0; i < lastn; ++i)
        shared[i] = myLastDOF - myFirstDOF + i;

#ifdef ESYS_HAVE_PASO
    paso::SharedComponents_ptr rcv_shcomp(new paso::SharedComponents(
                                    myLastDOF - myFirstDOF, neighbour, shared,
                                    offsetInShared));
#endif

    /////////////////////////////////
    //   now we build the sender   //
    /////////////////////////////////
#ifdef ESYS_MPI
    std::vector<MPI_Request> mpi_requests(mpiSize * 2);
    std::vector<MPI_Status> mpi_stati(mpiSize * 2);
    MPI_Alltoall(&rcv_len[0], 1, MPI_DIM_T, &snd_len[0], 1, MPI_DIM_T,
                 MPIInfo->comm);
    int count = 0;
    for (int p = 0; p < neighbour.size(); p++) {
        MPI_Isend(&wanted_DOFs[offsetInShared[p]],
                offsetInShared[p+1] - offsetInShared[p],
                MPI_DIM_T, neighbour[p], MPIInfo->counter() + myRank,
                MPIInfo->comm, &mpi_requests[count]);
        count++;
    }
#else
    snd_len[0] = rcv_len[0];
#endif
    n = 0;
    neighbour.clear();
    offsetInShared.clear();
#ifdef ESYS_MPI
    for (int p = 0; p < mpiSize; p++) {
        if (snd_len[p] > 0) {
            MPI_Irecv(&shared[n], snd_len[p], MPI_DIM_T, p,
                      MPIInfo->counter() + p, MPIInfo->comm,
                      &mpi_requests[count]);
            count++;
            neighbour.push_back(p);
            offsetInShared.push_back(n);
            n += snd_len[p];
        }
    }
    MPIInfo->incCounter(MPIInfo->size);
    MPI_Waitall(count, &mpi_requests[0], &mpi_stati[0]);
#endif
    offsetInShared.push_back(n);

    // map global IDs to local IDs
#pragma omp parallel for
    for (index_t i = 0; i < n; ++i) {
        shared[i] = locDOFMask[shared[i] - min_DOF];
    }

#ifdef ESYS_HAVE_PASO
    paso::SharedComponents_ptr snd_shcomp(new paso::SharedComponents(
                                    myLastDOF - myFirstDOF, neighbour, shared,
                                    offsetInShared));
    degreesOfFreedomConnector.reset(new paso::Connector(snd_shcomp, rcv_shcomp));
#endif

    delete[] wanted_DOFs;
    delete[] nodeMask;
    delete[] shared;
    delete[] locDOFMask;
}

void NodeFile::createNodeMappings(const IndexVector& dofDist,
                                  const IndexVector& nodeDist)
{
    // ==== distribution of Nodes ====
    nodesDistribution.reset(new escript::Distribution(MPIInfo, nodeDist));

    // ==== distribution of DOFs ====
    dofDistribution.reset(new escript::Distribution(MPIInfo, dofDist));

    index_t* nodeMask = new index_t[numNodes];
    const index_t UNUSED = -1;

    // ==== nodes mapping (dummy) ====
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i)
        nodeMask[i] = i;
    nodesMapping.assign(nodeMask, numNodes, UNUSED);

    // ==== mapping between nodes and DOFs + DOF connector ====
    createDOFMappingAndCoupling();

    // get the IDs for DOFs
#pragma omp parallel for
    for (index_t i = 0; i < degreesOfFreedomMapping.numTargets; ++i)
        degreesOfFreedomId[i] = Id[degreesOfFreedomMapping.map[i]];

    delete[] nodeMask;
}

} // namespace dudley

