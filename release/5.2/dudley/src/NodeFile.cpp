
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

namespace dudley {

static std::pair<index_t,index_t> getGlobalRange(dim_t n, const index_t* id,
                                                 escript::JMPI mpiInfo)
{
    std::pair<index_t,index_t> result(util::getMinMaxInt(1, n, id));

#ifdef ESYS_MPI
    index_t global_id_range[2];
    index_t id_range[2] = { -result.first, result.second };
    MPI_Allreduce(id_range, global_id_range, 2, MPI_DIM_T, MPI_MAX,
                  mpiInfo->comm);
    result.first = -global_id_range[0];
    result.second = global_id_range[1];
#endif
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

NodeFile::NodeFile(int nDim, escript::JMPI mpiInfo) :
    numNodes(0),
    MPIInfo(mpiInfo),
    numDim(nDim),
    Id(NULL),
    Tag(NULL),
    globalDegreesOfFreedom(NULL),
    Coordinates(NULL),
    globalNodesIndex(NULL),
    degreesOfFreedomId(NULL),
    status(DUDLEY_INITIAL_STATUS)
{
}

NodeFile::~NodeFile()
{
    freeTable();
}

void NodeFile::allocTable(dim_t NN)
{
    if (numNodes > 0)
        freeTable();

    Id = new index_t[NN];
    Coordinates = new escript::DataTypes::real_t[NN*numDim];
    Tag = new int[NN];
    globalDegreesOfFreedom = new index_t[NN];
    globalNodesIndex = new index_t[NN];
    degreesOfFreedomId = new index_t[NN];
    numNodes = NN;

    // this initialization makes sure that data are located on the right
    // processor
#pragma omp parallel for
    for (index_t n=0; n<numNodes; n++) {
        Id[n] = -1;
        for (int i=0; i<numDim; i++)
            Coordinates[INDEX2(i,n,numDim)] = 0.;
        Tag[n] = -1;
        globalDegreesOfFreedom[n] = -1;
        globalNodesIndex[n] = -1;
        degreesOfFreedomId[n] = -1;
    }
}

void NodeFile::freeTable()
{
    delete[] Id;
    delete[] Coordinates;
    delete[] globalDegreesOfFreedom;
    delete[] globalNodesIndex;
    delete[] Tag;
    delete[] degreesOfFreedomId;
    nodesMapping.clear();
    degreesOfFreedomMapping.clear();
    nodesDistribution.reset();
    dofDistribution.reset();
#ifdef ESYS_HAVE_PASO
    degreesOfFreedomConnector.reset();
#endif
    numNodes = 0;
}

void NodeFile::print() const
{
    std::cout << "=== " << numDim << "D-Nodes:\nnumber of nodes=" << numNodes
        << std::endl;
    std::cout << "Id,Tag,globalDegreesOfFreedom,degreesOfFreedom,node,Coordinates" << std::endl;
    for (index_t i=0; i<numNodes; i++) {
        std::cout << Id[i] << "," << Tag[i] << "," << globalDegreesOfFreedom[i]
            << "," << degreesOfFreedomMapping.target[i]
            << "," << nodesMapping.target[i] << " ";
        std::cout.precision(15);
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        for (int j=0; j<numDim; j++)
            std:: cout << Coordinates[INDEX2(j,i,numDim)];
        std::cout << std::endl;
    }
}

void NodeFile::copyTable(index_t offset, index_t idOffset, index_t dofOffset,
                         const NodeFile* in)
{
    // check number of dimensions and table size
    if (numDim != in->numDim)
        throw escript::ValueError("NodeFile::copyTable: dimensions of node files don't match");

    if (numNodes < in->numNodes + offset)
        throw escript::ValueError("NodeFile::copyTable: node table is too small.");

#pragma omp parallel for
    for (index_t n = 0; n < in->numNodes; n++) {
        Id[offset + n] = in->Id[n] + idOffset;
        Tag[offset + n] = in->Tag[n];
        globalDegreesOfFreedom[offset + n] = in->globalDegreesOfFreedom[n] + dofOffset;
        for (int i = 0; i < numDim; i++)
            Coordinates[INDEX2(i, offset + n, numDim)] =
                                    in->Coordinates[INDEX2(i, n, in->numDim)];
    }
}

std::pair<index_t,index_t> NodeFile::getDOFRange() const
{
    std::pair<index_t,index_t> result(util::getMinMaxInt(
                                        1, numNodes, globalDegreesOfFreedom));
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

std::pair<index_t,index_t> NodeFile::getGlobalIdRange() const
{
    return getGlobalRange(numNodes, Id, MPIInfo);
}

std::pair<index_t,index_t> NodeFile::getGlobalDOFRange() const
{
    return getGlobalRange(numNodes, globalDegreesOfFreedom, MPIInfo);
}

std::pair<index_t,index_t> NodeFile::getGlobalNodeIDIndexRange() const
{
    return getGlobalRange(numNodes, globalNodesIndex, MPIInfo);
}

void NodeFile::setCoordinates(const escript::Data& newX)
{
    if (newX.isComplex())
    {
        throw escript::ValueError("NodeFile::setCoordinates: argument can not be complex.");
    }
    if (newX.getDataPointSize() != numDim) {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of dimensions of new "
            "coordinates has to be " << numDim;
        throw escript::ValueError(ss.str());
    } else if (newX.getNumDataPointsPerSample() != 1 ||
            newX.getNumSamples() != numNodes) {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of given nodes must be "
            << numNodes;
        throw escript::ValueError(ss.str());
    } else {
        const size_t numDim_size = numDim * sizeof(double);
        ++status;
	escript::DataTypes::real_t wantreal=0;
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            memcpy(&Coordinates[INDEX2(0, n, numDim)],
                    newX.getSampleDataRO(n, wantreal), numDim_size);
        }
    }
}

void NodeFile::setTags(int newTag, const escript::Data& mask)
{
  
    if (1 != mask.getDataPointSize()) {
        throw escript::ValueError("NodeFile::setTags: number of components of mask must be 1.");
    } else if (mask.getNumDataPointsPerSample() != 1 ||
            mask.getNumSamples() != numNodes) {
        throw escript::ValueError("NodeFile::setTags: illegal number of samples of mask Data object");
    }
    
    escript::DataTypes::real_t wantreal=0;
#pragma omp parallel for
    for (index_t n = 0; n < numNodes; n++) {
        if (mask.getSampleDataRO(n, wantreal)[0] > 0)
            Tag[n] = newTag;
    }
    updateTagList();
}

void NodeFile::assignMPIRankToDOFs(int* mpiRankOfDOF,
                                   const IndexVector& distribution)
{
    int p_min = MPIInfo->size, p_max = -1;
    // first we calculate the min and max DOF on this processor to reduce
    // costs for searching
    const std::pair<index_t,index_t> dofRange(getDOFRange());

    for (int p = 0; p < MPIInfo->size; ++p) {
        if (distribution[p] <= dofRange.first)
            p_min = p;
        if (distribution[p] <= dofRange.second)
            p_max = p;
    }
#pragma omp parallel for
    for (index_t n = 0; n < numNodes; ++n) {
        const index_t k = globalDegreesOfFreedom[n];
        for (int p = p_min; p <= p_max; ++p) {
            if (k < distribution[p + 1]) {
                mpiRankOfDOF[n] = p;
                break;
            }
        }
    }
}

} // namespace dudley

