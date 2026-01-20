
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

#include "NodeFile.h"

#include <escript/Data.h>
#include <escript/index.h>

#include <limits>
#include <sstream>
#include <iostream>

namespace finley {

// helper function
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

// helper function
static void scatterEntries(dim_t n, const index_t* index, index_t min_index,
                           index_t max_index, index_t* Id_out,
                           const index_t* Id_in,
                           int* Tag_out, const int* Tag_in,
                           index_t* globalDegreesOfFreedom_out,
                           const index_t* globalDegreesOfFreedom_in,
                           int numDim, double* Coordinates_out,
                           const double* Coordinates_in)
{
    const index_t range = max_index-min_index;
    const size_t numDim_size = numDim*sizeof(double);

#pragma omp parallel for
    for (index_t i=0; i<n; i++) {
        const index_t k = index[i]-min_index;
        if (k>=0 && k<range) {
            Id_out[k] = Id_in[i];
            Tag_out[k] = Tag_in[i];
            globalDegreesOfFreedom_out[k] = globalDegreesOfFreedom_in[i];
            memcpy(&Coordinates_out[INDEX2(0,k,numDim)],
                    &Coordinates_in[INDEX2(0,i,numDim)], numDim_size);
        }
    }
}

// helper function
static void gatherEntries(dim_t n, const index_t* index,
                          index_t min_index, index_t max_index,
                          index_t* Id_out, const index_t* Id_in,
                          int* Tag_out, const int* Tag_in,
                          index_t* globalDegreesOfFreedom_out,
                          const index_t* globalDegreesOfFreedom_in,
                          int numDim, double* Coordinates_out,
                          const double* Coordinates_in)
{
    const index_t range = max_index-min_index;
    const size_t numDim_size = numDim*sizeof(double);

#pragma omp parallel for
    for (index_t i=0; i<n; i++) {
        const index_t k = index[i]-min_index;
        if (k>=0 && k<range) {
            Id_out[i] = Id_in[k];
            Tag_out[i] = Tag_in[k];
            globalDegreesOfFreedom_out[i] = globalDegreesOfFreedom_in[k];
            memcpy(&Coordinates_out[INDEX2(0,i,numDim)],
                    &Coordinates_in[INDEX2(0,k,numDim)], numDim_size);
        }
    }
}

/// constructor
/// use NodeFile::allocTable to allocate the node table (Id,Coordinates)
NodeFile::NodeFile(int nDim, escript::JMPI mpiInfo) :
    numNodes(0),
    MPIInfo(mpiInfo),
    numDim(nDim),
    Id(NULL),
    Tag(NULL),
    globalDegreesOfFreedom(NULL),
    Coordinates(NULL),
    globalReducedDOFIndex(NULL),
    globalReducedNodesIndex(NULL),
    globalNodesIndex(NULL),
    reducedNodesId(NULL),
    degreesOfFreedomId(NULL),
    reducedDegreesOfFreedomId(NULL),
    status(FINLEY_INITIAL_STATUS)
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
    globalReducedDOFIndex = new index_t[NN];
    globalReducedNodesIndex = new index_t[NN];
    globalNodesIndex = new index_t[NN];
    reducedNodesId = new index_t[NN];
    degreesOfFreedomId = new index_t[NN];
    reducedDegreesOfFreedomId = new index_t[NN];
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
        globalReducedDOFIndex[n] = -1;
        globalReducedNodesIndex[n] = -1;
        globalNodesIndex[n] = -1;
        reducedNodesId[n] = -1;
        degreesOfFreedomId[n] = -1;
        reducedDegreesOfFreedomId[n] = -1;
    }
}

void NodeFile::freeTable()
{
    delete[] Id;
    delete[] Coordinates;
    delete[] globalDegreesOfFreedom;
    delete[] globalReducedDOFIndex;
    delete[] globalReducedNodesIndex;
    delete[] globalNodesIndex;
    delete[] Tag;
    delete[] reducedNodesId;
    delete[] degreesOfFreedomId;
    delete[] reducedDegreesOfFreedomId;
    tagsInUse.clear();
    nodesMapping.clear();
    reducedNodesMapping.clear();
    degreesOfFreedomMapping.clear();
    reducedDegreesOfFreedomMapping.clear();
    nodesDistribution.reset();
    reducedNodesDistribution.reset();
    degreesOfFreedomDistribution.reset();
    reducedDegreesOfFreedomDistribution.reset();
#ifdef ESYS_HAVE_PASO
    degreesOfFreedomConnector.reset();
    reducedDegreesOfFreedomConnector.reset();
#endif
#ifdef ESYS_HAVE_TRILINOS
    trilinosRowMap.reset();
    trilinosReducedRowMap.reset();
    trilinosColMap.reset();
    trilinosReducedColMap.reset();
#endif
    numNodes = 0;
}

void NodeFile::print() const
{
    std::cout << "=== " << numDim << "D-Nodes:\nnumber of nodes=" << numNodes
        << std::endl;
    std::cout << "Id,Tag,globalDegreesOfFreedom,degreesOfFreedom,reducedDegreesOfFeedom,node,reducedNode,Coordinates" << std::endl;
    for (index_t i=0; i<numNodes; i++) {
        std::cout << Id[i] << "," << Tag[i] << "," << globalDegreesOfFreedom[i]
            << "," << degreesOfFreedomMapping.target[i]
            << "," << reducedDegreesOfFreedomMapping.target[i]
            << "," << nodesMapping.target[i] << reducedNodesMapping.target[i]
            << " ";
        std::cout.precision(15);
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        for (int j=0; j<numDim; j++)
            std:: cout << Coordinates[INDEX2(j,i,numDim)];
        std::cout << std::endl;
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

/// copies the array newX into this->coordinates
void NodeFile::setCoordinates(const escript::Data& newX)
{
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
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            memcpy(&Coordinates[INDEX2(0, n, numDim)],
                    newX.getSampleDataRO(n), numDim_size);
        }
    }
}

/// sets tags to newTag where mask>0
void NodeFile::setTags(int newTag, const escript::Data& mask)
{
    if (1 != mask.getDataPointSize()) {
        throw escript::ValueError("NodeFile::setTags: number of components of mask must be 1.");
    } else if (mask.getNumDataPointsPerSample() != 1 ||
            mask.getNumSamples() != numNodes) {
        throw escript::ValueError("NodeFile::setTags: illegal number of samples of mask Data object");
    }

#pragma omp parallel for
    for (index_t n = 0; n < numNodes; n++) {
        if (mask.getSampleDataRO(n)[0] > 0)
            Tag[n] = newTag;
    }
    updateTagList();
}


void NodeFile::copyTable(index_t offset, index_t idOffset, index_t dofOffset,
                         const NodeFile* in)
{
    // check number of dimensions and table size
    if (numDim != in->numDim) {
        throw escript::ValueError("NodeFile::copyTable: dimensions of node files don't match");
    }
    if (numNodes < in->numNodes+offset) {
        throw escript::ValueError("NodeFile::copyTable: node table is too small.");
    }

#pragma omp parallel for
    for (index_t n=0; n<in->numNodes; n++) {
        Id[offset+n]=in->Id[n]+idOffset;
        Tag[offset+n]=in->Tag[n];
        globalDegreesOfFreedom[offset+n]=in->globalDegreesOfFreedom[n]+dofOffset;
        for(int i=0; i<numDim; i++)
            Coordinates[INDEX2(i, offset+n, numDim)] =
                                    in->Coordinates[INDEX2(i, n, in->numDim)];
    }
}

/// scatters the NodeFile in into this NodeFile using index[0:in->numNodes-1].
/// index has to be between 0 and numNodes-1.
/// colouring is chosen for the worst case
void NodeFile::scatter(const index_t* index, const NodeFile* in)
{
    scatterEntries(numNodes, index, 0, in->numNodes, Id, in->Id, Tag, in->Tag,
                   globalDegreesOfFreedom, in->globalDegreesOfFreedom,
                   numDim, Coordinates, in->Coordinates);
}

/// gathers this NodeFile from the NodeFile 'in' using the entries in
/// index[0:out->numNodes-1] which are between 0 (and in->numNodes)
/// (exclusive)
// WARNING: This does not wotj for MPI!!!
void NodeFile::gather(const index_t* index, const NodeFile* in)
{
    gatherEntries(numNodes, index, 0, in->getNumNodes(), Id, in->Id,
            Tag, in->Tag, globalDegreesOfFreedom, in->globalDegreesOfFreedom,
            numDim, Coordinates, in->Coordinates);
    
}

void NodeFile::gather_global(const index_t* index, const NodeFile* in)
{
    // get the global range of node ids
    const std::pair<index_t,index_t> id_range(in->getGlobalIdRange());
    const index_t undefined_node = id_range.first-1;
    std::vector<index_t> distribution(in->MPIInfo->size+1);

    // distribute the range of node ids
    index_t buffer_len = in->MPIInfo->setDistribution(id_range.first, id_range.second, &distribution[0]);

    // allocate buffers
    index_t* Id_buffer = new index_t[buffer_len];
    int* Tag_buffer = new int[buffer_len];
    index_t* globalDegreesOfFreedom_buffer = new index_t[buffer_len];
    double* Coordinates_buffer = new double[buffer_len*numDim];

    // fill Id_buffer by the undefined_node marker to check if nodes
    // are defined
#pragma omp parallel for
    for (index_t n = 0; n < buffer_len; n++)
        Id_buffer[n] = undefined_node;

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest = in->MPIInfo->mod_rank(in->MPIInfo->rank+1);
    int source = in->MPIInfo->mod_rank(in->MPIInfo->rank-1);
#endif
    int buffer_rank = in->MPIInfo->rank;
    for (int p=0; p<in->MPIInfo->size; ++p) {
        if (p>0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_DIM_T, dest,
                    in->MPIInfo->counter(), source,
                    in->MPIInfo->counter(), in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->counter()+1, source,
                    in->MPIInfo->counter()+1, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                    MPI_DIM_T, dest, in->MPIInfo->counter()+2, source,
                    in->MPIInfo->counter()+2, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*numDim,
                    MPI_DOUBLE, dest, in->MPIInfo->counter()+3, source,
                    in->MPIInfo->counter()+3, in->MPIInfo->comm, &status);
	        in->MPIInfo->incCounter(4);
#endif
        }
        buffer_rank=in->MPIInfo->mod_rank(buffer_rank-1);
        scatterEntries(in->numNodes, in->Id, distribution[buffer_rank],
                distribution[buffer_rank+1], Id_buffer, in->Id,
                Tag_buffer, in->Tag, globalDegreesOfFreedom_buffer,
                in->globalDegreesOfFreedom, numDim, Coordinates_buffer,
                in->Coordinates);
    }
    // now entries are collected from the buffer again by sending the
    // entries around in a circle
#ifdef ESYS_MPI
    dest = in->MPIInfo->mod_rank(in->MPIInfo->rank+1);
    source = in->MPIInfo->mod_rank(in->MPIInfo->rank-1);
#endif
    buffer_rank=in->MPIInfo->rank;
    for (int p=0; p<in->MPIInfo->size; ++p) {
        gatherEntries(numNodes, index, distribution[buffer_rank],
                distribution[buffer_rank+1], Id, Id_buffer, Tag, Tag_buffer,
                globalDegreesOfFreedom, globalDegreesOfFreedom_buffer, numDim,
                Coordinates, Coordinates_buffer);
        if (p < in->MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_DIM_T, dest,
                    in->MPIInfo->counter(), source,
                    in->MPIInfo->counter(), in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->counter()+1, source,
                    in->MPIInfo->counter()+1, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                    MPI_DIM_T, dest, in->MPIInfo->counter()+2, source,
                    in->MPIInfo->counter()+2, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*numDim,
                    MPI_DOUBLE, dest, in->MPIInfo->counter()+3, source,
                    in->MPIInfo->counter()+3, in->MPIInfo->comm, &status);
            in->MPIInfo->incCounter(4);
#endif
        }
        buffer_rank=in->MPIInfo->mod_rank(buffer_rank-1);
    }
#if DOASSERT
    // check if all nodes are set:
    index_t err=-1;
#pragma omp parallel for
    for (index_t n=0; n<numNodes; ++n) {
        if (Id[n] == undefined_node) {
#pragma omp critical
            err=n;
        }
    }
    if (err>=0) {
        std::stringstream ss;
        ss << "NodeFile::gather_global: Node id " << Id[err]
            << " at position " << err << " is referenced but not defined.";
        const std::string errorMsg(ss.str());
        throw escript::AssertException(errorMsg);
    }
#endif // DOASSERT
    delete[] Id_buffer;
    delete[] Tag_buffer;
    delete[] globalDegreesOfFreedom_buffer;
    delete[] Coordinates_buffer;
}

void NodeFile::assignMPIRankToDOFs(std::vector<int>& mpiRankOfDOF,
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

dim_t NodeFile::prepareLabeling(const std::vector<short>& mask,
                                IndexVector& buffer,
                                IndexVector& distribution,
                                bool useNodes)
{
    const index_t UNSET_ID=-1,SET_ID=1;

    // get the global range of DOF/node ids
    std::pair<index_t,index_t> idRange(useNodes ?
            getGlobalNodeIDIndexRange() : getGlobalDOFRange());
    const index_t* indexArray = (useNodes ? globalNodesIndex : globalDegreesOfFreedom);
    // distribute the range of node ids
    distribution.assign(MPIInfo->size+1, 0);
    int buffer_len = MPIInfo->setDistribution(idRange.first,
            idRange.second, &distribution[0]);
    const dim_t myCount = distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];

    // fill buffer by the UNSET_ID marker to check if nodes are defined
    buffer.assign(buffer_len, UNSET_ID);

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        if (p>0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&buffer[0], buffer.size(), MPI_DIM_T, dest,
                    MPIInfo->counter(), source, MPIInfo->counter(),
                    MPIInfo->comm, &status);
            MPIInfo->incCounter();
#endif
        }
        buffer_rank = MPIInfo->mod_rank(buffer_rank-1);
        const index_t id0 = distribution[buffer_rank];
        const index_t id1 = distribution[buffer_rank+1];
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            if (mask.size() < numNodes || mask[n] > -1) {
                const index_t k = indexArray[n];
                if (id0 <= k && k < id1) {
                    buffer[k - id0] = SET_ID;
                }
            }
        }
    }
    // count the entries in the buffer
    // TODO: OMP parallel
    index_t myNewCount = 0;
    for (index_t n = 0; n < myCount; ++n) {
        if (buffer[n] == SET_ID) {
            buffer[n] = myNewCount;
            myNewCount++;
        }
    }
    return myNewCount;
}

dim_t NodeFile::createDenseDOFLabeling()
{
    std::vector<index_t> DOF_buffer;
    std::vector<index_t> distribution;
    std::vector<index_t> loc_offsets(MPIInfo->size);
    std::vector<index_t> offsets(MPIInfo->size);
    index_t new_numGlobalDOFs = 0;

    // retrieve the number of own DOFs and fill buffer
    loc_offsets[MPIInfo->rank] = prepareLabeling(std::vector<short>(),
            DOF_buffer, distribution, false);
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], MPIInfo->size, MPI_DIM_T,
                  MPI_SUM, MPIInfo->comm);
    for (int n=0; n<MPIInfo->size; ++n) {
        loc_offsets[n]=new_numGlobalDOFs;
        new_numGlobalDOFs+=offsets[n];
    }
#else
    new_numGlobalDOFs = loc_offsets[0];
    loc_offsets[0] = 0;
#endif

    const dim_t myDOFs = distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];
#pragma omp parallel for
    for (index_t n = 0; n < myDOFs; ++n)
        DOF_buffer[n] += loc_offsets[MPIInfo->rank];

    std::vector<unsigned char> set_new_DOF(numNodes, true);

    // now entries are collected from the buffer again by sending them around
    // in a circle
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank = MPIInfo->rank;
    for (int p = 0; p < MPIInfo->size; ++p) {
        const index_t dof0 = distribution[buffer_rank];
        const index_t dof1 = distribution[buffer_rank+1];
#pragma omp parallel for
        for (index_t n = 0; n < numNodes; n++) {
            const index_t k = globalDegreesOfFreedom[n];
            if (set_new_DOF[n] && dof0<=k && k<dof1) {
                globalDegreesOfFreedom[n]=DOF_buffer[k-dof0];
                set_new_DOF[n]=false;
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&DOF_buffer[0], DOF_buffer.size(), MPI_DIM_T,
                    dest, MPIInfo->counter(), source,
                    MPIInfo->counter(), MPIInfo->comm, &status);
	        MPIInfo->incCounter();
#endif
        }
        buffer_rank = MPIInfo->mod_rank(buffer_rank-1);
    }

    return new_numGlobalDOFs;
}

dim_t NodeFile::createDenseNodeLabeling(IndexVector& nodeDistribution,
                                        const IndexVector& dofDistribution)
{
    const index_t UNSET_ID=-1, SET_ID=1;
    const index_t myFirstDOF = dofDistribution[MPIInfo->rank];
    const index_t myLastDOF = dofDistribution[MPIInfo->rank+1];

    // find the range of node ids controlled by me
    index_t min_id = std::numeric_limits<index_t>::max();
    index_t max_id = std::numeric_limits<index_t>::min();
#pragma omp parallel
    {
        index_t loc_max_id = max_id;
        index_t loc_min_id = min_id;
#pragma omp for
        for (index_t n = 0; n < numNodes; n++) {
            const dim_t dof = globalDegreesOfFreedom[n];
            if (myFirstDOF <= dof && dof < myLastDOF) {
                loc_max_id = std::max(loc_max_id, Id[n]);
                loc_min_id = std::min(loc_min_id, Id[n]);
            }
        }
#pragma omp critical
        {
            max_id = std::max(loc_max_id, max_id);
            min_id = std::min(loc_min_id, min_id);
        }
    }
    index_t my_buffer_len = (max_id>=min_id ? max_id-min_id+1 : 0);
    index_t buffer_len;

#ifdef ESYS_MPI
    MPI_Allreduce(&my_buffer_len, &buffer_len, 1, MPI_DIM_T, MPI_MAX,
                  MPIInfo->comm);
#else
    buffer_len=my_buffer_len;
#endif

    const int header_len=2;
    std::vector<index_t> Node_buffer(buffer_len+header_len, UNSET_ID);
    // extra storage for these IDs
    Node_buffer[0]=min_id;
    Node_buffer[1]=max_id;

    // mark and count the nodes in use
#pragma omp parallel for
    for (index_t n = 0; n < numNodes; n++) {
        globalNodesIndex[n] = -1;
        const index_t dof = globalDegreesOfFreedom[n];
        if (myFirstDOF <= dof && dof < myLastDOF)
            Node_buffer[Id[n]-min_id+header_len] = SET_ID;
    }
    index_t myNewNumNodes = 0;
    for (index_t n = 0; n < my_buffer_len; n++) {
        if (Node_buffer[header_len+n] == SET_ID) {
            Node_buffer[header_len+n] = myNewNumNodes;
            myNewNumNodes++;
        }
    }
    // make the local number of nodes globally available
#ifdef ESYS_MPI
    MPI_Allgather(&myNewNumNodes, 1, MPI_DIM_T, &nodeDistribution[0], 1,
                  MPI_DIM_T, MPIInfo->comm);
#else
    nodeDistribution[0] = myNewNumNodes;
#endif

    dim_t globalNumNodes = 0;
    for (int p = 0; p < MPIInfo->size; ++p) {
        const dim_t itmp = nodeDistribution[p];
        nodeDistribution[p] = globalNumNodes;
        globalNumNodes += itmp;
    }
    nodeDistribution[MPIInfo->size] = globalNumNodes;

    // offset node buffer
#pragma omp parallel for
    for (index_t n = 0; n < my_buffer_len; n++)
        Node_buffer[n+header_len] += nodeDistribution[MPIInfo->rank];

    // now we send this buffer around to assign global node index
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const index_t nodeID_0 = Node_buffer[0];
        const index_t nodeID_1 = Node_buffer[1];
        const index_t dof0 = dofDistribution[buffer_rank];
        const index_t dof1 = dofDistribution[buffer_rank+1];
        if (nodeID_0 <= nodeID_1) {
#pragma omp parallel for
            for (index_t n = 0; n < numNodes; n++) {
                const index_t dof = globalDegreesOfFreedom[n];
                const index_t id = Id[n]-nodeID_0;
                if (dof0 <= dof && dof < dof1 && id>=0 && id<=nodeID_1-nodeID_0)
                    globalNodesIndex[n] = Node_buffer[id+header_len];
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&Node_buffer[0], Node_buffer.size(), MPI_DIM_T,
                    dest, MPIInfo->counter(), source,
                    MPIInfo->counter(), MPIInfo->comm, &status);
	        MPIInfo->incCounter();
#endif
        }
        buffer_rank = MPIInfo->mod_rank(buffer_rank-1);
    }
    return globalNumNodes;
}

dim_t NodeFile::createDenseReducedLabeling(const std::vector<short>& reducedMask,
                                           bool useNodes)
{
    std::vector<index_t> buffer;
    std::vector<index_t> distribution;
    std::vector<index_t> loc_offsets(MPIInfo->size);
    std::vector<index_t> offsets(MPIInfo->size);
    dim_t new_numGlobalReduced=0;

    // retrieve the number of own DOFs/nodes and fill buffer
    loc_offsets[MPIInfo->rank]=prepareLabeling(reducedMask, buffer,
                                               distribution, useNodes);
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], MPIInfo->size, MPI_DIM_T,
                  MPI_SUM, MPIInfo->comm);
    for (int n=0; n<MPIInfo->size; ++n) {
        loc_offsets[n]=new_numGlobalReduced;
        new_numGlobalReduced+=offsets[n];
    }
#else
    new_numGlobalReduced=loc_offsets[0];
    loc_offsets[0]=0;
#endif

    const dim_t myCount=distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];
#pragma omp parallel for
    for (index_t n=0; n<myCount; ++n)
        buffer[n]+=loc_offsets[MPIInfo->rank];

    const index_t* denseArray =
        (useNodes ? globalNodesIndex : globalDegreesOfFreedom);
    index_t* reducedArray =
        (useNodes ? globalReducedNodesIndex : globalReducedDOFIndex);

#pragma omp parallel for
    for (index_t n=0; n<numNodes; ++n)
        reducedArray[n]=loc_offsets[0]-1;

    // now entries are collected from the buffer by sending them around
    // in a circle
#ifdef ESYS_MPI
    int dest = MPIInfo->mod_rank(MPIInfo->rank + 1);
    int source = MPIInfo->mod_rank(MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const index_t id0=distribution[buffer_rank];
        const index_t id1=distribution[buffer_rank+1];
#pragma omp parallel for
        for (index_t n=0; n<numNodes; n++) {
            if (reducedMask[n] > -1) {
                const index_t k=denseArray[n];
                if (id0<=k && k<id1)
                    reducedArray[n]=buffer[k-id0];
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&buffer[0], buffer.size(), MPI_DIM_T, dest,
                    MPIInfo->counter(), source,
                    MPIInfo->counter(), MPIInfo->comm, &status);
	        MPIInfo->incCounter();
#endif
        }
        buffer_rank = MPIInfo->mod_rank(buffer_rank-1);
    }
    return new_numGlobalReduced;
}

void NodeFile::createDOFMappingAndCoupling(bool use_reduced_elements) 
{
    escript::Distribution_ptr dofDistribution;
    const index_t* globalDOFIndex;
    if (use_reduced_elements) {
        dofDistribution = reducedDegreesOfFreedomDistribution;
        globalDOFIndex = globalReducedDOFIndex;
    } else {
        dofDistribution = degreesOfFreedomDistribution;
        globalDOFIndex = globalDegreesOfFreedom;
    }
    NodeMapping& mapping = (use_reduced_elements ?
                     reducedDegreesOfFreedomMapping : degreesOfFreedomMapping);

    const index_t myFirstDOF = dofDistribution->getFirstComponent();
    const index_t myLastDOF = dofDistribution->getLastComponent();
    const int mpiSize = MPIInfo->size;
    const int myRank = MPIInfo->rank;

    index_t min_DOF, max_DOF;
    std::pair<index_t,index_t> DOF_range(util::getFlaggedMinMaxInt(
                                            numNodes, globalDOFIndex, -1));

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
        throw FinleyException(gmsg);
    }

    const index_t UNUSED = -1;
    const dim_t len_loc_dof = max_DOF - min_DOF + 1;
    std::vector<index_t> shared(numNodes * (p_max - p_min + 1));
    std::vector<index_t> locDOFMask(len_loc_dof, UNUSED);

#ifdef BOUNDS_CHECK
    ESYS_ASSERT(myLastDOF-min_DOF <= len_loc_dof, "BOUNDS_CHECK");
#endif

#pragma omp parallel
    {
#pragma omp for
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = globalDOFIndex[i];
            if (k > -1) {
#ifdef BOUNDS_CHECK
                ESYS_ASSERT(k - min_DOF < len_loc_dof, "BOUNDS_CHECK");
#endif
                locDOFMask[k - min_DOF] = UNUSED - 1;
            }
       }
#pragma omp for
        for (index_t i = myFirstDOF - min_DOF; i < myLastDOF - min_DOF; ++i) {
            locDOFMask[i] = i - myFirstDOF + min_DOF;
        }
    }

    std::vector<index_t> wanted_DOFs(numNodes);
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
            ESYS_ASSERT(lastDOF - min_DOF <= len_loc_dof, "BOUNDS_CHECK");
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
    std::vector<index_t> nodeMask(numNodes, UNUSED);
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i) {
        const index_t k = globalDOFIndex[i];
        if (k > -1)
            nodeMask[i] = locDOFMask[k - min_DOF];
    }

    // now we can set the mapping from nodes to local DOFs
    mapping.assign(nodeMask, UNUSED);

    // define how to get DOF values for controlled but other processors
#ifdef BOUNDS_CHECK
    ESYS_ASSERT(numNodes == 0 || offsetInShared.back() < numNodes * (p_max - p_min + 1), "BOUNDS_CHECK");
#endif
#pragma omp parallel for
    for (index_t i = 0; i < lastn; ++i)
        shared[i] = myLastDOF - myFirstDOF + i;

#ifdef ESYS_HAVE_PASO
    index_t* p = shared.empty() ? NULL : &shared[0];
    paso::SharedComponents_ptr rcv_shcomp(new paso::SharedComponents(
            myLastDOF - myFirstDOF, neighbour, p, offsetInShared));
#endif

    /////////////////////////////////
    //   now we build the sender   //
    /////////////////////////////////
#ifdef ESYS_MPI
    std::vector<MPI_Request> mpi_requests(mpiSize * 2);
    std::vector<MPI_Status> mpi_stati(mpiSize * 2);
    MPI_Alltoall(&rcv_len[0], 1, MPI_DIM_T, &snd_len[0], 1, MPI_DIM_T, MPIInfo->comm);
    int count = 0;
    for (int p = 0; p < neighbour.size(); p++) {
        MPI_Isend(&wanted_DOFs[offsetInShared[p]],
                offsetInShared[p+1] - offsetInShared[p],
                MPI_DIM_T, neighbour[p], MPIInfo->counter() + myRank,
                MPIInfo->comm, &mpi_requests[count]);
        count++;
    }
    n = 0;
    neighbour.clear();
    offsetInShared.clear();
    for (int p = 0; p < mpiSize; p++) {
        if (snd_len[p] > 0) {
            MPI_Irecv(&shared[n], snd_len[p], MPI_DIM_T, p,
                      MPIInfo->counter()+p, MPIInfo->comm,
                      &mpi_requests[count]);
            count++;
            neighbour.push_back(p);
            offsetInShared.push_back(n);
            n += snd_len[p];
        }
    }
    MPIInfo->incCounter(MPIInfo->size);
    MPI_Waitall(count, &mpi_requests[0], &mpi_stati[0]);
    offsetInShared.push_back(n);

    // map global IDs to local IDs
#pragma omp parallel for
    for (index_t i = 0; i < n; ++i) {
        shared[i] = locDOFMask[shared[i] - min_DOF];
    }
#endif // ESYS_MPI

#ifdef ESYS_HAVE_PASO
    paso::SharedComponents_ptr snd_shcomp(new paso::SharedComponents(
            myLastDOF - myFirstDOF, neighbour, p, offsetInShared));

    if (use_reduced_elements) {
        reducedDegreesOfFreedomConnector.reset(new paso::Connector(snd_shcomp, rcv_shcomp));
    } else {
        degreesOfFreedomConnector.reset(new paso::Connector(snd_shcomp, rcv_shcomp));
    }
#endif // ESYS_HAVE_PASO

#ifdef ESYS_HAVE_TRILINOS
    using namespace esys_trilinos;

    const dim_t myNumTargets = myLastDOF - myFirstDOF;
    const dim_t numTargets = mapping.getNumTargets();
    IndexVector myRows(myNumTargets);
    IndexVector columns(numTargets);
    const IndexVector& dofMap = mapping.map;

#pragma omp parallel
    {
#pragma omp for nowait
        for (size_t i = 0; i < myNumTargets; i++) {
            myRows[i] = globalDOFIndex[dofMap[i]];
        }
#pragma omp for
        for (size_t i = 0; i < numTargets; i++) {
            columns[i] = globalDOFIndex[dofMap[i]];
        }
    } // end parallel section

    const dim_t numTotal = dofDistribution->getGlobalNumComponents();
    if (use_reduced_elements) {
        trilinosReducedRowMap.reset(new MapType(numTotal, myRows, 0,
                                      TeuchosCommFromEsysComm(MPIInfo->comm)));
        trilinosReducedColMap.reset(new MapType(numTotal, columns, 0,
                                      TeuchosCommFromEsysComm(MPIInfo->comm)));
    } else {
        trilinosRowMap.reset(new MapType(numTotal, myRows, 0,
                                      TeuchosCommFromEsysComm(MPIInfo->comm)));
        trilinosColMap.reset(new MapType(numTotal, columns, 0,
                                      TeuchosCommFromEsysComm(MPIInfo->comm)));
    }
#endif // ESYS_HAVE_TRILINOS
}

void NodeFile::createNodeMappings(const IndexVector& indexReducedNodes,
                                  const IndexVector& dofDist,
                                  const IndexVector& nodeDist)
{
    const int mpiSize = MPIInfo->size;
    const int myRank = MPIInfo->rank;

    const index_t myFirstDOF = dofDist[myRank];
    const index_t myLastDOF = dofDist[myRank+1];
    const index_t myNumDOF = myLastDOF-myFirstDOF;

    const index_t myFirstNode = nodeDist[myRank];
    const index_t myLastNode = nodeDist[myRank+1];
    const index_t myNumNodes = myLastNode-myFirstNode;

    std::vector<short> maskMyReducedDOF(myNumDOF, -1);
    std::vector<short> maskMyReducedNodes(myNumNodes, -1);
    const index_t iRNsize = indexReducedNodes.size();

    // mark the nodes used by the reduced mesh
#pragma omp parallel for
    for (index_t i = 0; i < iRNsize; ++i) {
        index_t k = globalNodesIndex[indexReducedNodes[i]];
        if (k >= myFirstNode && myLastNode > k)
            maskMyReducedNodes[k - myFirstNode] = 1;
        k = globalDegreesOfFreedom[indexReducedNodes[i]];
        if (k >= myFirstDOF && myLastDOF > k) {
            maskMyReducedDOF[k - myFirstDOF] = 1;
        }
    }
    IndexVector indexMyReducedDOF = util::packMask(maskMyReducedDOF);
    index_t myNumReducedDOF = indexMyReducedDOF.size();
    IndexVector indexMyReducedNodes = util::packMask(maskMyReducedNodes);
    index_t myNumReducedNodes = indexMyReducedNodes.size();

    IndexVector rdofDist(mpiSize+1);
    IndexVector rnodeDist(mpiSize+1);
#ifdef ESYS_MPI
    MPI_Allgather(&myNumReducedNodes, 1, MPI_DIM_T, &rnodeDist[0], 1, MPI_DIM_T, MPIInfo->comm);
    MPI_Allgather(&myNumReducedDOF, 1, MPI_DIM_T, &rdofDist[0], 1, MPI_DIM_T, MPIInfo->comm);
#else
    rnodeDist[0] = myNumReducedNodes;
    rdofDist[0] = myNumReducedDOF;
#endif
    index_t globalNumReducedNodes = 0;
    index_t globalNumReducedDOF = 0;
    for (int i = 0; i < mpiSize; ++i) {
        index_t k = rnodeDist[i];
        rnodeDist[i] = globalNumReducedNodes;
        globalNumReducedNodes += k;

        k = rdofDist[i];
        rdofDist[i] = globalNumReducedDOF;
        globalNumReducedDOF += k;
    }
    rnodeDist[mpiSize] = globalNumReducedNodes;
    rdofDist[mpiSize] = globalNumReducedDOF;

    // ==== distribution of Nodes ====
    nodesDistribution.reset(new escript::Distribution(MPIInfo, nodeDist));

    // ==== distribution of DOFs ====
    degreesOfFreedomDistribution.reset(new escript::Distribution(MPIInfo, dofDist));

    // ==== distribution of reduced Nodes ====
    reducedNodesDistribution.reset(new escript::Distribution(MPIInfo, rnodeDist));

    // ==== distribution of reduced DOF ====
    reducedDegreesOfFreedomDistribution.reset(new escript::Distribution(
                                                MPIInfo, rdofDist));

    IndexVector nodeMask(numNodes);
    const index_t UNUSED = -1;

    // ==== nodes mapping (dummy) ====
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i)
        nodeMask[i] = i;
    nodesMapping.assign(nodeMask, UNUSED);

    // ==== mapping between nodes and reduced nodes ====
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; ++i)
        nodeMask[i] = UNUSED;
#pragma omp parallel for
    for (index_t i = 0; i < iRNsize; ++i)
        nodeMask[indexReducedNodes[i]] = i;
    reducedNodesMapping.assign(nodeMask, UNUSED);

    // ==== mapping between nodes and DOFs + DOF connector
    createDOFMappingAndCoupling(false);
    // ==== mapping between nodes and reduced DOFs + reduced DOF connector
    createDOFMappingAndCoupling(true);

    // get the Ids for DOFs and reduced nodes
    const index_t rnTargets = reducedNodesMapping.getNumTargets();
    const index_t dofTargets = degreesOfFreedomMapping.getNumTargets();
    const index_t rdofTargets = reducedDegreesOfFreedomMapping.getNumTargets();
#pragma omp parallel
    {
#pragma omp for nowait
        for (index_t i = 0; i < rnTargets; ++i)
            reducedNodesId[i] = Id[reducedNodesMapping.map[i]];
#pragma omp for nowait
        for (index_t i = 0; i < dofTargets; ++i)
            degreesOfFreedomId[i] = Id[degreesOfFreedomMapping.map[i]];
#pragma omp for
        for (index_t i = 0; i < rdofTargets; ++i)
            reducedDegreesOfFreedomId[i] = Id[reducedDegreesOfFreedomMapping.map[i]];
    }
}

} // namespace finley

