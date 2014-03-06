
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

  Finley: NodeFile

*****************************************************************************/

#include "NodeFile.h"
#include <escript/Data.h>
#include <paso/Coupler.h>

#include <limits>
#include <sstream>

namespace finley {

// helper function
static void scatterEntries(int n, int* index, int min_index, int max_index,
                           int* Id_out, int* Id_in, int* Tag_out, int* Tag_in,
                           int* globalDegreesOfFreedom_out,
                           int* globalDegreesOfFreedom_in,
                           int numDim, double* Coordinates_out,
                           double* Coordinates_in)
{
    const int range = max_index-min_index;
    const size_t numDim_size = numDim*sizeof(double);

#pragma omp parallel for
    for (int i=0; i<n; i++) {
        const int k=index[i]-min_index;
        if ((k>=0) && (k<range)) {
            Id_out[k]=Id_in[i];
            Tag_out[k]=Tag_in[i];
            globalDegreesOfFreedom_out[k]=globalDegreesOfFreedom_in[i];
            memcpy(&(Coordinates_out[INDEX2(0,k,numDim)]),
                    &(Coordinates_in[INDEX2(0,i,numDim)]), numDim_size);
        }
    }
}

// helper function
static void gatherEntries(int n, const int* index, int min_index, int max_index,
                          int* Id_out, int* Id_in, int* Tag_out, int* Tag_in,
                          int* globalDegreesOfFreedom_out,
                          int* globalDegreesOfFreedom_in,
                          int numDim, double* Coordinates_out,
                          double* Coordinates_in)
{
    const int range = max_index-min_index;
    const size_t numDim_size = numDim*sizeof(double);

#pragma omp parallel for
    for (int i=0; i<n; i++) {
        const int k=index[i]-min_index;
        if ((k>=0) && (k<range)) {
            Id_out[i]=Id_in[k];
            Tag_out[i]=Tag_in[k];
            globalDegreesOfFreedom_out[i]=globalDegreesOfFreedom_in[k];
            memcpy(&(Coordinates_out[INDEX2(0,i,numDim)]),
                    &(Coordinates_in[INDEX2(0,k,numDim)]), numDim_size);
        }
    }
}

/// constructor
/// use NodeFile::allocTable to allocate the node table (Id,Coordinates)
NodeFile::NodeFile(int nDim, Esys_MPIInfo *mpiInfo) :
    numNodes(0),
    numDim(nDim),
    Id(NULL),
    Tag(NULL),
    globalDegreesOfFreedom(NULL),
    Coordinates(NULL),
    globalReducedDOFIndex(NULL),
    globalReducedNodesIndex(NULL),
    globalNodesIndex(NULL),
    nodesDistribution(NULL),
    reducedNodesDistribution(NULL),
    degreesOfFreedomDistribution(NULL),
    reducedDegreesOfFreedomDistribution(NULL),
    degreesOfFreedomConnector(NULL),
    reducedDegreesOfFreedomConnector(NULL),
    reducedNodesId(NULL),
    degreesOfFreedomId(NULL),
    reducedDegreesOfFreedomId(NULL),
    status(FINLEY_INITIAL_STATUS)
{
    MPIInfo = Esys_MPIInfo_getReference(mpiInfo);
}

/// destructor
NodeFile::~NodeFile()
{
    freeTable();
    Esys_MPIInfo_free(MPIInfo);
}

/// allocates the node table within this node file to hold NN nodes.
void NodeFile::allocTable(int NN)
{
    if (numNodes>0)
        freeTable();

    Id=new int[NN];
    Coordinates=new double[NN*numDim];
    Tag=new int[NN];
    globalDegreesOfFreedom=new int[NN];
    globalReducedDOFIndex=new int[NN];
    globalReducedNodesIndex=new int[NN];
    globalNodesIndex=new int[NN];
    reducedNodesId=new int[NN];
    degreesOfFreedomId=new int[NN];
    reducedDegreesOfFreedomId=new int[NN];
    numNodes=NN;

    // this initialization makes sure that data are located on the right
    // processor
#pragma omp parallel for
    for (int n=0; n<numNodes; n++) {
        Id[n]=-1;
        for (int i=0; i<numDim; i++)
            Coordinates[INDEX2(i,n,numDim)]=0.;
        Tag[n]=-1;
        globalDegreesOfFreedom[n]=-1;
        globalReducedDOFIndex[n]=-1;
        globalReducedNodesIndex[n]=-1;
        globalNodesIndex[n]=-1;
        reducedNodesId[n]=-1;
        degreesOfFreedomId[n]=-1;
        reducedDegreesOfFreedomId[n]=-1;
    }
}

/// frees the node table within this node file
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
    Paso_Distribution_free(nodesDistribution);
    nodesDistribution=NULL;
    Paso_Distribution_free(reducedNodesDistribution);
    nodesDistribution=NULL;
    Paso_Distribution_free(degreesOfFreedomDistribution);
    degreesOfFreedomDistribution=NULL;
    Paso_Distribution_free(reducedDegreesOfFreedomDistribution);
    reducedDegreesOfFreedomDistribution=NULL;
    Paso_Connector_free(degreesOfFreedomConnector);
    degreesOfFreedomConnector=NULL;
    Paso_Connector_free(reducedDegreesOfFreedomConnector);
    reducedDegreesOfFreedomConnector=NULL;

    numNodes=0;
}

void NodeFile::print() const
{
    std::cout << "=== " << numDim << "D-Nodes:\nnumber of nodes=" << numNodes
        << std::endl;
    std::cout << "Id,Tag,globalDegreesOfFreedom,degreesOfFreedom,reducedDegreesOfFeedom,node,reducedNode,Coordinates" << std::endl;
    for (int i=0; i<numNodes; i++) {
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

/// copies the array newX into this->coordinates
void NodeFile::setCoordinates(const escript::Data& newX)
{
    if (newX.getDataPointSize() != numDim)  {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of dimensions of new "
            "coordinates has to be " << numDim;
        const std::string errorMsg(ss.str());
        setError(VALUE_ERROR, errorMsg.c_str());
    } else if (newX.getNumDataPointsPerSample() != 1 ||
            newX.getNumSamples() != numNodes) {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of given nodes must be "
            << numNodes;
        const std::string errorMsg(ss.str());
        setError(VALUE_ERROR, errorMsg.c_str());
    } else {
        const size_t numDim_size=numDim*sizeof(double);
        ++status;
#pragma omp parallel for
        for (int n=0; n<numNodes; n++) {
            memcpy(&(Coordinates[INDEX2(0,n,numDim)]), newX.getSampleDataRO(n), numDim_size);
        }
    }
}

/// sets tags to newTag where mask>0
void NodeFile::setTags(const int newTag, const escript::Data& mask)
{
    resetError();

    if (1 != mask.getDataPointSize()) {
       setError(TYPE_ERROR, "NodeFile::setTags: number of components of mask must be 1.");
       return;
    } else if (mask.getNumDataPointsPerSample() != 1 ||
            mask.getNumSamples() != numNodes) {
       setError(TYPE_ERROR, "NodeFile::setTags: illegal number of samples of mask Data object");
       return;
    }

#pragma omp parallel for
    for (int n=0; n<numNodes; n++) {
         if (mask.getSampleDataRO(n)[0] > 0)
             Tag[n]=newTag;
    }
    updateTagList();
}

std::pair<int,int> NodeFile::getDOFRange() const
{
    std::pair<int,int> result(util::getMinMaxInt(
                                        1, numNodes, globalDegreesOfFreedom));
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

std::pair<int,int> NodeFile::getGlobalIdRange() const
{
    std::pair<int,int> result(util::getMinMaxInt(1, numNodes, Id));

#ifdef ESYS_MPI
    int global_id_range[2];
    int id_range[2] = { -result.first, result.second };
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, MPIInfo->comm);
    result.first = -global_id_range[0];
    result.second = global_id_range[1];
#endif
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

std::pair<int,int> NodeFile::getGlobalDOFRange() const
{
    std::pair<int,int> result(util::getMinMaxInt(
                                        1, numNodes, globalDegreesOfFreedom));

#ifdef ESYS_MPI
    int global_id_range[2];
    int id_range[2] = { -result.first, result.second };
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, MPIInfo->comm);
    result.first = -global_id_range[0];
    result.second = global_id_range[1];
#endif
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

std::pair<int,int> NodeFile::getGlobalNodeIDIndexRange() const
{
    std::pair<int,int> result(util::getMinMaxInt(1, numNodes, globalNodesIndex));

#ifdef ESYS_MPI
    int global_id_range[2];
    int id_range[2] = { -result.first, result.second };
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, MPIInfo->comm);
    result.first = -global_id_range[0];
    result.second = global_id_range[1];
#endif
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

void NodeFile::copyTable(int offset, int idOffset, int dofOffset,
                         const NodeFile* in)
{
    // check number of dimensions and table size
    if (numDim != in->numDim) {
        setError(TYPE_ERROR, "NodeFile::copyTable: dimensions of node files don't match");
        return;
    }
    if (numNodes < in->numNodes+offset) {
        setError(MEMORY_ERROR, "NodeFile::copyTable: node table is too small.");
        return;
    }

#pragma omp parallel for
    for (int n=0; n<in->numNodes; n++) {
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
void NodeFile::scatter(int* index, const NodeFile* in)
{
    scatterEntries(numNodes, index, 0, in->numNodes, Id, in->Id, Tag, in->Tag,
                   globalDegreesOfFreedom, in->globalDegreesOfFreedom,
                   numDim, Coordinates, in->Coordinates);
}

/// gathers this NodeFile from the NodeFile 'in' using the entries in
/// index[0:out->numNodes-1] which are between min_index and max_index
/// (exclusive)
void NodeFile::gather(int* index, const NodeFile* in)
{
    const std::pair<int,int> id_range(in->getGlobalIdRange());
    gatherEntries(numNodes, index, id_range.first, id_range.second, Id, in->Id,
            Tag, in->Tag, globalDegreesOfFreedom, in->globalDegreesOfFreedom,
            numDim, Coordinates, in->Coordinates);
}

void NodeFile::gather_global(const std::vector<int>& index, const NodeFile* in)
{
    // get the global range of node ids
    const std::pair<int,int> id_range(in->getGlobalIdRange());
    const int undefined_node=id_range.first-1;
    std::vector<int> distribution(in->MPIInfo->size+1);

    // distribute the range of node ids
    int buffer_len=Esys_MPIInfo_setDistribution(in->MPIInfo,
            id_range.first, id_range.second, &distribution[0]);

    // allocate buffers
    int *Id_buffer=new int[buffer_len];
    int *Tag_buffer=new int[buffer_len];
    int *globalDegreesOfFreedom_buffer=new int[buffer_len];
    double *Coordinates_buffer=new double[buffer_len*numDim];

    // fill Id_buffer by the undefined_node marker to check if nodes
    // are defined
#pragma omp parallel for
    for (int n=0; n<buffer_len; n++)
        Id_buffer[n]=undefined_node;

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest=Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank+1);
    int source=Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank-1);
#endif
    int buffer_rank=in->MPIInfo->rank;
    for (int p=0; p<in->MPIInfo->size; ++p) {
        if (p>0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->msg_tag_counter, source,
                    in->MPIInfo->msg_tag_counter, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->msg_tag_counter+1, source,
                    in->MPIInfo->msg_tag_counter+1, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                    MPI_INT, dest, in->MPIInfo->msg_tag_counter+2, source,
                    in->MPIInfo->msg_tag_counter+2, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*numDim,
                    MPI_DOUBLE, dest, in->MPIInfo->msg_tag_counter+3, source,
                    in->MPIInfo->msg_tag_counter+3, in->MPIInfo->comm, &status);
#endif
	    ESYS_MPI_INC_COUNTER(*(in->MPIInfo), 4)
        }
        buffer_rank=Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
        scatterEntries(in->numNodes, in->Id, distribution[buffer_rank],
                distribution[buffer_rank+1], Id_buffer, in->Id,
                Tag_buffer, in->Tag, globalDegreesOfFreedom_buffer,
                in->globalDegreesOfFreedom, numDim, Coordinates_buffer,
                in->Coordinates);
    }
    // now entries are collected from the buffer again by sending the
    // entries around in a circle
#ifdef ESYS_MPI
    dest=Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank+1);
    source=Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank-1);
#endif
    buffer_rank=in->MPIInfo->rank;
    for (int p=0; p<in->MPIInfo->size; ++p) {
        gatherEntries(numNodes, &index[0], distribution[buffer_rank],
                distribution[buffer_rank+1], Id, Id_buffer, Tag, Tag_buffer,
                globalDegreesOfFreedom, globalDegreesOfFreedom_buffer, numDim,
                Coordinates, Coordinates_buffer);
        if (p < in->MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->msg_tag_counter, source,
                    in->MPIInfo->msg_tag_counter, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT, dest,
                    in->MPIInfo->msg_tag_counter+1, source,
                    in->MPIInfo->msg_tag_counter+1, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len,
                    MPI_INT, dest, in->MPIInfo->msg_tag_counter+2, source,
                    in->MPIInfo->msg_tag_counter+2, in->MPIInfo->comm, &status);
            MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*numDim,
                    MPI_DOUBLE, dest, in->MPIInfo->msg_tag_counter+3, source,
                    in->MPIInfo->msg_tag_counter+3, in->MPIInfo->comm, &status);
#endif
            ESYS_MPI_INC_COUNTER(*(in->MPIInfo), 4)
        }
        buffer_rank=Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
    }
    // check if all nodes are set:
#pragma omp parallel for
    for (int n=0; n<numNodes; ++n) {
        if (Id[n] == undefined_node) {
            std::stringstream ss;
            ss << "NodeFile::gather_global: Node id " << Id[n]
                << " at position " << n << " is referenced but not defined.";
            const std::string errorMsg(ss.str());
            setError(VALUE_ERROR, errorMsg.c_str());
        }
    }
    delete[] Id_buffer;
    delete[] Tag_buffer;
    delete[] globalDegreesOfFreedom_buffer;
    delete[] Coordinates_buffer;
    // make sure that the error is global
    Esys_MPIInfo_noError(in->MPIInfo);
}

void NodeFile::assignMPIRankToDOFs(std::vector<int>& mpiRankOfDOF,
                                   const std::vector<int>& distribution)
{
    Esys_MPI_rank p_min=MPIInfo->size, p_max=-1;
    // first we retrieve the min and max DOF on this processor to reduce
    // costs for searching
    const std::pair<int,int> dof_range(getDOFRange());

    for (int p=0; p<MPIInfo->size; ++p) {
        if (distribution[p]<=dof_range.first) p_min=p;
        if (distribution[p]<=dof_range.second) p_max=p;
    }
#pragma omp parallel for
    for (int n=0; n<numNodes; ++n) {
        const int k=globalDegreesOfFreedom[n];
        for (int p=p_min; p<=p_max; ++p) {
            if (k < distribution[p+1]) {
                mpiRankOfDOF[n]=p;
                break;
            }
        }
    }
}

int NodeFile::prepareLabeling(const std::vector<short>& mask,
                              std::vector<int>& buffer,
                              std::vector<int>& distribution, bool useNodes)
{
    const int UNSET_ID=-1,SET_ID=1;

    // get the global range of DOF/node ids
    std::pair<int,int> idRange(useNodes ?
            getGlobalNodeIDIndexRange() : getGlobalDOFRange());
    const int* indexArray = (useNodes ? globalNodesIndex : globalDegreesOfFreedom);
    // distribute the range of node ids
    distribution.assign(MPIInfo->size+1, 0);
    int buffer_len=Esys_MPIInfo_setDistribution(MPIInfo, idRange.first,
            idRange.second, &distribution[0]);
    const int myCount=distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];

    // fill buffer by the UNSET_ID marker to check if nodes are defined
    buffer.assign(buffer_len, UNSET_ID);

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    int dest=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank + 1);
    int source=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        if (p>0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&buffer[0], buffer.size(), MPI_INT, dest,
                    MPIInfo->msg_tag_counter, source, MPIInfo->msg_tag_counter,
                    MPIInfo->comm, &status);
#endif
            MPIInfo->msg_tag_counter++;
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
        const int id0=distribution[buffer_rank];
        const int id1=distribution[buffer_rank+1];
#pragma omp parallel for
        for (int n=0; n<numNodes; n++) {
            if (mask.size()<numNodes || mask[n]>-1) {
                const int k=indexArray[n];
                if (id0<=k && k<id1) {
                    buffer[k-id0] = SET_ID;
                }
            }
        }
    }
    // count the entries in the buffer
    // TODO: OMP parallel
    int myNewCount=0;
    for (int n=0; n<myCount; ++n) {
        if (buffer[n] == SET_ID) {
            buffer[n]=myNewCount;
            myNewCount++;
        }
    }
    return myNewCount;
}

int NodeFile::createDenseDOFLabeling()
{
    std::vector<int> DOF_buffer;
    std::vector<int> distribution;
    std::vector<int> loc_offsets(MPIInfo->size);
    std::vector<int> offsets(MPIInfo->size);
    int new_numGlobalDOFs=0;

    // retrieve the number of own DOFs and fill buffer
    loc_offsets[MPIInfo->rank]=prepareLabeling(std::vector<short>(),
            DOF_buffer, distribution, false);
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], MPIInfo->size, MPI_INT,
                  MPI_SUM, MPIInfo->comm);
    for (int n=0; n<MPIInfo->size; ++n) {
        loc_offsets[n]=new_numGlobalDOFs;
        new_numGlobalDOFs+=offsets[n];
    }
#else
    new_numGlobalDOFs=loc_offsets[0];
    loc_offsets[0]=0;
#endif

    const int myDOFs=distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];
#pragma omp parallel for
    for (int n=0; n<myDOFs; ++n)
        DOF_buffer[n]+=loc_offsets[MPIInfo->rank];

    std::vector<unsigned char> set_new_DOF(numNodes, true);

    // now entries are collected from the buffer again by sending them around
    // in a circle
#ifdef ESYS_MPI
    int dest=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank + 1);
    int source=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int dof0=distribution[buffer_rank];
        const int dof1=distribution[buffer_rank+1];
#pragma omp parallel for
        for (int n=0; n<numNodes; n++) {
            const int k=globalDegreesOfFreedom[n];
            if (set_new_DOF[n] && dof0<=k && k<dof1) {
                globalDegreesOfFreedom[n]=DOF_buffer[k-dof0];
                set_new_DOF[n]=false;
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&DOF_buffer[0], DOF_buffer.size(), MPI_INT,
                    dest, MPIInfo->msg_tag_counter, source,
                    MPIInfo->msg_tag_counter, MPIInfo->comm, &status);
#endif
	    ESYS_MPI_INC_COUNTER(*MPIInfo, 1)
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }

    return new_numGlobalDOFs;
}

int NodeFile::createDenseNodeLabeling(std::vector<int>& nodeDistribution,
                                      const std::vector<int>& dofDistribution)
{
    const int UNSET_ID=-1, SET_ID=1;
    const int myFirstDOF=dofDistribution[MPIInfo->rank];
    const int myLastDOF=dofDistribution[MPIInfo->rank+1];

    // find the range of node ids controlled by me
    int min_id=std::numeric_limits<int>::max();
    int max_id=std::numeric_limits<int>::min();
#pragma omp parallel
    {
        int loc_max_id=max_id;
        int loc_min_id=min_id;
#pragma omp for
        for (int n=0; n<numNodes; n++) {
            const int dof=globalDegreesOfFreedom[n];
            if (myFirstDOF<=dof && dof<myLastDOF) {
                loc_max_id=std::max(loc_max_id, Id[n]);
                loc_min_id=std::min(loc_min_id, Id[n]);
            }
        }
#pragma omp critical
        {
            max_id=std::max(loc_max_id, max_id);
            min_id=std::min(loc_min_id, min_id);
        }
    }
    int my_buffer_len = (max_id>=min_id ? max_id-min_id+1 : 0);
    int buffer_len;

#ifdef ESYS_MPI
    MPI_Allreduce(&my_buffer_len, &buffer_len, 1, MPI_INT, MPI_MAX,
                  MPIInfo->comm);
#else
    buffer_len=my_buffer_len;
#endif

    const int header_len=2;
    std::vector<int> Node_buffer(buffer_len+header_len, UNSET_ID);
    // extra storage for these IDs
    Node_buffer[0]=min_id;
    Node_buffer[1]=max_id;

    // mark and count the nodes in use
#pragma omp parallel for
    for (int n=0; n<numNodes; n++) {
        globalNodesIndex[n]=-1;
        const int dof=globalDegreesOfFreedom[n];
        if (myFirstDOF<=dof && dof<myLastDOF)
            Node_buffer[Id[n]-min_id+header_len]=SET_ID;
    }
    int myNewNumNodes=0;
    for (int n=0; n<my_buffer_len; n++) {
        if (Node_buffer[header_len+n]==SET_ID) {
            Node_buffer[header_len+n]=myNewNumNodes;
            myNewNumNodes++;
        }
    }
    // make the local number of nodes globally available
#ifdef ESYS_MPI
    MPI_Allgather(&myNewNumNodes, 1, MPI_INT, &nodeDistribution[0], 1, MPI_INT,
                  MPIInfo->comm);
#else
    nodeDistribution[0]=myNewNumNodes;
#endif

    int globalNumNodes=0;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int itmp=nodeDistribution[p];
        nodeDistribution[p]=globalNumNodes;
        globalNumNodes+=itmp;
    }
    nodeDistribution[MPIInfo->size]=globalNumNodes;

    // offset node buffer
#pragma omp parallel for
    for (int n=0; n<my_buffer_len; n++)
        Node_buffer[n+header_len]+=nodeDistribution[MPIInfo->rank];

    // now we send this buffer around to assign global node index
#ifdef ESYS_MPI
    int dest=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank + 1);
    int source=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int nodeID_0=Node_buffer[0];
        const int nodeID_1=Node_buffer[1];
        const int dof0=dofDistribution[buffer_rank];
        const int dof1=dofDistribution[buffer_rank+1];
        if (nodeID_0 <= nodeID_1) {
#pragma omp parallel for
            for (int n=0; n<numNodes; n++) {
                const int dof=globalDegreesOfFreedom[n];
                const int id=Id[n]-nodeID_0;
                if (dof0<=dof && dof<dof1 && id>=0 && id<=nodeID_1-nodeID_0)
                    globalNodesIndex[n]=Node_buffer[id+header_len];
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&Node_buffer[0], Node_buffer.size(), MPI_INT,
                    dest, MPIInfo->msg_tag_counter, source,
                    MPIInfo->msg_tag_counter, MPIInfo->comm, &status);
#endif
	    ESYS_MPI_INC_COUNTER(*MPIInfo, 1)
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }
    return globalNumNodes;
}

int NodeFile::createDenseReducedLabeling(const std::vector<short>& reducedMask,
                                         bool useNodes)
{
    std::vector<int> buffer;
    std::vector<int> distribution;
    std::vector<int> loc_offsets(MPIInfo->size);
    std::vector<int> offsets(MPIInfo->size);
    int new_numGlobalReduced=0;

    // retrieve the number of own DOFs/nodes and fill buffer
    loc_offsets[MPIInfo->rank]=prepareLabeling(reducedMask, buffer,
                                               distribution, useNodes);
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], MPIInfo->size, MPI_INT,
                  MPI_SUM, MPIInfo->comm);
    for (int n=0; n<MPIInfo->size; ++n) {
        loc_offsets[n]=new_numGlobalReduced;
        new_numGlobalReduced+=offsets[n];
    }
#else
    new_numGlobalReduced=loc_offsets[0];
    loc_offsets[0]=0;
#endif

    const int myCount=distribution[MPIInfo->rank+1]-distribution[MPIInfo->rank];
#pragma omp parallel for
    for (int n=0; n<myCount; ++n)
        buffer[n]+=loc_offsets[MPIInfo->rank];

    const int* denseArray =
        (useNodes ? globalNodesIndex : globalDegreesOfFreedom);
    int* reducedArray =
        (useNodes ? globalReducedNodesIndex : globalReducedDOFIndex);

#pragma omp parallel for
    for (int n=0; n<numNodes; ++n)
        reducedArray[n]=loc_offsets[0]-1;

    // now entries are collected from the buffer by sending them around
    // in a circle
#ifdef ESYS_MPI
    int dest=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank + 1);
    int source=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int id0=distribution[buffer_rank];
        const int id1=distribution[buffer_rank+1];
#pragma omp parallel for
        for (int n=0; n<numNodes; n++) {
            if (reducedMask[n] > -1) {
                const int k=denseArray[n];
                if (id0<=k && k<id1)
                    reducedArray[n]=buffer[k-id0];
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&buffer[0], buffer.size(), MPI_INT, dest,
                    MPIInfo->msg_tag_counter, source,
                    MPIInfo->msg_tag_counter, MPIInfo->comm, &status);
#endif
	    ESYS_MPI_INC_COUNTER(*MPIInfo, 1)
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }
    return new_numGlobalReduced;
}

void NodeFile::createDOFMappingAndCoupling(bool use_reduced_elements) 
{
    Paso_Distribution* dof_distribution;
    const int* globalDOFIndex;
    if (use_reduced_elements) {
        dof_distribution=reducedDegreesOfFreedomDistribution;
        globalDOFIndex=globalReducedDOFIndex;
    } else {
        dof_distribution=degreesOfFreedomDistribution;
        globalDOFIndex=globalDegreesOfFreedom;
    }
    const int myFirstDOF=Paso_Distribution_getFirstComponent(dof_distribution);
    const int myLastDOF=Paso_Distribution_getLastComponent(dof_distribution);
    const int mpiSize=MPIInfo->size;
    const int myRank=MPIInfo->rank;

    int min_DOF, max_DOF;
    std::pair<int,int> DOF_range(util::getFlaggedMinMaxInt(
                                            numNodes, globalDOFIndex, -1));

    if (DOF_range.second < DOF_range.first) {
        min_DOF=myFirstDOF;
        max_DOF=myLastDOF-1;
    } else {
        min_DOF=DOF_range.first;
        max_DOF=DOF_range.second;
    }

    int p_min=mpiSize;
    int p_max=-1;
    if (max_DOF >= min_DOF) {
        for (int p=0; p<mpiSize; ++p) {
            if (dof_distribution->first_component[p]<=min_DOF) p_min=p;
            if (dof_distribution->first_component[p]<=max_DOF) p_max=p;
        }
    }

    if (!((min_DOF<=myFirstDOF) && (myLastDOF-1<=max_DOF))) {
        setError(SYSTEM_ERROR, "Local elements do not span local degrees of freedom.");
        return;
    }
    const int UNUSED = -1;
    const int len_loc_dof=max_DOF-min_DOF+1;
    std::vector<int> shared(numNodes*(p_max-p_min+1));
    std::vector<int> offsetInShared(mpiSize+1);
    std::vector<int> locDOFMask(len_loc_dof, UNUSED);

#pragma omp parallel 
    {
#pragma omp for
        for (int i=0;i<numNodes;++i) {
            const int k=globalDOFIndex[i];
            if (k > -1) {
#ifdef BOUNDS_CHECK
                if ((k-min_DOF)>=len_loc_dof) {
                    printf("BOUNDS_CHECK %s %d i=%d k=%d min_DOF=%d\n", __FILE__, __LINE__, i, k, min_DOF);
                    exit(1);
                }
#endif
                locDOFMask[k-min_DOF]=UNUSED-1;
            }
       }
#ifdef BOUNDS_CHECK
       if (myLastDOF-min_DOF > len_loc_dof) {
           printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
           exit(1);
       }
#endif
#pragma omp for
       for (int i=myFirstDOF-min_DOF; i<myLastDOF-min_DOF; ++i) {
            locDOFMask[i]=i-myFirstDOF+min_DOF;
        }
    }

    std::vector<int> wanted_DOFs(numNodes);
    std::vector<int> rcv_len(mpiSize);
    std::vector<int> snd_len(mpiSize);
    std::vector<int> neighbor(mpiSize);
    int numNeighbors=0;
    int n=0;
    int lastn=n;
    for (int p=p_min; p<=p_max; ++p) {
        if (p != myRank) {
            const int firstDOF=std::max(min_DOF, dof_distribution->first_component[p]);
            const int lastDOF=std::min(max_DOF+1, dof_distribution->first_component[p+1]);
#ifdef BOUNDS_CHECK
            if (firstDOF-min_DOF<0 || lastDOF-min_DOF>len_loc_dof) {
                printf("BOUNDS_CHECK %s %d p=%d\n", __FILE__, __LINE__, p);
                exit(1);
            }
#endif
            for (int i=firstDOF-min_DOF; i<lastDOF-min_DOF; ++i) {
                if (locDOFMask[i] == UNUSED-1) {
                   locDOFMask[i]=myLastDOF-myFirstDOF+n;
                   wanted_DOFs[n]=i+min_DOF;
                   ++n;
                }
            }
            if (n > lastn) {
                rcv_len[p]=n-lastn;
#ifdef BOUNDS_CHECK
                if (numNeighbors >= mpiSize+1) {
                    printf("BOUNDS_CHECK %s %d p=%d numNeighbors=%d n=%d\n", __FILE__, __LINE__, p, numNeighbors, n);
                    exit(1);
                }
#endif
                neighbor[numNeighbors]=p;
                offsetInShared[numNeighbors]=lastn;
                numNeighbors++;
                lastn=n;
            }
        } // if p!=myRank
    } // for p

#ifdef BOUNDS_CHECK
    if (numNeighbors >= mpiSize+1) {
        printf("BOUNDS_CHECK %s %d numNeighbors=%d\n", __FILE__, __LINE__, numNeighbors);
        exit(1);
    }
#endif
    offsetInShared[numNeighbors]=lastn;

    // assign new DOF labels to nodes
    std::vector<int> nodeMask(numNodes, UNUSED);
#pragma omp parallel for
    for (int i=0; i<numNodes; ++i) {
        const int k=globalDOFIndex[i];
        if (k > -1)
            nodeMask[i]=locDOFMask[k-min_DOF];
    }

    // now we can set the mapping from nodes to local DOFs
    if (use_reduced_elements) {
        reducedDegreesOfFreedomMapping.assign(nodeMask, UNUSED);
    } else {
        degreesOfFreedomMapping.assign(nodeMask, UNUSED);
    }

    // define how to get DOF values for controlled but other processors
#ifdef BOUNDS_CHECK
    if (offsetInShared[numNeighbors] >= numNodes*(p_max-p_min+1)) {
        printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
        exit(1);
    }
#endif
#pragma omp parallel for
    for (int i=0; i<offsetInShared[numNeighbors]; ++i)
        shared[i]=myLastDOF-myFirstDOF+i;

    Paso_SharedComponents *rcv_shcomp=Paso_SharedComponents_alloc(
            myLastDOF-myFirstDOF, numNeighbors, &neighbor[0], &shared[0],
            &offsetInShared[0], 1, 0, MPIInfo);

    /////////////////////////////////
    //   now we build the sender   //
    /////////////////////////////////
#ifdef ESYS_MPI
    std::vector<MPI_Request> mpi_requests(mpiSize*2);
    std::vector<MPI_Status> mpi_stati(mpiSize*2);
    MPI_Alltoall(&rcv_len[0], 1, MPI_INT, &snd_len[0], 1, MPI_INT, MPIInfo->comm);
    int count=0;
#else
    snd_len[0]=rcv_len[0];
#endif

    for (int p=0; p<rcv_shcomp->numNeighbors; p++) {
#ifdef ESYS_MPI
        MPI_Isend(&(wanted_DOFs[rcv_shcomp->offsetInShared[p]]),
                rcv_shcomp->offsetInShared[p+1]-rcv_shcomp->offsetInShared[p],
                MPI_INT, rcv_shcomp->neighbor[p],
                MPIInfo->msg_tag_counter+myRank, MPIInfo->comm,
                &mpi_requests[count]);
        count++;
#endif
    }
    n=0;
    numNeighbors=0;
    for (int p=0; p<mpiSize; p++) {
        if (snd_len[p] > 0) {
#ifdef ESYS_MPI
            MPI_Irecv(&shared[n], snd_len[p], MPI_INT, p,
                    MPIInfo->msg_tag_counter+p, MPIInfo->comm,
                    &mpi_requests[count]);
            count++;
#endif
            neighbor[numNeighbors]=p;
            offsetInShared[numNeighbors]=n;
            numNeighbors++;
            n+=snd_len[p];
        }
    }
    ESYS_MPI_INC_COUNTER(*MPIInfo, MPIInfo->size)
    offsetInShared[numNeighbors]=n;
#ifdef ESYS_MPI
    MPI_Waitall(count, &mpi_requests[0], &mpi_stati[0]);
#endif
    // map global ids to local id's
#pragma omp parallel for
    for (int i=0; i<offsetInShared[numNeighbors]; ++i) {
        shared[i]=locDOFMask[shared[i]-min_DOF];
    }

    Paso_SharedComponents* snd_shcomp=Paso_SharedComponents_alloc(
            myLastDOF-myFirstDOF, numNeighbors, &neighbor[0], &shared[0],
            &offsetInShared[0], 1, 0, MPIInfo);

    if (noError()) {
        if (use_reduced_elements) {
            reducedDegreesOfFreedomConnector=Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
        } else {
            degreesOfFreedomConnector=Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
        }
    }

    Paso_SharedComponents_free(rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);
}

void NodeFile::createNodeMappings(const std::vector<int>& indexReducedNodes,
                                  const std::vector<int>& dofDist,
                                  const std::vector<int>& nodeDist)
{
    const int mpiSize=MPIInfo->size;
    const int myRank=MPIInfo->rank;

    const int myFirstDOF=dofDist[myRank];
    const int myLastDOF=dofDist[myRank+1];
    const int myNumDOF=myLastDOF-myFirstDOF;

    const int myFirstNode=nodeDist[myRank];
    const int myLastNode=nodeDist[myRank+1];
    const int myNumNodes=myLastNode-myFirstNode;

    std::vector<short> maskMyReducedDOF(myNumDOF, -1);
    std::vector<short> maskMyReducedNodes(myNumNodes, -1);

    // mark the nodes used by the reduced mesh
#pragma omp parallel for
    for (int i=0; i<indexReducedNodes.size(); ++i) {
        int k=globalNodesIndex[indexReducedNodes[i]];
        if (k>=myFirstNode && myLastNode>k)
            maskMyReducedNodes[k-myFirstNode]=1;
        k=globalDegreesOfFreedom[indexReducedNodes[i]];
        if (k>=myFirstDOF && myLastDOF>k) {
            maskMyReducedDOF[k-myFirstDOF]=1;
        }
    }
    std::vector<int> indexMyReducedDOF = util::packMask(maskMyReducedDOF);
    int myNumReducedDOF=indexMyReducedDOF.size();
    std::vector<int> indexMyReducedNodes = util::packMask(maskMyReducedNodes);
    int myNumReducedNodes=indexMyReducedNodes.size();

    std::vector<int> rdofDist(mpiSize+1);
    std::vector<int> rnodeDist(mpiSize+1);
#ifdef ESYS_MPI
    MPI_Allgather(&myNumReducedNodes, 1, MPI_INT, &rnodeDist[0], 1, MPI_INT, MPIInfo->comm);
    MPI_Allgather(&myNumReducedDOF, 1, MPI_INT, &rdofDist[0], 1, MPI_INT, MPIInfo->comm);
#else
    rnodeDist[0]=myNumReducedNodes;
    rdofDist[0]=myNumReducedDOF;
#endif
    int globalNumReducedNodes=0;
    int globalNumReducedDOF=0;
    for (int i=0; i<mpiSize;++i) {
        int k=rnodeDist[i];
        rnodeDist[i]=globalNumReducedNodes;
        globalNumReducedNodes+=k;

        k=rdofDist[i];
        rdofDist[i]=globalNumReducedDOF;
        globalNumReducedDOF+=k;
    }
    rnodeDist[mpiSize]=globalNumReducedNodes;
    rdofDist[mpiSize]=globalNumReducedDOF;

    // ==== distribution of Nodes ===============================
    nodesDistribution=Paso_Distribution_alloc(MPIInfo, &nodeDist[0], 1, 0);
    // ==== distribution of DOFs ================================
    degreesOfFreedomDistribution=Paso_Distribution_alloc(MPIInfo, &dofDist[0], 1,0);
    // ==== distribution of reduced Nodes =======================
    reducedNodesDistribution=Paso_Distribution_alloc(MPIInfo, &rnodeDist[0], 1, 0);
    // ==== distribution of reduced DOF =========================
    reducedDegreesOfFreedomDistribution=Paso_Distribution_alloc(MPIInfo, &rdofDist[0], 1, 0);

    std::vector<int> nodeMask(numNodes);

    if (noError()) {
        const int UNUSED = -1;
        // ==== nodes mapping which is a dummy structure ========
#pragma omp parallel for
        for (int i=0; i<numNodes; ++i)
            nodeMask[i]=i;
        nodesMapping.assign(nodeMask, UNUSED);

        // ==== mapping between nodes and reduced nodes ==========
#pragma omp parallel for
        for (int i=0; i<numNodes; ++i)
            nodeMask[i]=UNUSED;
#pragma omp parallel for
        for (int i=0; i<indexReducedNodes.size(); ++i)
            nodeMask[indexReducedNodes[i]]=i;
        reducedNodesMapping.assign(nodeMask, UNUSED);
    }
    // ==== mapping between nodes and DOFs + DOF connector
    if (noError())
        createDOFMappingAndCoupling(false);
    // ==== mapping between nodes and reduced DOFs + reduced DOF connector
    if (noError())
        createDOFMappingAndCoupling(true);

    // get the Ids for DOFs and reduced nodes
    if (noError()) {
#pragma omp parallel
        {
#pragma omp for
         for (int i=0; i<reducedNodesMapping.getNumTargets(); ++i)
             reducedNodesId[i]=Id[reducedNodesMapping.map[i]];
#pragma omp for
         for (int i=0; i<degreesOfFreedomMapping.getNumTargets(); ++i)
             degreesOfFreedomId[i]=Id[degreesOfFreedomMapping.map[i]];
#pragma omp for
         for (int i=0; i<reducedDegreesOfFreedomMapping.getNumTargets(); ++i)
             reducedDegreesOfFreedomId[i]=Id[reducedDegreesOfFreedomMapping.map[i]];
        }
    } else {
        Paso_Distribution_free(nodesDistribution);
        Paso_Distribution_free(reducedNodesDistribution);
        Paso_Distribution_free(degreesOfFreedomDistribution);
        Paso_Distribution_free(reducedDegreesOfFreedomDistribution);
        Paso_Connector_free(degreesOfFreedomConnector);
        Paso_Connector_free(reducedDegreesOfFreedomConnector);
        nodesDistribution=NULL;
        reducedNodesDistribution=NULL;
        degreesOfFreedomDistribution=NULL;
        reducedDegreesOfFreedomDistribution=NULL;
        degreesOfFreedomConnector=NULL;
        reducedDegreesOfFreedomConnector=NULL;
    }
}

} // namespace finley

