
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Finley: NodeFile

*****************************************************************************/

#include "NodeFile.h"
#include "Util.h"
#include <escript/Data.h>

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
static void gatherEntries(int n, int* index, int min_index, int max_index,
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
    nodesMapping(NULL),
    reducedNodesMapping(NULL),
    degreesOfFreedomMapping(NULL),
    reducedDegreesOfFreedomMapping(NULL),
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
    Finley_NodeMapping_free(nodesMapping);
    nodesMapping=NULL;
    Finley_NodeMapping_free(reducedNodesMapping);
    reducedNodesMapping=NULL;
    Finley_NodeMapping_free(degreesOfFreedomMapping);
    degreesOfFreedomMapping=NULL;
    Finley_NodeMapping_free(reducedDegreesOfFreedomMapping);
    reducedDegreesOfFreedomMapping=NULL;
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

void NodeFile::updateTagList()
{
    Finley_Util_setValuesInUse(Tag, numNodes, tagsInUse, MPIInfo);
}

/// copies the array newX into this->coordinates
void NodeFile::setCoordinates(const escript::Data& cNewX)
{
    if (cNewX.getDataPointSize() != numDim)  {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of dimensions of new "
            "coordinates has to be " << numDim;
        const std::string errorMsg(ss.str());
        Finley_setError(VALUE_ERROR, errorMsg.c_str());
    } else if (cNewX.getNumDataPointsPerSample() != 1 ||
            cNewX.getNumSamples() != numNodes) {
        std::stringstream ss;
        ss << "NodeFile::setCoordinates: number of given nodes must be "
            << numNodes;
        const std::string errorMsg(ss.str());
        Finley_setError(VALUE_ERROR, errorMsg.c_str());
    } else {
        const size_t numDim_size=numDim*sizeof(double);
        Finley_increaseStatus(this);
        escript::Data& newX = *const_cast<escript::Data*>(&cNewX);
#pragma omp parallel for
        for (int n=0; n<numNodes; n++) {
            memcpy(&(Coordinates[INDEX2(0,n,numDim)]), newX.getSampleDataRO(n), numDim_size);
        }
    }
}

/// sets tags to newTag where mask>0
void NodeFile::setTags(const int newTag, const escript::Data& cMask)
{
    Finley_resetError();

    if (1 != cMask.getDataPointSize()) {
       Finley_setError(TYPE_ERROR, "NodeFile::setTags: number of components of mask must be 1.");
       return;
    } else if (cMask.getNumDataPointsPerSample() != 1 ||
            cMask.getNumSamples() != numNodes) {
       Finley_setError(TYPE_ERROR, "NodeFile::setTags: illegal number of samples of mask Data object");
       return;
    }

    escript::Data& mask = *const_cast<escript::Data*>(&cMask);
#pragma omp parallel for
    for (int n=0; n<numNodes; n++) {
         if (mask.getSampleDataRO(n)[0] > 0)
             Tag[n]=newTag;
    }
    updateTagList();
}

std::pair<int,int> NodeFile::getDOFRange() const
{
    std::pair<int,int> result(
            Finley_Util_getMinMaxInt(1, numNodes, globalDegreesOfFreedom));
    if (result.second < result.first) {
        result.first = -1;
        result.second = 0;
    }
    return result;
}

std::pair<int,int> NodeFile::getGlobalIdRange() const
{
    std::pair<int,int> result(Finley_Util_getMinMaxInt(1, numNodes, Id));

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
    std::pair<int,int> result(
            Finley_Util_getMinMaxInt(1, numNodes, globalDegreesOfFreedom));

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
    std::pair<int,int> result(
            Finley_Util_getMinMaxInt(1, numNodes, globalNodesIndex));

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
        Finley_setError(TYPE_ERROR, "NodeFile::copyTable: dimensions of node files don't match");
        return;
    }
    if (numNodes < in->numNodes+offset) {
        Finley_setError(MEMORY_ERROR, "NodeFile::copyTable: node table is too small.");
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

void NodeFile::gather_global(int* index, const NodeFile* in)
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
            in->MPIInfo->msg_tag_counter+=4;
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
        gatherEntries(numNodes, index, distribution[buffer_rank],
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
             in->MPIInfo->msg_tag_counter+=4;
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
            Finley_setError(VALUE_ERROR, errorMsg.c_str());
        }
    }
    delete[] Id_buffer;
    delete[] Tag_buffer;
    delete[] globalDegreesOfFreedom_buffer;
    delete[] Coordinates_buffer;
    // make sure that the error is global
    Esys_MPIInfo_noError(in->MPIInfo);
}

void NodeFile::assignMPIRankToDOFs(Esys_MPI_rank* mpiRankOfDOF,
                                   int *distribution)
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

int NodeFile::prepareLabeling(int* mask, std::vector<int>& buffer,
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
            if (!mask || mask[n]>-1) {
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
    loc_offsets[MPIInfo->rank]=prepareLabeling(NULL, DOF_buffer, distribution,
                                               false);
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

    std::vector<bool_t> set_new_DOF(numNodes, TRUE);

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
                set_new_DOF[n]=FALSE;
            }
        }
        if (p<MPIInfo->size-1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&DOF_buffer[0], DOF_buffer.size(), MPI_INT,
                    dest, MPIInfo->msg_tag_counter, source,
                    MPIInfo->msg_tag_counter, MPIInfo->comm, &status);
#endif
            MPIInfo->msg_tag_counter+=1;
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }

    return new_numGlobalDOFs;
}

int NodeFile::createDenseNodeLabeling(int* node_distribution,
                                      const int* dof_distribution) 
{
    const int UNSET_ID=-1, SET_ID=1;
    const int myFirstDOF=dof_distribution[MPIInfo->rank];
    const int myLastDOF=dof_distribution[MPIInfo->rank+1];

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
    MPI_Allgather(&myNewNumNodes, 1, MPI_INT, node_distribution, 1, MPI_INT,
                  MPIInfo->comm);
#else
    node_distribution[0]=myNewNumNodes;
#endif

    int globalNumNodes=0;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int itmp=node_distribution[p];
        node_distribution[p]=globalNumNodes;
        globalNumNodes+=itmp;
    }
    node_distribution[MPIInfo->size]=globalNumNodes;

    // offset node buffer
#pragma omp parallel for
    for (int n=0; n<my_buffer_len; n++)
        Node_buffer[n+header_len]+=node_distribution[MPIInfo->rank];

    // now we send this buffer around to assign global node index
#ifdef ESYS_MPI
    int dest=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank + 1);
    int source=Esys_MPIInfo_mod(MPIInfo->size, MPIInfo->rank - 1);
#endif
    int buffer_rank=MPIInfo->rank;
    for (int p=0; p<MPIInfo->size; ++p) {
        const int nodeID_0=Node_buffer[0];
        const int nodeID_1=Node_buffer[1];
        const int dof0=dof_distribution[buffer_rank];
        const int dof1=dof_distribution[buffer_rank+1];
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
            MPIInfo->msg_tag_counter+=1;
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }
    return globalNumNodes;
}

int NodeFile::createDenseReducedLabeling(int* reducedMask, bool useNodes) 
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
            MPIInfo->msg_tag_counter+=1;
        }
        buffer_rank=Esys_MPIInfo_mod(MPIInfo->size, buffer_rank-1);
    }
    return new_numGlobalReduced;
}

} // namespace finley

