
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/NodeFile.h>
extern "C" {
#include <paso/Coupler.h>
#include <paso/Distribution.h>
}
#include <escript/Data.h>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif


using namespace std;

namespace ripley {

NodeFile::NodeFile(dim_t numDim, Esys_MPIInfo* mpiInfo) :
    m_numDim(numDim),
    m_status(0),
    m_nodesDistribution(NULL),
    m_reducedNodesDistribution(NULL),
    m_degreesOfFreedomDistribution(NULL),
    m_reducedDegreesOfFreedomDistribution(NULL),
    m_degreesOfFreedomConnector(NULL),
    m_reducedDegreesOfFreedomConnector(NULL)
{
    m_mpiInfo = Esys_MPIInfo_getReference(mpiInfo);
}

NodeFile::~NodeFile()
{
    Esys_MPIInfo_free(m_mpiInfo);
}

#ifdef USE_NETCDF
void NodeFile::readFromNetCDF(NcFile &dataFile, dim_t numNodes)
{
    // clear and resize data vectors
    m_id.assign(numNodes, -1);
    m_tag.assign(numNodes, -1);
    m_globalDegreesOfFreedom.assign(numNodes, -1);
    m_globalNodesIndex.assign(numNodes, -1);
    m_globalReducedDOFIndex.assign(numNodes, -1);
    m_globalReducedNodesIndex.assign(numNodes, -1);
    m_coordinates.assign(numNodes*m_numDim, 0.);

    NcVar *ncVar;

    // Id
    if (! ( ncVar = dataFile.get_var("Nodes_Id")) )
        throw RipleyException("get_var(Nodes_Id)");
    if (!ncVar->get(&m_id[0], numNodes))
        throw RipleyException("get(Nodes_Id)");
    // Tag
    if (! ( ncVar = dataFile.get_var("Nodes_Tag")) )
        throw RipleyException("get_var(Nodes_Tag)");
    if (!ncVar->get(&m_tag[0], numNodes))
        throw RipleyException("get(Nodes_Tag)");
    // gDOF
    if (! ( ncVar = dataFile.get_var("Nodes_gDOF")) )
        throw RipleyException("get_var(Nodes_gDOF)");
    if (!ncVar->get(&m_globalDegreesOfFreedom[0], numNodes))
        throw RipleyException("get(Nodes_gDOF)");
    // gNI
    if (! ( ncVar = dataFile.get_var("Nodes_gNI")) )
        throw RipleyException("get_var(Nodes_gNI)");
    if (!ncVar->get(&m_globalNodesIndex[0], numNodes))
        throw RipleyException("get(Nodes_gNI)");
    // grDfI
    if (! ( ncVar = dataFile.get_var("Nodes_grDfI")) )
        throw RipleyException("get_var(Nodes_grDfI)");
    if (!ncVar->get(&m_globalReducedDOFIndex[0], numNodes))
        throw RipleyException("get(Nodes_grDfI)");
    // grNI
    if (! ( ncVar = dataFile.get_var("Nodes_grNI")) )
        throw RipleyException("get_var(Nodes_grNI)");
    if (!ncVar->get(&m_globalReducedNodesIndex[0], numNodes))
        throw RipleyException("get(Nodes_grNI)");
    // Coordinates
    if (!(ncVar = dataFile.get_var("Nodes_Coordinates")))
        throw RipleyException("get_var(Nodes_Coordinates)");
    if (!ncVar->get(&m_coordinates[0], numNodes, m_numDim))
        throw RipleyException("get(Nodes_Coordinates)");
}

void NodeFile::dumpToNetCDF(NcFile &dataFile)
{
    vector<NcDim*> ncdims(3);
    NcVar *ncVar;
    const string msgPrefix("dump: netCDF operation failed - ");
    const dim_t numNodes = getNumNodes();

    if (!dataFile.add_att("numNodes", numNodes))
        throw RipleyException(msgPrefix+"add_att(numNodes)");
    if (! (ncdims[0] = dataFile.add_dim("numNodes", numNodes)) )
        throw RipleyException(msgPrefix+"add_dim(numNodes)");
    if (! (ncdims[1] = dataFile.add_dim("numDim", m_numDim)) )
        throw RipleyException(msgPrefix+"add_dim(numDim)");
    if (! (ncdims[2] = dataFile.add_dim("mpi_size_plus_1", m_mpiInfo->size+1)) )
        throw RipleyException(msgPrefix+"add_dim(mpi_size)");

    // nodeDistribution
    if (!(ncVar = dataFile.add_var("Nodes_NodeDistribution", ncInt, ncdims[2])))
        throw RipleyException(msgPrefix+"add_var(Nodes_NodeDistribution)");
    if (!ncVar->put(&m_nodesDistribution->first_component[0], m_mpiInfo->size+1))
        throw RipleyException(msgPrefix+"put(Nodes_NodeDistribution)");

    // degreesOfFreedomDistribution
    if (! ( ncVar = dataFile.add_var("Nodes_DofDistribution", ncInt, ncdims[2])) )
        throw RipleyException(msgPrefix+"add_var(Nodes_DofDistribution)");
    if (!ncVar->put(&m_degreesOfFreedomDistribution->first_component[0], m_mpiInfo->size+1))
        throw RipleyException(msgPrefix+"put(Nodes_DofDistribution)");

    // Only write data vectors if non-empty because netCDF doesn't like empty
    // arrays (it treats them as NC_UNLIMITED)
    if (numNodes > 0) {
        // Id
        if (! (ncVar = dataFile.add_var("Nodes_Id", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_Id)");
        if (!ncVar->put(&m_id[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_Id)");

        // Tag
        if (! (ncVar = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_Tag)");
        if (!ncVar->put(&m_tag[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_Tag)");

        // gDOF
        if (! (ncVar = dataFile.add_var("Nodes_gDOF", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_gDOF)");
        if (!ncVar->put(&m_globalDegreesOfFreedom[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_gDOF)");

        // global node index
        if (! (ncVar = dataFile.add_var("Nodes_gNI", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_gNI)");
        if (!ncVar->put(&m_globalNodesIndex[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_gNI)");

        // grDof
        if (! (ncVar = dataFile.add_var("Nodes_grDfI", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_grDfI)");
        if (!ncVar->put(&m_globalReducedDOFIndex[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_grDfI)");

        // grNI
        if (! (ncVar = dataFile.add_var("Nodes_grNI", ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_grNI)");
        if (!ncVar->put(&m_globalReducedNodesIndex[0], numNodes))
            throw RipleyException(msgPrefix+"put(Nodes_grNI)");

        // Coordinates
        if (! (ncVar = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1])) )
            throw RipleyException(msgPrefix+"add_var(Nodes_Coordinates)");
        if (!ncVar->put(&(m_coordinates[INDEX2(0,0,m_numDim)]), numNodes, m_numDim))
            throw RipleyException(msgPrefix+"put(Nodes_Coordinates)");
    }
}
#endif

dim_t NodeFile::getGlobalNumNodes() const
{
    return Paso_Distribution_getGlobalNumComponents(m_nodesDistribution);
}

dim_t NodeFile::getGlobalNumReducedNodes() const
{
    return Paso_Distribution_getGlobalNumComponents(m_reducedNodesDistribution);
}

dim_t NodeFile::getNumDegreesOfFreedom() const
{
    return Paso_Distribution_getMyNumComponents(m_degreesOfFreedomDistribution);
}

dim_t NodeFile::getNumReducedDegreesOfFreedom() const
{
    return Paso_Distribution_getMyNumComponents(m_reducedDegreesOfFreedomDistribution);
}

index_t NodeFile::getFirstNode() const
{
    return Paso_Distribution_getFirstComponent(m_nodesDistribution);
}

index_t NodeFile::getLastNode() const
{
    return Paso_Distribution_getLastComponent(m_nodesDistribution);
}

index_t NodeFile::getFirstReducedNode() const
{
    return Paso_Distribution_getFirstComponent(m_reducedNodesDistribution);
}

index_t NodeFile::getLastReducedNode() const
{
    return Paso_Distribution_getLastComponent(m_reducedNodesDistribution);
}

void NodeFile::assembleCoordinates(escript::Data &arg) const
{
    dim_t numNodes = getNumNodes();
    dim_t numDim = getNumDim();
    escriptDataC x = arg.getDataC();
    if (!arg.isExpanded())
        throw RipleyException("setToX: Expanded Data object expected");
    if (!isDataPointShapeEqual(&x, 1, &numDim))
        throw RipleyException("setToX: Invalid Data object shape");
    if (!numSamplesEqual(&x, 1, numNodes))
        throw RipleyException("setToX: Illegal number of samples in Data object");

    size_t dimSize = getNumDim() * sizeof(double);
    dim_t n;
    requireWrite(&x);
#pragma omp parallel for private(n)
        for (n = 0; n < numNodes; n++)
            memcpy(getSampleDataRWFast(&x, n), &(m_coordinates[INDEX2(0, n, getNumDim())]), dimSize);
}

void NodeFile::setCoordinates(const escript::Data &coords)
{
    dim_t numNodes = getNumNodes();
    escriptDataC newX = coords.getDataC();
    if (getDataPointSize(&newX) != getNumDim())
        throw RipleyException("setCoordinates: Wrong dimensionality of new coordinates");
    if (!numSamplesEqual(&newX, 1, numNodes))
        throw RipleyException("setCoordinates: Illegal number of samples in Data object");

    size_t dimSize = getNumDim() * sizeof(double);
    m_status++;
    dim_t n;
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < numNodes; n++)
        memcpy(&m_coordinates[INDEX2(0, n, getNumDim())], getSampleDataROFast(&newX, n), dimSize);
}

void NodeFile::updateTagsInUse()
{
    m_tagsInUse = getUniqueValues(m_tag, m_mpiInfo);
}

void NodeFile::setTags(int newTag, const escript::Data &mask)
{
    dim_t numNodes = getNumNodes();
    register dim_t n;
    escriptDataC maskC = mask.getDataC();

    if (1 != getDataPointSize(&maskC))
        throw RipleyException("setTags: Number of components in mask must be 1");
    if (!numSamplesEqual(&maskC, 1, numNodes))
        throw RipleyException("setTags: Illegal number of samples in mask Data object");

#pragma omp parallel for schedule(static) private(n)
    for (n = 0; n < numNodes; n++) {
        const double *maskArray = getSampleDataRO(&maskC, n);
        if (maskArray[0] > 0)
            m_tag[n] = newTag;
    }
    updateTagsInUse();
}

RankVector NodeFile::getOwnerOfDOFs(const IndexVector &distribution)
{
    // first we obtain the min and max DOF on this processor to reduce
    // costs for searching
    IndexPair dofRange = getDOFRange();

    Esys_MPI_rank p, pMin = m_mpiInfo->size, pMax = -1;
    for (p = 0; p < m_mpiInfo->size; ++p) {
        if (distribution[p] <= dofRange.first)
            pMin = p;
        if (distribution[p] <= dofRange.second)
            pMax = p;
    }
    dim_t n;
    RankVector mpiRankOfDOF(getNumNodes(), 0);

#pragma omp parallel for private(n,p) schedule(static)
    for (n = 0; n < getNumNodes(); ++n) {
        dim_t k = m_globalDegreesOfFreedom[n];
        for (p = pMin; p <= pMax; ++p) {
            if (k < distribution[p + 1]) {
                mpiRankOfDOF[n] = p;
                break;
            }
        }
    }
    return mpiRankOfDOF;
}

IndexPair NodeFile::getDOFRange() const
{
    if (m_globalDegreesOfFreedom.size() == 0)
        return IndexPair(-1, 0);

    return getMinMax(m_globalDegreesOfFreedom);
}

IndexPair NodeFile::getGlobalNodeIndexRange() const
{
    IndexPair minMax = getGlobalMinMax(m_globalNodesIndex, m_mpiInfo);
    if (minMax.first > minMax.second)
        return IndexPair(-1, 0);
    return minMax;
}

IndexPair NodeFile::getGlobalIdRange() const
{
    IndexPair minMax = getGlobalMinMax(m_id, m_mpiInfo);
    if (minMax.first > minMax.second)
        return IndexPair(-1, 0);
    return minMax;
}

IndexPair NodeFile::getGlobalDOFRange() const
{
    IndexPair minMax = getGlobalMinMax(m_globalDegreesOfFreedom, m_mpiInfo);
    if (minMax.first > minMax.second)
        return IndexPair(-1, 0);
    return minMax;
}

dim_t NodeFile::createDenseDOFLabeling()
{
    const index_t UNSET = -1;
    const index_t SET = 1;

    // get the global range of node IDs
    IndexPair minMaxDof = getGlobalDOFRange();

    RankVector distribution(m_mpiInfo->size+1);
    vector<dim_t> offsets(m_mpiInfo->size);
    vector<dim_t> loc_offsets(m_mpiInfo->size, 0);

    // distribute the range of node IDs
    dim_t bufferLen = Esys_MPIInfo_setDistribution(m_mpiInfo,
            minMaxDof.first, minMaxDof.second, &distribution[0]);
    dim_t myDOFs = distribution[m_mpiInfo->rank+1]-distribution[m_mpiInfo->rank];
    // fill DOF_buffer by the UNSET marker to check if nodes are defined
    IndexVector DOF_buffer(bufferLen, UNSET);

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank + 1);
    Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank - 1);
#endif
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        if (p > 0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&DOF_buffer[0], bufferLen, MPI_INT, dest,
                    m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
#endif
            m_mpiInfo->msg_tag_counter++;
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
        const index_t dof_0 = distribution[curRank];
        const index_t dof_1 = distribution[curRank + 1];
#pragma omp parallel for schedule(static)
        for (dim_t n = 0; n < getNumNodes(); n++) {
            const index_t k = m_globalDegreesOfFreedom[n];
            if ((dof_0 <= k) && (k < dof_1))
                DOF_buffer[k-dof_0] = SET;
        }
    }
    // count the entries in the DOF_buffer
    dim_t myNewDOFs = 0;
    for (dim_t n = 0; n < myDOFs; ++n) {
        if (DOF_buffer[n] == SET) {
            DOF_buffer[n] = myNewDOFs;
            myNewDOFs++;
        }
    }
    loc_offsets[m_mpiInfo->rank] = myNewDOFs;
    dim_t new_numGlobalDOFs = 0;

    if (m_mpiInfo->size > 1) {
#ifdef ESYS_MPI
        MPI_Allreduce(&loc_offsets[0], &offsets[0], m_mpiInfo->size, MPI_INT, MPI_SUM, m_mpiInfo->comm);
        for (dim_t n = 0; n < m_mpiInfo->size; ++n) {
            loc_offsets[n] = new_numGlobalDOFs;
            new_numGlobalDOFs += offsets[n];
        }
#endif
    } else {
        new_numGlobalDOFs = loc_offsets[0];
        loc_offsets[0] = 0;
    }

#pragma omp parallel for schedule(static)
    for (dim_t n = 0; n < myDOFs; ++n)
        DOF_buffer[n] += loc_offsets[m_mpiInfo->rank];

    // now entries are collected from the buffer again by sending them around
    // in a circle
    resetGlobalDegreesOfFreedom(DOF_buffer, distribution);

    return new_numGlobalDOFs;
}

void NodeFile::resetGlobalDegreesOfFreedom(IndexVector &newGlobalDOFid,
                                           const IndexVector &distribution)
{
    const dim_t numNodes = getNumNodes();
    vector<bool> setNewDOF(numNodes, true);
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        const index_t firstDOF = distribution[curRank];
        const index_t lastDOF = distribution[curRank+1];
#pragma omp parallel for schedule(static)
        for (dim_t n = 0; n < numNodes; n++) {
            const index_t k = m_globalDegreesOfFreedom[n];
            if (setNewDOF[n] && (firstDOF <= k) && (k < lastDOF)) {
                m_globalDegreesOfFreedom[n] = newGlobalDOFid[k-firstDOF];
                setNewDOF[n] = false;
            }
        }
        if (p < m_mpiInfo->size - 1) { // the last send can be skipped
#ifdef ESYS_MPI
            Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank+1);
            Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank-1);
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFid[0], newGlobalDOFid.size(),
                    MPI_INT, dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
            m_mpiInfo->msg_tag_counter++;
#endif
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
    }
}

dim_t NodeFile::createDenseNodeLabeling(IndexVector &nodeDistribution,
                                        const IndexVector &dofDistribution)
{
    const index_t UNSET = -1;
    const index_t SET = 1;
    const dim_t HEADER_LEN = 2;

    // find the range of node IDs controlled by me
    index_t myFirstDOF = dofDistribution[m_mpiInfo->rank];
    index_t myLastDOF = dofDistribution[m_mpiInfo->rank + 1];
    index_t maxId = -INDEX_T_MAX;
    index_t minId = INDEX_T_MAX;
#pragma omp parallel
    {
        index_t loc_maxId = maxId;
        index_t loc_minId = minId;
#pragma omp for schedule(static)
        for (dim_t n = 0; n < getNumNodes(); n++) {
            const index_t dof = m_globalDegreesOfFreedom[n];
            if ((myFirstDOF <= dof) && (dof < myLastDOF)) {
                loc_maxId = max(loc_maxId, m_id[n]);
                loc_minId = min(loc_minId, m_id[n]);
            }
        }
#pragma omp critical
        {
            maxId = max(loc_maxId, maxId);
            minId = min(loc_minId, minId);
        }
    }

    dim_t myBufferLen = (maxId >= minId ? maxId-minId+1 : 0);
    dim_t bufferLen;

#ifdef ESYS_MPI
    MPI_Allreduce(&myBufferLen, &bufferLen, 1, MPI_INT, MPI_MAX, m_mpiInfo->comm);
#else
    bufferLen = myBufferLen;
#endif

    IndexVector nodeBuffer(bufferLen+HEADER_LEN, UNSET);
    // mark and count the nodes in use
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (dim_t n = 0; n < getNumNodes(); n++) {
            const index_t dof = m_globalDegreesOfFreedom[n];
            if ((myFirstDOF <= dof) && (dof < myLastDOF))
                nodeBuffer[m_id[n] - minId + HEADER_LEN] = SET;
            m_globalNodesIndex[n] = -1;
        }
    }
    dim_t myNewNumNodes = 0;
    for (dim_t n = 0; n < myBufferLen; n++) {
        if (nodeBuffer[HEADER_LEN + n] == SET) {
            nodeBuffer[HEADER_LEN + n] = myNewNumNodes;
            myNewNumNodes++;
        }
    }
    // make the local number of nodes globally available
#ifdef ESYS_MPI
    MPI_Allgather(&myNewNumNodes, 1, MPI_INT, &nodeDistribution[0], 1, MPI_INT, m_mpiInfo->comm);
#else
    nodeDistribution[0] = myNewNumNodes;
#endif

    dim_t globalNumNodes = 0;
    for (Esys_MPI_rank p = 0; p < m_mpiInfo->size; ++p) {
        index_t itmp = nodeDistribution[p];
        nodeDistribution[p] = globalNumNodes;
        globalNumNodes += itmp;
    }
    nodeDistribution[m_mpiInfo->size] = globalNumNodes;

    // offset node buffer
#pragma omp for schedule(static)
    for (dim_t n = 0; n < myBufferLen; n++)
        nodeBuffer[n + HEADER_LEN] += nodeDistribution[m_mpiInfo->rank];

    // now we send this buffer around to assign a global node index
    nodeBuffer[0] = minId;
    nodeBuffer[1] = maxId;
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (Esys_MPI_rank p = 0; p < m_mpiInfo->size; ++p) {
        const index_t nodeID_0 = nodeBuffer[0];
        const index_t nodeID_1 = nodeBuffer[1];
        const index_t dof_0 = dofDistribution[curRank];
        const index_t dof_1 = dofDistribution[curRank + 1];
        if (nodeID_0 <= nodeID_1) {
#pragma omp for schedule(static)
            for (dim_t n = 0; n < getNumNodes(); n++) {
                const index_t dof = m_globalDegreesOfFreedom[n];
                const index_t id_n = m_id[n] - nodeID_0;
                if ((dof_0 <= dof) && (dof < dof_1) && (id_n >= 0) && (id_n <= nodeID_1 - nodeID_0))
                    m_globalNodesIndex[n] = nodeBuffer[id_n + HEADER_LEN];
            }
        }
        if (p < m_mpiInfo->size - 1) { // the last send can be skipped
#ifdef ESYS_MPI
            Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank + 1);
            Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank - 1);
            MPI_Status status;
            MPI_Sendrecv_replace(&nodeBuffer[0], bufferLen + HEADER_LEN,
                    MPI_INT, dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
            m_mpiInfo->msg_tag_counter++;
#endif
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank - 1);
    }
    return globalNumNodes;
}

dim_t NodeFile::createDenseReducedDOFLabeling(const IndexVector &mask)
{
    const index_t UNSET = -1;
    const index_t SET = 1;
    dim_t n, globalNumReducedDOFs;

    // get the global range of node IDs
    IndexPair minMaxDof = getGlobalDOFRange();

    RankVector distribution(m_mpiInfo->size+1);
    vector<dim_t> offsets(m_mpiInfo->size);
    vector<dim_t> loc_offsets(m_mpiInfo->size, 0);

    // distribute the range of node IDs
    dim_t bufferLen = Esys_MPIInfo_setDistribution(m_mpiInfo, minMaxDof.first, minMaxDof.second, &distribution[0]);
    dim_t myDOFs = distribution[m_mpiInfo->rank+1]-distribution[m_mpiInfo->rank];
    // fill DOF_buffer by the UNSET marker to check if nodes are defined    
    IndexVector DOF_buffer(bufferLen, UNSET);

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank+1);
    Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank-1);
#endif
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        if (p > 0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&DOF_buffer[0], bufferLen, MPI_INT, dest,
                    m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
            m_mpiInfo->msg_tag_counter++;
#endif
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
        index_t dof_0 = distribution[curRank];
        index_t dof_1 = distribution[curRank + 1];
#pragma omp parallel for private(n) schedule(static)
        for (n = 0; n < getNumNodes(); n++) {
            if (mask[n] > -1) {
                index_t k = m_globalDegreesOfFreedom[n];
                if ((dof_0 <= k) && (k < dof_1))
                    DOF_buffer[k - dof_0] = SET;
            }
        }
    }
    // count the entries in the DOF_buffer
    // TODO: OMP parallel
    dim_t myNewDOFs = 0;
    for (n = 0; n < myDOFs; ++n) {
        if (DOF_buffer[n] == SET) {
            DOF_buffer[n] = myNewDOFs;
            myNewDOFs++;
        }
    }
    loc_offsets[m_mpiInfo->rank] = myNewDOFs;
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], m_mpiInfo->size, MPI_INT, MPI_SUM, m_mpiInfo->comm);
    globalNumReducedDOFs = 0;
    for (n = 0; n < m_mpiInfo->size; ++n) {
        loc_offsets[n] = globalNumReducedDOFs;
        globalNumReducedDOFs += offsets[n];
    }
#else
    globalNumReducedDOFs = loc_offsets[0];
    loc_offsets[0] = 0;
#endif
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < myDOFs; ++n)
        DOF_buffer[n] += loc_offsets[m_mpiInfo->rank];
    // now entries are collected from the buffer again by sending the
    // entries around in a circle
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < getNumNodes(); ++n)
        m_globalReducedDOFIndex[n] = loc_offsets[0] - 1;
#ifdef ESYS_MPI
    dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank + 1);
    source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank - 1);
#endif
    curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        index_t dof_0 = distribution[curRank];
        index_t dof_1 = distribution[curRank + 1];
#pragma omp parallel for private(n) schedule(static)
        for (n = 0; n < getNumNodes(); n++) {
            if (mask[n] > -1) {
                index_t k = m_globalDegreesOfFreedom[n];
                if ((dof_0 <= k) && (k < dof_1))
                    m_globalReducedDOFIndex[n] = DOF_buffer[k-dof_0];
            }
        }
        if (p < m_mpiInfo->size - 1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&DOF_buffer[0], bufferLen, MPI_INT, dest,
                    m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
#endif
            m_mpiInfo->msg_tag_counter += 1;
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
    }

    return globalNumReducedDOFs;
}

dim_t NodeFile::createDenseReducedNodeLabeling(const IndexVector &mask)
{
    const index_t UNSET = -1;
    const index_t SET = 1;
    dim_t n, globalNumReducedNodes;

    // get the global range of node IDs
    IndexPair minMaxNode = getGlobalNodeIndexRange();

    RankVector distribution(m_mpiInfo->size+1);
    vector<dim_t> offsets(m_mpiInfo->size);
    vector<dim_t> loc_offsets(m_mpiInfo->size, 0);

    // distribute the range of node IDs
    dim_t bufferLen = Esys_MPIInfo_setDistribution(m_mpiInfo, minMaxNode.first, minMaxNode.second, &distribution[0]);
    dim_t myNodes = distribution[m_mpiInfo->rank+1]-distribution[m_mpiInfo->rank];
    // fill Nodes_buffer by the UNSET marker to check if nodes are defined
    IndexVector Nodes_buffer(bufferLen, UNSET);

    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank+1);
    Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank-1);
#endif
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        if (p > 0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&Nodes_buffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
#endif
            m_mpiInfo->msg_tag_counter++;
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
        index_t node_0 = distribution[curRank];
        index_t node_1 = distribution[curRank + 1];
#pragma omp parallel for private(n) schedule(static)
        for (n = 0; n < getNumNodes(); n++) {
            if (mask[n] > -1) {
                index_t k = m_globalNodesIndex[n];
                if ((node_0 <= k) && (k < node_1))
                    Nodes_buffer[k - node_0] = SET;
            }
        }
    }
    // count the entries in the Nodes_buffer
    // TODO: OMP parallel
    dim_t myNewNodes = 0;
    for (n = 0; n < myNodes; ++n) {
        if (Nodes_buffer[n] == SET) {
            Nodes_buffer[n] = myNewNodes;
            myNewNodes++;
        }
    }
    loc_offsets[m_mpiInfo->rank] = myNewNodes;
#ifdef ESYS_MPI
    MPI_Allreduce(&loc_offsets[0], &offsets[0], m_mpiInfo->size, MPI_INT, MPI_SUM, m_mpiInfo->comm);
    globalNumReducedNodes = 0;
    for (n = 0; n < m_mpiInfo->size; ++n) {
        loc_offsets[n] = globalNumReducedNodes;
        globalNumReducedNodes += offsets[n];
    }
#else
    globalNumReducedNodes = loc_offsets[0];
    loc_offsets[0] = 0;
#endif
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < myNodes; ++n)
        Nodes_buffer[n] += loc_offsets[m_mpiInfo->rank];
    // now entries are collected from the buffer again by sending the
    // entries around in a circle
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < getNumNodes(); ++n)
        m_globalReducedNodesIndex[n] = loc_offsets[0] - 1;
#ifdef ESYS_MPI
    dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank + 1);
    source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank - 1);
#endif
    curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        index_t node_0 = distribution[curRank];
        index_t node_1 = distribution[curRank + 1];
#pragma omp parallel for private(n) schedule(static)
        for (n = 0; n < getNumNodes(); n++) {
            if (mask[n] > -1) {
                index_t k = m_globalNodesIndex[n];
                if ((node_0 <= k) && (k < node_1))
                    m_globalReducedNodesIndex[n] = Nodes_buffer[k - node_0];
            }
        }
        if (p < m_mpiInfo->size - 1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&Nodes_buffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm, &status);
#endif
            m_mpiInfo->msg_tag_counter += 1;
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
    }

    return globalNumReducedNodes;
}

NodeFile_ptr NodeFile::gatherGlobal(IndexVector &index)
{
    dim_t n;

    // get the global range of node IDs
    IndexPair minMaxId = getGlobalIdRange();
    const index_t UNDEFINED_NODE = minMaxId.first - 1;

    RankVector distribution(m_mpiInfo->size + 1);

    // distribute the range of node IDs
    dim_t bufferLen = Esys_MPIInfo_setDistribution(m_mpiInfo, minMaxId.first, minMaxId.second, &distribution[0]);

    // fill idBuffer by the UNDEFINED_NODE marker to check if nodes are defined
    IndexVector idBuffer(bufferLen, UNDEFINED_NODE);
    IndexVector tagBuffer(bufferLen);
    IndexVector gDOFbuffer(bufferLen);
    vector<double> coordBuffer(bufferLen*getNumDim());
    // fill the buffer by sending portions around in a circle
#ifdef ESYS_MPI
    MPI_Status status;
    Esys_MPI_rank dest = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank+1);
    Esys_MPI_rank source = Esys_MPIInfo_mod(m_mpiInfo->size, m_mpiInfo->rank-1);
#endif
    Esys_MPI_rank curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        if (p > 0) { // the initial send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&idBuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm,
                    &status);
            MPI_Sendrecv_replace(&tagBuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter + 1, source,
                    m_mpiInfo->msg_tag_counter + 1,
                    m_mpiInfo->comm, &status);
            MPI_Sendrecv_replace(&gDOFbuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter + 2, source,
                    m_mpiInfo->msg_tag_counter + 2,
                    m_mpiInfo->comm, &status);
            MPI_Sendrecv_replace(&coordBuffer[0],
                    bufferLen*getNumDim(), MPI_DOUBLE, dest,
                    m_mpiInfo->msg_tag_counter + 3, source,
                    m_mpiInfo->msg_tag_counter + 3,
                    m_mpiInfo->comm, &status);
            m_mpiInfo->msg_tag_counter += 4;
#endif
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank - 1);
        scatterEntries(m_id, distribution[curRank],
                distribution[curRank+1], idBuffer, tagBuffer, gDOFbuffer,
                coordBuffer);
    }

    NodeFile_ptr out(new NodeFile(getNumDim(), m_mpiInfo));
    out->m_id.assign(index.size(), -1);
    out->m_tag.assign(index.size(), -1);
    out->m_globalDegreesOfFreedom.assign(index.size(), -1);
    out->m_coordinates.assign(index.size()*getNumDim(), 0.);
    out->m_globalReducedDOFIndex.assign(index.size(), -1);
    out->m_globalNodesIndex.assign(index.size(), -1);
    out->m_globalReducedNodesIndex.assign(index.size(), -1);

    // now entries are collected from the buffer again by sending them around
    // in a circle
    curRank = m_mpiInfo->rank;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p) {
        out->gatherEntries(index, distribution[curRank],
                distribution[curRank+1], idBuffer, tagBuffer, gDOFbuffer,
                coordBuffer);
        if (p < m_mpiInfo->size - 1) { // the last send can be skipped
#ifdef ESYS_MPI
            MPI_Sendrecv_replace(&idBuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter, source,
                    m_mpiInfo->msg_tag_counter, m_mpiInfo->comm,
                    &status);
            MPI_Sendrecv_replace(&tagBuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter + 1, source,
                    m_mpiInfo->msg_tag_counter + 1,
                    m_mpiInfo->comm, &status);
            MPI_Sendrecv_replace(&gDOFbuffer[0], bufferLen, MPI_INT,
                    dest, m_mpiInfo->msg_tag_counter + 2, source,
                    m_mpiInfo->msg_tag_counter + 2,
                    m_mpiInfo->comm, &status);
            MPI_Sendrecv_replace(&coordBuffer[0],
                    bufferLen*getNumDim(), MPI_DOUBLE, dest,
                    m_mpiInfo->msg_tag_counter + 3, source,
                    m_mpiInfo->msg_tag_counter + 3,
                    m_mpiInfo->comm, &status);
            m_mpiInfo->msg_tag_counter += 4;
#endif
        }
        curRank = Esys_MPIInfo_mod(m_mpiInfo->size, curRank-1);
    }
    // check if all nodes are set
#pragma omp parallel for private(n) schedule(static)
    for (n = 0; n < out->getNumNodes(); ++n) {
        if (out->m_id[n] == UNDEFINED_NODE) {
            stringstream msg;
            msg << "gatherGlobal: Node id " << out->m_id[n] <<
                " at position " << n << " is referenced but not defined";
            throw RipleyException(msg.str());
        }
    }
    return out;
}

void NodeFile::gatherEntries(const IndexVector &index, index_t minIndex,
        index_t maxIndex, IndexVector &idIn, IndexVector &tagIn,
        IndexVector &gDOFin, vector<double> &coordIn)
{
    const dim_t numNodes = index.size();
    const index_t range = maxIndex-minIndex;
    const size_t numDimSize = m_numDim*sizeof(double);

#pragma omp parallel for schedule(static)
    for (dim_t i = 0; i < numNodes; i++) {
        const index_t k = index[i] - minIndex;
        if ((k >= 0) && (k < range)) {
            m_id[i] = idIn[k];
            m_tag[i] = tagIn[k];
            m_globalDegreesOfFreedom[i] = gDOFin[k];
            memcpy(&(m_coordinates[INDEX2(0, i, m_numDim)]),
                    &(coordIn[INDEX2(0, k, m_numDim)]), numDimSize);
        }
    }
}

void NodeFile::scatterEntries(const IndexVector &index, index_t minIndex,
        index_t maxIndex, IndexVector &idOut, IndexVector &tagOut,
        IndexVector &gDOFout, vector<double> &coordOut)
{
    dim_t i;
    const dim_t numNodes = getNumNodes();
    const index_t range = maxIndex-minIndex;
    const size_t numDimSize = m_numDim*sizeof(double);

#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < numNodes; i++) {
        const index_t k = index[i] - minIndex;
        if (k >= 0 && k < range) {
            idOut[k] = m_id[i];
            tagOut[k] = m_tag[i];
            gDOFout[k] = m_globalDegreesOfFreedom[i];
            memcpy(&(coordOut[INDEX2(0, k, m_numDim)]),
                    &(m_coordinates[INDEX2(0, i, m_numDim)]), numDimSize);
        }
    }
}

void NodeFile::swapEntries(IndexVector &idIn, IndexVector &tagIn,
        IndexVector &gDOFin, vector<double> &coordIn)
{
    m_id.swap(idIn);
    m_tag.swap(tagIn);
    m_globalDegreesOfFreedom.swap(gDOFin);
    m_coordinates.swap(coordIn);
    m_globalReducedDOFIndex.assign(idIn.size(), -1);
    m_globalNodesIndex.assign(idIn.size(), -1);
    m_globalReducedNodesIndex.assign(idIn.size(), -1);
}


void NodeFile::createDOFMappingAndCoupling(bool useReduced)
{
    IndexVector *globalDOFIndex;
    Paso_Distribution *dofDistribution;
    if (useReduced) {
        dofDistribution = &m_reducedDegreesOfFreedomDistribution[0];
        globalDOFIndex = &m_globalReducedDOFIndex;
    } else {
        dofDistribution = &m_degreesOfFreedomDistribution[0];
        globalDOFIndex = &m_globalDegreesOfFreedom;
    }
    index_t myFirstDOF = Paso_Distribution_getFirstComponent(dofDistribution);
    index_t myLastDOF = Paso_Distribution_getLastComponent(dofDistribution);
    index_t minDOF, maxDOF;

    if (globalDOFIndex->size() == 0) {
        minDOF = myFirstDOF;
        maxDOF = myLastDOF - 1;
    } else {
        IndexPair minMaxDOF = getFlaggedMinMax(*globalDOFIndex, -1);
        minDOF = minMaxDOF.first;
        maxDOF = minMaxDOF.second;
    }

    const dim_t UNUSED=-1;
    Esys_MPI_rank p;
    dim_t mpiSize = m_mpiInfo->size;
    Esys_MPI_rank p_min = mpiSize;
    Esys_MPI_rank p_max = -1;
    if (maxDOF >= minDOF) {
        for (p = 0; p < mpiSize; ++p) {
            if (dofDistribution->first_component[p] <= minDOF)
                p_min = p;
            if (dofDistribution->first_component[p] <= maxDOF)
                p_max = p;
        }
    }

    if (!((minDOF <= myFirstDOF) && (myLastDOF - 1 <= maxDOF)))
        throw RipleyException("Local elements do not span local degrees of freedom");

    const dim_t len_loc_dof = maxDOF - minDOF + 1;
    const dim_t numNodes = getNumNodes();
    IndexVector wanted_DOFs(numNodes);
    IndexVector nodeMask(numNodes, UNUSED);
    IndexVector shared(numNodes * (p_max - p_min + 1));
    IndexVector offsetInShared(mpiSize + 1);
    IndexVector locDOFMask(len_loc_dof, UNUSED);

#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (index_t i = 0; i < numNodes; ++i) {
            const index_t k = (*globalDOFIndex)[i];
            if (k > -1) {
                locDOFMask[k - minDOF] = UNUSED - 1;
#ifdef BOUNDS_CHECK
                if ((k - minDOF) >= len_loc_dof) {
                    printf("BOUNDS_CHECK %s %d i=%d k=%d minDOF=%d\n", __FILE__, __LINE__, i, k, minDOF);
                    exit(1);
                }
#endif
            }
        }

#pragma omp for schedule(static)
        for (index_t i = myFirstDOF-minDOF; i < myLastDOF-minDOF; ++i) {
            locDOFMask[i] = i - myFirstDOF + minDOF;
#ifdef BOUNDS_CHECK
            if (i < 0 || i >= len_loc_dof) {
                printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i);
                exit(1);
            }
#endif
        }
    }

    RankVector neighbor(mpiSize);
    dim_t numNeighbors = 0, n = 0, lastn = 0;
    vector<dim_t> rcv_len(mpiSize, 0);
    for (p = p_min; p <= p_max; ++p) {
        index_t firstDOF = MAX(minDOF, dofDistribution->first_component[p]);
        index_t lastDOF = MIN(maxDOF+1, dofDistribution->first_component[p+1]);
        if (p != m_mpiInfo->rank) {
            for (index_t i = firstDOF-minDOF; i < lastDOF-minDOF; ++i) {
#ifdef BOUNDS_CHECK
                if (i < 0 || i >= len_loc_dof) {
                    printf("BOUNDS_CHECK %s %d p=%d i=%d\n", __FILE__, __LINE__, p, i);
                    exit(1);
                }
#endif
                if (locDOFMask[i] == UNUSED - 1) {
                    locDOFMask[i] = myLastDOF - myFirstDOF + n;
                    wanted_DOFs[n] = i + minDOF;
                    ++n;
                }
            }
            if (n > lastn) {
                rcv_len[p] = n - lastn;
                neighbor[numNeighbors] = p;
#ifdef BOUNDS_CHECK
                if (numNeighbors < 0 || numNeighbors >= mpiSize + 1) {
                    printf("BOUNDS_CHECK %s %d p=%d numNeighbors=%d n=%d\n", __FILE__, __LINE__, p, numNeighbors,
                           n);
                    exit(1);
                }
#endif
                offsetInShared[numNeighbors] = lastn;
                numNeighbors++;
                lastn = n;
            }
        }
    }
#ifdef BOUNDS_CHECK
    if (numNeighbors < 0 || numNeighbors >= mpiSize + 1) {
        printf("BOUNDS_CHECK %s %d numNeighbors=%d\n", __FILE__, __LINE__, numNeighbors);
        exit(1);
    }
#endif
    offsetInShared[numNeighbors] = lastn;

    // assign new DOF labels to nodes
#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < numNodes; ++i) {
        const index_t k = (*globalDOFIndex)[i];
        if (k > -1)
            nodeMask[i] = locDOFMask[k - minDOF];
    }

#ifdef BOUNDS_CHECK
    for (index_t i = 0; i < offsetInShared[numNeighbors]; ++i) {
        if (i < 0 || i >= numNodes * (p_max - p_min + 1)) {
            printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i);
            exit(1);
        }
    }
#endif

    // define how to get DOF values for controlled but other processors
#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < offsetInShared[numNeighbors]; ++i)
        shared[i] = myLastDOF - myFirstDOF + i;

    Paso_SharedComponents *rcv_shcomp = Paso_SharedComponents_alloc(
            myLastDOF - myFirstDOF, numNeighbors, &neighbor[0],
            &shared[0], &offsetInShared[0], 1, 0, m_mpiInfo);

     // now we build the sender
    vector<dim_t> snd_len(mpiSize);
#ifdef ESYS_MPI
    vector<MPI_Request> mpiRequests(mpiSize*2);
    MPI_Alltoall(&rcv_len[0], 1, MPI_INT, &snd_len[0], 1, MPI_INT, m_mpiInfo->comm);
#else
    for (p = 0; p < mpiSize; ++p)
        snd_len[p] = rcv_len[p];
#endif
    dim_t count = 0;
    for (p = 0; p < rcv_shcomp->numNeighbors; p++) {
#ifdef ESYS_MPI
        MPI_Isend(&(wanted_DOFs[rcv_shcomp->offsetInShared[p]]),
                rcv_shcomp->offsetInShared[p + 1] - rcv_shcomp->offsetInShared[p],
                MPI_INT, rcv_shcomp->neighbor[p], m_mpiInfo->msg_tag_counter + m_mpiInfo->rank,
                m_mpiInfo->comm, &mpiRequests[count]);
#endif
        count++;
    }
    n = 0;
    numNeighbors = 0;
    for (p = 0; p < mpiSize; p++) {
        if (snd_len[p] > 0) {
#ifdef ESYS_MPI
            MPI_Irecv(&(shared[n]), snd_len[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + p, m_mpiInfo->comm,
                    &mpiRequests[count]);
#endif
            count++;
            neighbor[numNeighbors] = p;
            offsetInShared[numNeighbors] = n;
            numNeighbors++;
            n += snd_len[p];
        }
    }
    m_mpiInfo->msg_tag_counter += mpiSize;
    offsetInShared[numNeighbors] = n;
#ifdef ESYS_MPI
    vector<MPI_Status> mpiStati(mpiSize*2);
    MPI_Waitall(count, &mpiRequests[0], &mpiStati[0]);
#endif
    // map global ids to local ids
#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < offsetInShared[numNeighbors]; ++i)
        shared[i] = locDOFMask[shared[i] - minDOF];

    Paso_SharedComponents *snd_shcomp = Paso_SharedComponents_alloc(
            myLastDOF - myFirstDOF, numNeighbors, &neighbor[0],
            &shared[0], &offsetInShared[0], 1, 0, m_mpiInfo);

    Paso_Connector *connector = Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
    Paso_SharedComponents_free(rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);

    if (useReduced) {
        m_reducedDegreesOfFreedomMapping.reset(new NodeMapping(nodeMask, UNUSED));
        m_reducedDegreesOfFreedomConnector = connector;
    } else {
        m_degreesOfFreedomMapping.reset(new NodeMapping(nodeMask, UNUSED));
        m_degreesOfFreedomConnector = connector;
    }
}

void NodeFile::createNodeFileMappings(const IndexVector &indexReducedNodes,
        const IndexVector &dof_first_component,
        const IndexVector &nodes_first_component)
{
    index_t myFirstDOF, myLastDOF, myFirstNode, myLastNode;
    dim_t myNumDOF, myNumNodes, globalNumReducedNodes, globalNumReducedDOF, i;
    const dim_t numReducedNodes = indexReducedNodes.size();
    const dim_t UNUSED=-1;

    // mark the nodes used by the reduced mesh

    myFirstDOF = dof_first_component[m_mpiInfo->rank];
    myLastDOF = dof_first_component[m_mpiInfo->rank + 1];
    myNumDOF = myLastDOF - myFirstDOF;

    myFirstNode = nodes_first_component[m_mpiInfo->rank];
    myLastNode = nodes_first_component[m_mpiInfo->rank + 1];
    myNumNodes = myLastNode - myFirstNode;

    IndexVector maskMyReducedDOF(myNumDOF, -1);
    IndexVector maskMyReducedNodes(myNumNodes, -1);

#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < numReducedNodes; ++i) {
        index_t k = m_globalNodesIndex[indexReducedNodes[i]];
        if ((k >= myFirstNode) && (myLastNode > k))
            maskMyReducedNodes[k - myFirstNode] = i;
        k = m_globalDegreesOfFreedom[indexReducedNodes[i]];
        if ((k >= myFirstDOF) && (myLastDOF > k))
            maskMyReducedDOF[k - myFirstDOF] = i;
    }
    const IndexVector indexMyReducedNodes = packMask(maskMyReducedNodes);
    const IndexVector indexMyReducedDOF = packMask(maskMyReducedDOF);
    IndexVector reduced_dof_first_component(m_mpiInfo->size + 1);
    IndexVector reduced_nodes_first_component(m_mpiInfo->size + 1);

#ifdef ESYS_MPI
    dim_t myNumReducedNodes = indexMyReducedNodes.size();
    dim_t myNumReducedDOF = indexMyReducedDOF.size();
    MPI_Allgather(&myNumReducedNodes, 1, MPI_INT,
            &reduced_nodes_first_component[0], 1, MPI_INT, m_mpiInfo->comm);
    MPI_Allgather(&myNumReducedDOF, 1, MPI_INT,
            &reduced_dof_first_component[0], 1, MPI_INT, m_mpiInfo->comm);
#else
    reduced_nodes_first_component[0] = indexMyReducedNodes.size();
    reduced_dof_first_component[0] = indexMyReducedDOF.size();
#endif
    globalNumReducedNodes = 0;
    globalNumReducedDOF = 0;
    for (i = 0; i < m_mpiInfo->size; ++i) {
        index_t k = reduced_nodes_first_component[i];
        reduced_nodes_first_component[i] = globalNumReducedNodes;
        globalNumReducedNodes += k;

        k = reduced_dof_first_component[i];
        reduced_dof_first_component[i] = globalNumReducedDOF;
        globalNumReducedDOF += k;
    }
    reduced_nodes_first_component[m_mpiInfo->size] = globalNumReducedNodes;
    reduced_dof_first_component[m_mpiInfo->size] = globalNumReducedDOF;
    // distribution of nodes
    const index_t *ptr = &nodes_first_component[0];
    m_nodesDistribution = Paso_Distribution_alloc(m_mpiInfo, const_cast<index_t *>(ptr), 1, 0);

    // distribution of DOFs
    ptr = &dof_first_component[0];
    m_degreesOfFreedomDistribution = Paso_Distribution_alloc(m_mpiInfo,
            const_cast<index_t *>(ptr), 1, 0);

    // distribution of reduced nodes
    ptr = &reduced_nodes_first_component[0];
    m_reducedNodesDistribution = Paso_Distribution_alloc( m_mpiInfo,
            const_cast<index_t *>(ptr), 1, 0);

    // distribution of reduced DOFs
    ptr = &reduced_dof_first_component[0];
    m_reducedDegreesOfFreedomDistribution = Paso_Distribution_alloc(m_mpiInfo,
            const_cast<index_t *>(ptr), 1, 0);

    // dummy node mapping
    IndexVector nodeMask(getNumNodes());
#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < getNumNodes(); ++i)
        nodeMask[i] = i;
    m_nodesMapping.reset(new NodeMapping(nodeMask, UNUSED));

    // mapping between nodes and reduced nodes
#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < getNumNodes(); ++i)
        nodeMask[i] = UNUSED;
#pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < numReducedNodes; ++i)
        nodeMask[indexReducedNodes[i]] = i;

    m_reducedNodesMapping.reset(new NodeMapping(nodeMask, UNUSED));

    // mapping between nodes and DOFs + DOF connector
    createDOFMappingAndCoupling(false);
    // mapping between nodes and reduced DOFs + reduced DOF connector
    createDOFMappingAndCoupling(true);

    // get the Ids for DOFs and reduced nodes
    m_reducedNodesId.resize(m_reducedNodesMapping->numTargets);
    m_degreesOfFreedomId.resize(m_degreesOfFreedomMapping->numTargets);
    m_reducedDegreesOfFreedomId.resize(m_reducedDegreesOfFreedomMapping->numTargets);
#pragma omp parallel private(i)
    {
#pragma omp for
        for (i = 0; i < m_reducedNodesMapping->numTargets; ++i)
            m_reducedNodesId[i] = m_id[m_reducedNodesMapping->map[i]];
#pragma omp for
        for (i = 0; i < m_degreesOfFreedomMapping->numTargets; ++i)
            m_degreesOfFreedomId[i] = m_id[m_degreesOfFreedomMapping->map[i]];
#pragma omp for
        for (i = 0; i < m_reducedDegreesOfFreedomMapping->numTargets; ++i)
            m_reducedDegreesOfFreedomId[i] =
                m_id[m_reducedDegreesOfFreedomMapping->map[i]];
    }
}

} // end of namespace ripley

