
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

#include <ripley/ElementFile.h>
#include <ripley/RipleyException.h>
#include <ripley/Util.h>

#include <escript/Data.h>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif


using namespace std;

namespace ripley {

ElementTypeId ElementFile::ripleyTypeFromString(const string &s)
{
    if (s == "Point1")
        return Point1;
    if (s == "Line2")
        return Line2;
    if (s == "Rec4")
        return Rec4;
    if (s == "Hex8")
        return Hex8;
    if (s == "Line2Face")
        return Line2Face;
    if (s == "Rec4Face")
        return Rec4Face;
    if (s == "Hex8Face")
        return Hex8Face;
    return InvalidElementType;
}

ElementFile::ElementFile(ElementTypeId type, Esys_MPIInfo* mpiInfo) :
    m_type(type),
    m_minColor(0),
    m_maxColor(-1)
{
    m_mpiInfo = Esys_MPIInfo_getReference(mpiInfo);
    switch (m_type) {
        case Hex8:
            m_numDim=m_localDim=3;
            break;
        case Rec4:
            m_numDim=m_localDim=2;
            break;
        case Line2:
            m_numDim=m_localDim=1;
            break;
        case Hex8Face:
            m_numDim=3;
            m_localDim=2;
            break;
        case Rec4Face:
            m_numDim=2;
            m_localDim=1;
            break;
        case Line2Face:
            m_numDim=1;
            m_localDim=0;
            break;
        default:
            m_numDim=m_localDim=0;
    }
}

ElementFile::~ElementFile()
{
    Esys_MPIInfo_free(m_mpiInfo);
}

string ElementFile::getName() const
{
    switch (m_type) {
        case Point1:
            return "Point1";
        case Line2:
            return "Line2";
        case Rec4:
            return "Rec4";
        case Hex8:
            return "Hex8";
        case Line2Face:
            return "Line2Face";
        case Rec4Face:
            return "Rec4Face";
        case Hex8Face:
            return "Hex8Face";
        default:
            return "NoType";
    }
}

#ifdef USE_NETCDF
void ElementFile::readFromNetCDF(NcFile &dataFile, dim_t numElements, const string &name)
{
    m_minColor=0;
    m_maxColor=numElements-1;

    // clear and resize data vectors
    m_id.assign(numElements, -1);
    m_tag.assign(numElements, -1);
    m_owner.assign(numElements, -1);
    m_color.assign(numElements, -1);
    m_nodes.assign(numElements*getNumNodes(), -1);

    NcVar *ncVar;
    const string msgPrefix("read: netCDF operation failed - ");

    if (numElements > 0) {
        // Id
        if (! (ncVar = dataFile.get_var((name+"_Id").c_str())) )
            throw RipleyException(msgPrefix+"get_var("+name+"_Id)");
        if (!ncVar->get(&m_id[0], numElements))
            throw RipleyException(msgPrefix+"get("+name+"_Id)");
        // Tag
        if (! (ncVar = dataFile.get_var((name+"_Tag").c_str())) )
            throw RipleyException(msgPrefix+"get_var("+name+"_Tag)");
        if (!ncVar->get(&m_tag[0], numElements))
            throw RipleyException(msgPrefix+"get("+name+"_Tag)");
        // Owner
        if (! (ncVar = dataFile.get_var((name+"_Owner").c_str())) )
            throw RipleyException(msgPrefix+"get_var("+name+"_Owner)");
        if (!ncVar->get(&m_owner[0], numElements))
            throw RipleyException(msgPrefix+"get("+name+"_Owner)");
        // Color
        if (! (ncVar = dataFile.get_var((name+"_Color").c_str())) )
            throw RipleyException(msgPrefix+"get_var("+name+"_Color)");
        if (! ncVar->get(&m_color[0], numElements) )
            throw RipleyException(msgPrefix+"get("+name+"_Color)");
        // Nodes
        if (!(ncVar = dataFile.get_var((name+"_Nodes").c_str())))
            throw RipleyException(msgPrefix+"get_var("+name+"_Nodes)");

        if (!ncVar->get(&(m_nodes[0]), numElements, getNumNodes()))
            throw RipleyException(msgPrefix+"get("+name+"_Nodes)");
    } // numElements > 0
}

void ElementFile::dumpToNetCDF(NcFile &dataFile, const string &name)
{
    vector<NcDim*> ncdims(2);
    NcVar *ncVar;
    const string msgPrefix("dump: netCDF operation failed - ");
    const dim_t numElements = getNumElements();

    if (!dataFile.add_att(("num_"+name).c_str(), numElements))
        throw RipleyException(msgPrefix+"add_att(num_"+name+")");
    if (!dataFile.add_att(("num_"+name+"_numNodes").c_str(), getNumNodes()))
        throw RipleyException(msgPrefix+"add_att(num_"+name+"_numNodes)");
    if (!dataFile.add_att((name+"_TypeId").c_str(), m_type))
        throw RipleyException(msgPrefix+"add_att("+name+"_TypeId)");

    if (numElements > 0) {
        if (! (ncdims[0] = dataFile.add_dim(("dim_"+name).c_str(), numElements)) )
            throw RipleyException(msgPrefix+"add_dim(dim_Elements)");
        if (! (ncdims[1] = dataFile.add_dim(("dim_"+name+"_Nodes").c_str(), getNumNodes())) )
            throw RipleyException(msgPrefix+"add_dim(dim_"+name+"_Nodes)");

        // Id
        if (! (ncVar = dataFile.add_var((name+"_Id").c_str(), ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var("+name+"_Id)");
        if (!ncVar->put(&m_id[0], numElements))
            throw RipleyException(msgPrefix+"put("+name+"_Id)");

        // Tag
        if (! (ncVar = dataFile.add_var((name+"_Tag").c_str(), ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var("+name+"_Tag)");
        if (!ncVar->put(&m_tag[0], numElements))
            throw RipleyException(msgPrefix+"put("+name+"_Tag)");

        // Owner
        if (! (ncVar = dataFile.add_var((name+"_Owner").c_str(), ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var("+name+"_Owner)");
        if (!ncVar->put(&m_owner[0], numElements))
            throw RipleyException(msgPrefix+"put("+name+"_Owner)");

        // Color
        if (! (ncVar = dataFile.add_var((name+"_Color").c_str(), ncInt, ncdims[0])) )
            throw RipleyException(msgPrefix+"add_var("+name+"_Color)");
        if (!ncVar->put(&m_color[0], numElements))
            throw RipleyException(msgPrefix+"put("+name+"_Color)");

        // Nodes
        if (! (ncVar = dataFile.add_var((name+"_Nodes").c_str(), ncInt, ncdims[0], ncdims[1]) ) )
            throw RipleyException(msgPrefix+"add_var("+name+"_Nodes)");
        if (!ncVar->put(&m_nodes[0], numElements, getNumNodes()))
            throw RipleyException(msgPrefix+"put("+name+"_Nodes)");
    }
}
#endif

void ElementFile::updateTagsInUse()
{
    m_tagsInUse = getUniqueValues(m_tag, m_mpiInfo);
}

void ElementFile::setTags(int newTag, const escript::Data &mask)
{
    dim_t numElements = getNumElements();
    register dim_t n;
    escriptDataC maskC = mask.getDataC();

    if (1 != getDataPointSize(&maskC))
        throw RipleyException("setTags: Number of components in mask must be 1");
    if (!numSamplesEqual(&maskC, 1, numElements))
        throw RipleyException("setTags: Illegal number of samples in mask Data object");

    dim_t numQuad;
    const int fs = getFunctionSpaceType(&maskC);
    if (fs == ReducedElements || fs == ReducedFaceElements)
        numQuad = 1;
    else
        numQuad = m_numDim+1;

    if (isExpanded(&maskC) || numQuad==1) {
#pragma omp parallel for schedule(static) private(n)
        for (n = 0; n < numElements; n++) {
            const double *maskArray = getSampleDataRO(&maskC, n);
            if (maskArray[0] > 0)
                m_tag[n] = newTag;
        }
    } else {
#pragma omp parallel for schedule(static) private(n)
        for (n = 0; n < numElements; n++) {
            const double *maskArray = getSampleDataRO(&maskC, n);
            bool setTag=false;
            for (dim_t q=0; q<numQuad; q++)
                setTag = setTag || maskArray[q];
            if (setTag)
                m_tag[n] = newTag;
        }
    }
    updateTagsInUse();
}

void ElementFile::swapEntries(IndexVector &idIn, IndexVector &tagIn,
        IndexVector &colorIn, RankVector &ownerIn, IndexVector &nodesIn)
{
    m_id.swap(idIn);
    m_tag.swap(tagIn);
    m_color.swap(colorIn);
    m_owner.swap(ownerIn);
    m_nodes.swap(nodesIn);
}

void ElementFile::relabelNodes(const IndexVector &newNode, index_t offset)
{
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();
#pragma omp parallel for schedule(static)
    for (dim_t j=0; j<NE; j++)
        for (dim_t i=0; i<NN; i++)
            m_nodes[INDEX2(i, j, NN)] = newNode[m_nodes[INDEX2(i, j, NN)]-offset];
}

void ElementFile::markNodes(IndexVector &mask, index_t offset)
{
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();
#pragma omp parallel for schedule(static)
    for (dim_t j=0; j<NE; j++)
        for (dim_t i=0; i<NN; i++)
            mask[m_nodes[INDEX2(i, j, NN)]-offset] = 1;
}

IndexPair ElementFile::getNodeRange() const
{
    return getMinMax(m_nodes);
}

void ElementFile::createColoring(const IndexVector &degreesOfFreedom)
{
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();
    if (NE == 0)
        return;

    IndexPair minMaxId = getMinMax(degreesOfFreedom);
    const dim_t len = minMaxId.second - minMaxId.first + 1;

#pragma omp parallel for schedule(static)
    for (dim_t e = 0; e < NE; e++)
        m_color[e] = -1;

    dim_t numUncoloredElements = NE;
    m_minColor = 0;
    m_maxColor = -1;
    while (numUncoloredElements > 0) {
        IndexVector maskDOF(len, -1);
        numUncoloredElements = 0;
        /* TODO: OMP ? */
        for (dim_t e = 0; e < NE; e++) {
            if (m_color[e] < 0) {
                // find out if element e is independent from the elements
                // already colored
                bool independent = true;
                for (dim_t i = 0; i < NN; i++) {
                    const index_t maskIdx = degreesOfFreedom[m_nodes[INDEX2(i,e,NN)]]-minMaxId.first;
#ifdef BOUNDS_CHECK
                    if (m_nodes[INDEX2(i, e, NN)] < 0 || m_nodes[INDEX2(i, e, NN)] >= degreesOfFreedom.size()) {
                        printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d minId=%d m_nodes[INDEX2...]=%d\n",
                               __FILE__, __LINE__, i, e, NN, minMaxId.first,
                               m_nodes[INDEX2(i, e, NN)]);
                        exit(1);
                    }
                    if (maskIdx >= len || maskIdx < 0) {
                        printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d minId=%d maskIdx=%d\n",
                               __FILE__, __LINE__, i, e, NN, minMaxId.first, maskIdx);
                        exit(1);
                    }
#endif
                    if (maskDOF[maskIdx] > 0) {
                        independent = false;
                        break;
                    }
                }
                // if e is independent a new color is assigned and the
                // nodes are marked as being used
                if (independent) {
                    for (dim_t i = 0; i < NN; i++)
                        maskDOF[degreesOfFreedom[m_nodes[INDEX2(i, e, NN)]] - minMaxId.first] = 1;
                    m_color[e] = m_maxColor + 1;
                } else {
                    numUncoloredElements++;
                }
            }
        }
        m_maxColor++;
    }
}

void ElementFile::distributeByRankOfDOF(const RankVector &mpiRankOfDOF, const IndexVector &id)
{
    const Esys_MPI_rank myRank = m_mpiInfo->rank;
    const dim_t size = m_mpiInfo->size;
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();

    if (size == 1) {
#pragma omp for schedule(static)
        for (dim_t e = 0; e < NE; e++) {
            m_owner[e] = myRank;
            for (dim_t i = 0; i < NN; i++)
                m_nodes[INDEX2(i, e, NN)] = id[m_nodes[INDEX2(i, e, NN)]];
        }
        return;
    }

    // size > 1
#ifdef ESYS_MPI
    // count the number of elements that have to be sent to each processor
    // (sendCount) and define a new element owner as the processor with
    // the largest number of DOFs and the smallest ID
    vector<dim_t> sendCount(size, 0);
    vector<dim_t> recvCount(size);
    RankVector newOwner(NE);
#pragma omp parallel
    {
        vector<dim_t> locSendCount(size, 0);
#pragma omp for schedule(static)
        for (dim_t e = 0; e < NE; e++) {
            if (m_owner[e] == myRank) {
                newOwner[e] = myRank;
                vector<dim_t> locProcMask(size, 0);
                for (dim_t j = 0; j < NN; j++) {
                    Esys_MPI_rank p = mpiRankOfDOF[m_nodes[INDEX2(j, e, NN)]];
                    locProcMask[p]++;
                }
                dim_t locProcMaskMax = 0;
                for (Esys_MPI_rank p = 0; p < size; ++p) {
                    if (locProcMask[p] > 0)
                        locSendCount[p]++;
                    if (locProcMask[p] > locProcMaskMax) {
                        newOwner[e] = p;
                        locProcMaskMax = locProcMask[p];
                    }
                }
            } else {
                newOwner[e] = -1;
            }
        }
#pragma omp critical
        for (Esys_MPI_rank p = 0; p < size; ++p)
            sendCount[p] += locSendCount[p];
    }
    MPI_Alltoall(&sendCount[0], 1, MPI_INT, &recvCount[0], 1, MPI_INT,
                 m_mpiInfo->comm);
    // get the new number of elements for this processor
    dim_t newNumElements = 0;
    for (Esys_MPI_rank p = 0; p < size; ++p)
        newNumElements += recvCount[p];

    dim_t numElementsInBuffer = 0;
    for (Esys_MPI_rank p = 0; p < size; ++p)
        numElementsInBuffer += sendCount[p];

    // allocate buffers
    IndexVector idBuffer(numElementsInBuffer);
    IndexVector tagBuffer(numElementsInBuffer);
    RankVector ownerBuffer(numElementsInBuffer);
    IndexVector nodesBuffer(numElementsInBuffer * NN);
    IndexVector sendOffset(size);
    IndexVector recvOffset(size);

    /* calculate the offsets for the processor buffers */
    recvOffset[0] = 0;
    sendOffset[0] = 0;
    for (Esys_MPI_rank p = 0; p < size - 1; ++p) {
        recvOffset[p + 1] = recvOffset[p] + recvCount[p];
        sendOffset[p + 1] = sendOffset[p] + sendCount[p];
    }

    sendCount.assign(size, 0);
    // copy elements into buffers. procMask makes sure that an
    // element is copied once only for each processor
    for (dim_t e = 0; e < NE; e++) {
        if (m_owner[e] == myRank) {
            vector<bool> procMask(size, true);
            for (dim_t j = 0; j < NN; j++) {
                Esys_MPI_rank p = mpiRankOfDOF[m_nodes[INDEX2(j, e, NN)]];
                if (procMask[p]) {
                    const index_t k = sendOffset[p] + sendCount[p];
                    idBuffer[k] = m_id[e];
                    tagBuffer[k] = m_tag[e];
                    ownerBuffer[k] = newOwner[e];
                    for (dim_t i = 0; i < NN; i++)
                        nodesBuffer[INDEX2(i, k, NN)] = id[m_nodes[INDEX2(i, e, NN)]];
                    sendCount[p]++;
                    procMask[p] = false;
                }
            }
        }
    }
    // allocate new tables
    m_id.resize(newNumElements);
    m_tag.resize(newNumElements);
    m_color.assign(newNumElements, -1);
    m_owner.resize(newNumElements);
    m_nodes.resize(newNumElements * NN);

    // start to receive new elements
    dim_t numRequests = 0;
    vector<MPI_Request> mpiRequests(8*size);
    for (Esys_MPI_rank p = 0; p < size; ++p) {
        if (recvCount[p] > 0) {
            MPI_Irecv(&(m_id[recvOffset[p]]), recvCount[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + myRank,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Irecv(&(m_tag[recvOffset[p]]), recvCount[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + size + myRank,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Irecv(&(m_owner[recvOffset[p]]), recvCount[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + 2*size + myRank,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Irecv(&(m_nodes[recvOffset[p]*NN]), recvCount[p]*NN, MPI_INT,
                    p, m_mpiInfo->msg_tag_counter + 3*size + myRank,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
        }
    }
    // now the buffers can be sent away
    for (Esys_MPI_rank p = 0; p < size; ++p) {
        if (sendCount[p] > 0) {
            MPI_Issend(&(idBuffer[sendOffset[p]]), sendCount[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + p,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Issend(&(tagBuffer[sendOffset[p]]), sendCount[p], MPI_INT, p,
                    m_mpiInfo->msg_tag_counter + size + p,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Issend(&(ownerBuffer[sendOffset[p]]), sendCount[p], MPI_INT,
                    p, m_mpiInfo->msg_tag_counter + 2*size + p,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
            MPI_Issend(&(nodesBuffer[sendOffset[p]*NN]), sendCount[p]*NN,
                    MPI_INT, p, m_mpiInfo->msg_tag_counter + 3*size + p,
                    m_mpiInfo->comm, &mpiRequests[numRequests]);
            numRequests++;
        }
    }
    m_mpiInfo->msg_tag_counter += 4*size;
    // wait for the requests to be finalized
    vector<MPI_Status> mpiStati(8*size);
    MPI_Waitall(numRequests, &mpiRequests[0], &mpiStati[0]);
#endif
}

static bool compareIndexPair(const IndexPair &a, const IndexPair &b)
{
    if (a.second < b.second)
        return true;
    if (a.second > b.second)
        return false;
    if (a.first < b.first)
        return true;
    if (a.first > b.first)
        return false;
    return true;
}

void ElementFile::optimizeOrdering()
{
    const dim_t NE = getNumElements();
    if (NE > 0) {
        const dim_t NN = getNumNodes();
        vector<IndexPair> itemList(NE);
#pragma omp parallel for schedule(static)
        for (dim_t e = 0; e < NE; e++) {
            index_t n = m_nodes[INDEX2(0, e, NN)];
            for (dim_t i = 1; i < NN; i++)
                n = MIN(n, m_nodes[INDEX2(i, e, NN)]);
            itemList[e] = IndexPair(e, n);
        }
        sort(itemList.begin(), itemList.end(), compareIndexPair);

        // now reorder vectors
        IndexVector newId(NE);
        IndexVector newTag(NE);
        IndexVector newColor(NE);
        IndexVector newNodes(NE*NN);
        RankVector newOwner(NE);
#pragma omp parallel for schedule(static)
        for (dim_t e = 0; e < NE; e++) {
            const index_t k = itemList[e].first;
            newId[e] = m_id[k];
            newTag[e] = m_tag[k];
            newOwner[e] = m_owner[k];
            newColor[e] = m_color[k];
            for (dim_t j=0; j<NN; j++)
                newNodes[INDEX2(j, e, NN)] = m_nodes[INDEX2(j, k, NN)];
        }
        swapEntries(newId, newTag, newColor, newOwner, newNodes);
    }
}

void ElementFile::insertIntoIndexMatrix(IndexMatrix &matrix,
        const IndexVector &rowMap, const IndexVector &colMap,
        index_t firstRow, index_t lastRow)
{
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();

    if (firstRow >=0 && lastRow >= 0) {
        for (index_t color = m_minColor; color <= m_maxColor; color++) {
#pragma omp parallel for schedule(static)
            for (dim_t e = 0; e < NE; e++) {
                if (m_color[e] == color) {
                    for (dim_t kr = 0; kr < NN; kr++) {
                        const dim_t irow = rowMap[m_nodes[INDEX2(kr, e, NN)]];
                        if ((firstRow <= irow) && (irow < lastRow)) {
                            for (dim_t kc = 0; kc < NN; kc++) {
                                const dim_t icol = colMap[m_nodes[INDEX2(kc, e, NN)]];
                                matrix[irow-firstRow].insert(icol);
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (index_t color = m_minColor; color <= m_maxColor; color++) {
#pragma omp parallel for schedule(static)
            for (dim_t e = 0; e < NE; e++) {
                if (m_color[e] == color) {
                    for (dim_t kr = 0; kr < NN; kr++) {
                        const dim_t irow = rowMap[m_nodes[INDEX2(kr, e, NN)]];
                        for (dim_t kc = 0; kc < NN; kc++) {
                            const dim_t icol = colMap[m_nodes[INDEX2(kc, e, NN)]];
                            matrix[irow].insert(icol);
                        }
                    }
                }
            }
        }
    }
}

void ElementFile::insertIntoIndexMatrixNoMainDiagonal(IndexMatrix &matrix,
        const IndexVector &rowMap, const IndexVector &colMap,
        index_t firstRow, index_t lastRow)
{
    const dim_t NE = getNumElements();
    const dim_t NN = getNumNodes();

    for (index_t color = m_minColor; color <= m_maxColor; color++) {
#pragma omp parallel for schedule(static)
        for (dim_t e = 0; e < NE; e++) {
            if (m_color[e] == color) {
                for (dim_t kr = 0; kr < NN; kr++) {
                    const dim_t irow = rowMap[m_nodes[INDEX2(kr, e, NN)]];
                    if ((firstRow <= irow) && (irow < lastRow)) {
                        for (dim_t kc = 0; kc < NN; kc++) {
                            const dim_t icol = colMap[m_nodes[INDEX2(kc, e, NN)]];
                            if (icol != irow)
                                matrix[irow-firstRow].insert(icol);
                        }
                    }
                }
            }
        }
    }
}


} // end of namespace ripley

