
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

#include <escript/index.h>

using escript::IOError;

namespace {

using namespace finley;

ElementFile* readElementFile(std::ifstream& fileHandle, int order,
                             int reducedOrder, escript::JMPI mpiInfo)
{
    dim_t numEle = 0;
    ElementTypeId typeID = NoRef;
    std::string elementType, line;

    // Read the element typeID and number of elements
    if (mpiInfo->rank == 0) {
        std::getline(fileHandle, line);
        if (!fileHandle.good())
            throw IOError("Mesh::read: Scan error while reading file - expected <ElementType> <numEle>");
        size_t pos = line.find(' ');
        if (pos == std::string::npos)
            throw IOError("Mesh::read: Scan error reading file - expected <ElementType> <numEle>");
        elementType = line.substr(0, pos);
        numEle = std::stol(line.substr(pos+1));
        typeID = ReferenceElement::getTypeId(elementType.c_str());
    }
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        dim_t temp1[2];
        temp1[0] = (dim_t)typeID;
        temp1[1] = numEle;
        int mpiError = MPI_Bcast(temp1, 2, MPI_DIM_T, 0, mpiInfo->comm);
        if (mpiError != MPI_SUCCESS) {
            throw FinleyException("Mesh::read: broadcast of element typeID failed");
        }
        typeID = static_cast<ElementTypeId>(temp1[0]);
        numEle = temp1[1];
    }
#endif
    if (typeID == NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << elementType;
        throw IOError(ss.str());
    }

    // Allocate the ElementFile
    const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(
                                                typeID, order, reducedOrder));
    ElementFile* out = new ElementFile(refElements, mpiInfo);
    const int numNodes = out->numNodes;

    /********************** Read the element data **************************/
    dim_t chunkSize = numEle / mpiInfo->size + 1;
    dim_t totalEle = 0;
    dim_t chunkEle = 0;
    int nextCPU = 1;
    /// Store Id + Tag + node list (+ one int at end for chunkEle)
    index_t* tempInts = new index_t[chunkSize * (2 + numNodes) + 1];
    // Elements are specified as a list of integers...only need one message
    // instead of two as with the nodes
    if (mpiInfo->rank == 0) { // Master
        for (;;) {            // Infinite loop
#pragma omp parallel for
            for (index_t i0 = 0; i0 < chunkSize * (2 + numNodes) + 1; i0++)
                tempInts[i0] = -1;
            chunkEle = 0;
            for (index_t i0 = 0; i0 < chunkSize; i0++) {
                if (totalEle >= numEle)
                    break; // End inner loop
                std::getline(fileHandle, line);
                if (!fileHandle.good())
                    throw IOError("Mesh::read: Scan error while reading element data");
                std::stringstream ss;
                ss << line;
                ss >> tempInts[i0 * (2 + numNodes) + 0]
                   >> tempInts[i0 * (2 + numNodes) + 1];
                for (int i1 = 0; i1 < numNodes; i1++) {
                    ss >> tempInts[i0 * (2 + numNodes) + 2 + i1];
                }
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunk of elements to each CPU except 0
            // itself, here goes one of them
            if (nextCPU < mpiInfo->size) {
                tempInts[chunkSize * (2 + numNodes)] = chunkEle;
                MPI_Send(tempInts, chunkSize * (2 + numNodes) + 1, MPI_DIM_T,
                         nextCPU, 81722, mpiInfo->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of the worker
            // nodes plus one more chunk for the master
            if (nextCPU > mpiInfo->size)
                break; // End infinite loop
        } // Infinite loop
    } // end master
    else { // Worker
#ifdef ESYS_MPI
        // Each worker receives one message
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize * (2 + numNodes) + 1, MPI_DIM_T, 0,
                 81722, mpiInfo->comm, &status);
        chunkEle = tempInts[chunkSize * (2 + numNodes)];
#endif
    } // Worker

    out->allocTable(chunkEle);

    // Copy Element data from tempInts to element file
    out->minColor = 0;
    out->maxColor = chunkEle - 1;
#pragma omp parallel for
    for (index_t i0 = 0; i0 < chunkEle; i0++) {
        out->Id[i0] = tempInts[i0 * (2 + numNodes) + 0];
        out->Tag[i0] = tempInts[i0 * (2 + numNodes) + 1];
        out->Owner[i0] = mpiInfo->rank;
        out->Color[i0] = i0;
        for (int i1 = 0; i1 < numNodes; i1++) {
            out->Nodes[INDEX2(i1, i0, numNodes)] =
                tempInts[i0 * (2 + numNodes) + 2 + i1];
        }
    }
    delete[] tempInts;
    return out;
}

} // anonymous

namespace finley {

escript::Domain_ptr FinleyDomain::read(escript::JMPI mpiInfo,
                                       const std::string& filename,
                                       int order, int reducedOrder,
                                       bool optimize)
{
    dim_t numNodes = 0;
    int numDim = 0;
    std::string name, line, token;
    std::ifstream fileHandle;

    if (mpiInfo->rank == 0) {
        // open file
        fileHandle.open(filename.c_str());
        if (!fileHandle.good()) {
            std::stringstream ss;
            ss << "Mesh::read: Opening file " << filename
               << " for reading failed.";
            throw IOError(ss.str());
        }

        // read header
        std::getline(fileHandle, name);
        if (!fileHandle.good())
            throw IOError("Mesh::read: Scan error while reading file header");

        // get the number of dimensions and nodes
        std::getline(fileHandle, line);
        if (!fileHandle.good())
            throw IOError("Mesh::read: Scan error while reading file - expected <?D-Nodes> <numNodes>");
        numDim = std::stoi(line.substr(0, 1));
        token = line.substr(line.find(' ')+1);
        numNodes = std::stoi(token);
    }

#ifdef ESYS_MPI
    // MPI Broadcast numDim, numNodes, name if there are multiple MPI procs
    if (mpiInfo->size > 1) {
        dim_t temp1[3];
        if (mpiInfo->rank == 0) {
            temp1[0] = numDim;
            temp1[1] = numNodes;
            temp1[2] = name.length() + 1;
        } else {
            temp1[0] = 0;
            temp1[1] = 0;
            temp1[2] = 1;
        }
        MPI_Bcast(temp1, 3, MPI_DIM_T, 0, mpiInfo->comm);
        numDim = temp1[0];
        numNodes = temp1[1];
        name.resize(temp1[2]);
        MPI_Bcast(&name[0], temp1[2], MPI_CHAR, 0, mpiInfo->comm);
    }
#endif

    // allocate domain
    FinleyDomain* domain = new FinleyDomain(name, numDim, mpiInfo);

    // Each CPU will get at most chunkSize nodes so the message has to be
    // sufficiently large
    dim_t chunkSize = numNodes / mpiInfo->size + 1;
    dim_t totalNodes = 0;
    dim_t chunkNodes = 0;
    int nextCPU = 1;
    // Stores the integer message data
    index_t* tempInts = new index_t[chunkSize * 3 + 1];
    // Stores the double message data
    double* tempCoords = new double[chunkSize * numDim];

    // Read chunkSize nodes, send it in a chunk to worker CPU which copies
    // chunk into its local domain.  It doesn't matter that a CPU has the wrong
    // nodes for its elements, this is sorted out later. First chunk sent to
    // CPU 1, second to CPU 2, ..., last chunk stays on CPU 0 (the master).
    // The three columns of integers (Id, gDOF, Tag) are gathered into a single
    // array tempInts and sent together in a single MPI message.
    if (mpiInfo->rank == 0) { // Master
        for (;;) {            // Infinite loop
#pragma omp parallel for
            for (index_t i0 = 0; i0 < chunkSize * 3 + 1; i0++)
                tempInts[i0] = -1;

#pragma omp parallel for
            for (index_t i0 = 0; i0 < chunkSize * numDim; i0++)
                tempCoords[i0] = -1.0;

            chunkNodes = 0;
            for (index_t i1 = 0; i1 < chunkSize; i1++) {
                if (totalNodes >= numNodes)
                    break;  // End of inner loop
                std::getline(fileHandle, line);
                if (!fileHandle.good())
                    throw IOError("Mesh::read: Scan error while reading node data");
                std::stringstream ss;
                ss << line;
                ss >> tempInts[0 + i1] >> tempInts[chunkSize + i1]
                   >> tempInts[chunkSize * 2 + i1];
                ss >> tempCoords[i1 * numDim];
                if (numDim > 1)
                    ss >> tempCoords[i1 * numDim + 1];
                if (numDim > 2)
                    ss >> tempCoords[i1 * numDim + 2];
                totalNodes++; // When do we quit the infinite loop?
                chunkNodes++; // How many nodes do we actually have in this chunk? It may be smaller than chunkSize.
            }
            if (chunkNodes > chunkSize) {
                throw FinleyException("Mesh::read: error reading chunks of domain, data too large for message size");
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunkSize nodes to each CPU numbered
            // 1 ... mpiInfo->size-1, here goes one of them
            if (nextCPU < mpiInfo->size) {
                // The message has one more int to send chunkNodes
                tempInts[chunkSize * 3] = chunkNodes;
                MPI_Send(tempInts, chunkSize * 3 + 1, MPI_DIM_T, nextCPU, 81720, mpiInfo->comm);
                MPI_Send(tempCoords, chunkSize * numDim, MPI_DOUBLE, nextCPU, 81721, mpiInfo->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of the worker
            // nodes plus one more chunk for the master
            if (nextCPU > mpiInfo->size)
                break; // End infinite loop
        } // Infinite loop
    } // End master
    else { // Worker
#ifdef ESYS_MPI
        // Each worker receives two messages
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize * 3 + 1, MPI_DIM_T, 0, 81720, mpiInfo->comm, &status);
        MPI_Recv(tempCoords, chunkSize * numDim, MPI_DOUBLE, 0, 81721, mpiInfo->comm, &status);
        // How many nodes are in this worker's chunk?
        chunkNodes = tempInts[chunkSize * 3];
#endif
    } // Worker

    // Copy node data from tempMem to domain
    NodeFile* nodes = domain->getNodes();
    nodes->allocTable(chunkNodes);

#pragma omp parallel for
    for (index_t i0 = 0; i0 < chunkNodes; i0++) {
        nodes->Id[i0] = tempInts[0 + i0];
        nodes->globalDegreesOfFreedom[i0] = tempInts[chunkSize + i0];
        nodes->Tag[i0] = tempInts[chunkSize * 2 + i0];
        for (int i1 = 0; i1 < numDim; i1++) {
            nodes->Coordinates[INDEX2(i1, i0, numDim)] = tempCoords[i0 * numDim + i1];
        }
    }
    delete[] tempInts;
    delete[] tempCoords;

    /*************************** read elements ******************************/
    domain->setElements(readElementFile(fileHandle, order, reducedOrder, mpiInfo));

    /************************ read face elements ****************************/
    domain->setFaceElements(readElementFile(fileHandle, order, reducedOrder, mpiInfo));

    /************************ read contact elements ****************************/
    domain->setContactElements(readElementFile(fileHandle, order, reducedOrder, mpiInfo));

    /************************ read nodal elements ***************************/
    domain->setPoints(readElementFile(fileHandle, order, reducedOrder, mpiInfo));

    /************************  get the name tags ****************************/
    std::string remainder;
    size_t len = 0;
    int tag_key;
    if (mpiInfo->rank == 0) { // Master
        // Read the word 'Tags'
        if (!fileHandle.eof()) {
            std::getline(fileHandle, name);
            if (!fileHandle.good())
                throw IOError("Mesh::read: Scan error while reading tag header");
        }
        // Read rest of file in one chunk, after using seek to find length
        std::ios::pos_type cur_pos = fileHandle.tellg();
        fileHandle.seekg(0, std::ios::end);
        std::ios::pos_type end_pos = fileHandle.tellg();
        fileHandle.seekg(cur_pos);
        remainder.resize(end_pos - cur_pos + 1);
        if (!fileHandle.eof()) {
            fileHandle.read(&remainder[0], end_pos-cur_pos);
            if (fileHandle.bad())
                throw IOError("Mesh::read: Error reading remainder");
            remainder[end_pos - cur_pos] = 0;
        }
        len = remainder.find_last_not_of(' ');
        remainder = remainder.substr(0, len+1);
    } // Master

#ifdef ESYS_MPI
    int len_i = static_cast<int>(len);
    MPI_Bcast(&len_i, 1, MPI_INT, 0, mpiInfo->comm);
    len = static_cast<size_t>(len_i);
    if (mpiInfo->rank != 0) {
        remainder.resize(len + 1);
    }
    if (MPI_Bcast(&remainder[0], len+1, MPI_CHAR, 0, mpiInfo->comm) != MPI_SUCCESS)
        throw FinleyException("Mesh::read: broadcast of remainder failed");
#endif

    std::stringstream rem;
    rem << remainder;
    while (std::getline(rem, line)) {
        size_t pos = line.find(' ');
        if (pos != std::string::npos) {
            name = line.substr(0, pos);
            tag_key = std::stoi(line.substr(pos+1));
            domain->setTagMap(name, tag_key);
        }
    }

    // close file
    if (mpiInfo->rank == 0)
        fileHandle.close();

    domain->resolveNodeIds();
    domain->prepare(optimize);
    return domain->getPtr();
}

} // namespace finley

