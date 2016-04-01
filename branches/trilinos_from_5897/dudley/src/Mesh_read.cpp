
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "Mesh.h"

using escript::IOError;

namespace {

using namespace dudley;

ElementFile* readElementFile(FILE* fileHandle, escript::JMPI mpiInfo)
{
    dim_t numEle;
    ElementTypeId typeID = Dudley_NoRef;
    char elementType[1024];
    int scan_ret;

    // Read the element typeID and number of elements
    if (mpiInfo->rank == 0) {
        scan_ret = fscanf(fileHandle, "%s %d\n", elementType, &numEle);
        if (scan_ret == EOF)
            throw IOError("Mesh::read: Scan error while reading file");
        typeID = eltTypeFromString(elementType);
    }
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        dim_t temp1[2];
        temp1[0] = (dim_t)typeID;
        temp1[1] = numEle;
        int mpiError = MPI_Bcast(temp1, 2, MPI_DIM_T, 0, mpiInfo->comm);
        if (mpiError != MPI_SUCCESS) {
            throw DudleyException("Mesh::read: broadcast of element typeID failed");
        }
        typeID = static_cast<ElementTypeId>(temp1[0]);
        numEle = temp1[1];
    }
#endif
    if (typeID == Dudley_NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << elementType;
        throw IOError(ss.str());
    }

    // Allocate the ElementFile
    ElementFile* out = new ElementFile(typeID, mpiInfo);
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
                scan_ret = fscanf(fileHandle, "%d %d",
                                  &tempInts[i0 * (2 + numNodes) + 0],
                                  &tempInts[i0 * (2 + numNodes) + 1]);
                if (scan_ret == EOF)
                    throw IOError("Mesh::read: Scan error while reading file");
                for (int i1 = 0; i1 < numNodes; i1++) {
                    scan_ret = fscanf(fileHandle, " %d",
                                      &tempInts[i0 * (2 + numNodes) + 2 + i1]);
                    if (scan_ret == EOF)
                        throw IOError("Mesh::read: Scan error while reading file");
                }
                scan_ret = fscanf(fileHandle, "\n");
                if (scan_ret == EOF)
                    throw IOError("Mesh::read: Scan error while reading file");
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

    // Copy Element data from tempInts to mesh
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

namespace dudley {

Mesh* Mesh::read(escript::JMPI mpiInfo, const std::string& filename,
                 bool optimize)
{
    dim_t numNodes;
    int numDim = 0;
    char name[1024], frm[20];
    FILE *fileHandle = NULL;
    int scan_ret;

    if (mpiInfo->rank == 0) {
        // open file
        fileHandle = fopen(filename.c_str(), "r");
        if (!fileHandle) {
            std::stringstream ss;
            ss << "Mesh::read: Opening file " << filename
               << " for reading failed.";
            throw DudleyException(ss.str());
        }

        // read header
        sprintf(frm, "%%%d[^\n]", 1023);
        scan_ret = fscanf(fileHandle, frm, name);
        if (scan_ret == EOF)
            throw IOError("Mesh::read: Scan error while reading file");

        // get the number of nodes
        scan_ret = fscanf(fileHandle, "%1d%*s %d\n", &numDim, &numNodes);
        if (scan_ret == EOF)
            throw IOError("Mesh::read: Scan error while reading file");
    }

#ifdef ESYS_MPI
    // MPI Broadcast numDim, numNodes, name if there are multiple MPI procs
    if (mpiInfo->size > 1) {
        dim_t temp1[3];
        if (mpiInfo->rank == 0) {
            temp1[0] = numDim;
            temp1[1] = numNodes;
            temp1[2] = strlen(name) + 1;
        } else {
            temp1[0] = 0;
            temp1[1] = 0;
            temp1[2] = 1;
        }
        MPI_Bcast(temp1, 3, MPI_DIM_T, 0, mpiInfo->comm);
        numDim = temp1[0];
        numNodes = temp1[1];
        MPI_Bcast(name, temp1[2], MPI_CHAR, 0, mpiInfo->comm);
    }
#endif

    // allocate mesh
    Mesh* mesh = new Mesh(name, numDim, mpiInfo);

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
    // chunk into its local mesh.  It doesn't matter that a CPU has the wrong
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
                if (1 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %d %d %le\n",
                                      &tempInts[0 + i1],
                                      &tempInts[chunkSize + i1],
                                      &tempInts[chunkSize * 2 + i1],
                                      &tempCoords[i1 * numDim + 0]);
                } else if (2 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %d %d %le %le\n",
                                      &tempInts[0 + i1],
                                      &tempInts[chunkSize + i1],
                                      &tempInts[chunkSize * 2 + i1],
                                      &tempCoords[i1 * numDim + 0],
                                      &tempCoords[i1 * numDim + 1]);
                } else if (3 == numDim) {
                    scan_ret = fscanf(fileHandle, "%d %d %d %le %le %le\n",
                                      &tempInts[0 + i1],
                                      &tempInts[chunkSize + i1],
                                      &tempInts[chunkSize * 2 + i1],
                                      &tempCoords[i1 * numDim + 0],
                                      &tempCoords[i1 * numDim + 1],
                                      &tempCoords[i1 * numDim + 2]);
                }
                if (scan_ret == EOF)
                    throw IOError("Mesh::read: Scan error while reading file");
                totalNodes++; // When do we quit the infinite loop?
                chunkNodes++; // How many nodes do we actually have in this chunk? It may be smaller than chunkSize.
            }
            if (chunkNodes > chunkSize) {
                throw DudleyException("Mesh::read: error reading chunks of mesh, data too large for message size");
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

    // Copy node data from tempMem to mesh
    mesh->Nodes->allocTable(chunkNodes);

#pragma omp parallel for
    for (index_t i0 = 0; i0 < chunkNodes; i0++) {
        mesh->Nodes->Id[i0] = tempInts[0 + i0];
        mesh->Nodes->globalDegreesOfFreedom[i0] = tempInts[chunkSize + i0];
        mesh->Nodes->Tag[i0] = tempInts[chunkSize * 2 + i0];
        for (int i1 = 0; i1 < numDim; i1++) {
            mesh->Nodes->Coordinates[INDEX2(i1, i0, numDim)] = tempCoords[i0 * numDim + i1];
        }
    }
    delete[] tempInts;
    delete[] tempCoords;

    /*************************** read elements ******************************/
    mesh->Elements = readElementFile(fileHandle, mpiInfo);

    /************************ read face elements ****************************/
    mesh->FaceElements = readElementFile(fileHandle, mpiInfo);

    /************************ read nodal elements ***************************/
    mesh->Points = readElementFile(fileHandle, mpiInfo);

    /************************  get the name tags ****************************/
    char *remainder = NULL, *ptr;
    size_t len = 0;
    int tag_key;
    if (mpiInfo->rank == 0) { // Master
        // Read the word 'Tag'
        if (!feof(fileHandle)) {
            scan_ret = fscanf(fileHandle, "%s\n", name);
            if (scan_ret == EOF)
                throw IOError("Mesh::read: Scan error while reading file");
        }
        // Read rest of file in one chunk, after using seek to find length
        long cur_pos = ftell(fileHandle);
        fseek(fileHandle, 0L, SEEK_END);
        long end_pos = ftell(fileHandle);
        fseek(fileHandle, (long)cur_pos, SEEK_SET);
        remainder = new char[end_pos - cur_pos + 1];
        if (!feof(fileHandle)) {
            scan_ret = fread(remainder, (size_t) end_pos - cur_pos,
                             sizeof(char), fileHandle);
            if (scan_ret == EOF)
                throw IOError("Mesh::read: Scan error while reading file");
            remainder[end_pos - cur_pos] = 0;
        }
        len = strlen(remainder);
        while (len > 1 && isspace(remainder[--len])) {
            remainder[len] = 0;
        }
        len = strlen(remainder);
    } // Master

#ifdef ESYS_MPI
    int len_i = static_cast<int>(len);
    MPI_Bcast(&len_i, 1, MPI_INT, 0, mpiInfo->comm);
    len = static_cast<size_t>(len_i);
    if (mpiInfo->rank != 0) {
        remainder = new char[len + 1];
        remainder[0] = 0;
    }
    if (MPI_Bcast(remainder, len+1, MPI_CHAR, 0, mpiInfo->comm) != MPI_SUCCESS)
        throw DudleyException("Mesh::read: broadcast of remainder failed");
#endif

    if (remainder[0]) {
        ptr = remainder;
        do {
            sscanf(ptr, "%s %d\n", name, &tag_key);
            if (*name)
                mesh->addTagMap(name, tag_key);
            ptr++;
        } while (NULL != (ptr = strchr(ptr, '\n')) && *ptr);
    }
    delete[] remainder;

    // close file
    if (mpiInfo->rank == 0)
        fclose(fileHandle);

    mesh->resolveNodeIds();
    mesh->prepare(optimize);
    return mesh;
}

} // namespace dudley

