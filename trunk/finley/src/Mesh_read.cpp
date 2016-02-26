
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

  Finley: read mesh from file

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include <ctype.h>
#include "Mesh.h"

namespace finley {

Mesh* Mesh::read(esysUtils::JMPI& mpi_info, const std::string fname,
                 int order, int reduced_order, bool optimize)
{
    int numNodes, numDim=0, numEle, i0, i1;
    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    char name[1024], element_type[1024], frm[20];
    FILE *fileHandle_p = NULL;
    ElementTypeId typeID=NoRef;
    int scan_ret;

    if (mpi_info->rank == 0) {
        // get file handle
        fileHandle_p = fopen(fname.c_str(), "r");
        if (fileHandle_p==NULL) {
            std::stringstream ss;
            ss << "Mesh::read: Opening file " << fname << " for reading failed.";
            const std::string msg(ss.str());
            throw escript::IOError(msg);
        }

        // read header
        sprintf(frm,"%%%d[^\n]",1023);
        scan_ret = fscanf(fileHandle_p, frm, name);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");

        // get the number of nodes
        scan_ret = fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");
    }

#ifdef ESYS_MPI
    // Broadcast numDim, numNodes, name if there are multiple MPI procs
    if (mpi_info->size > 1) {
        int temp1[3];
        if (mpi_info->rank == 0) {
            temp1[0] = numDim;
            temp1[1] = numNodes;
            temp1[2] = strlen(name) + 1;
        } else {
            temp1[0] = 0;
            temp1[1] = 0;
            temp1[2] = 1;
        }
        MPI_Bcast(temp1, 3, MPI_INT,  0, mpi_info->comm);
        numDim = temp1[0];
        numNodes = temp1[1];
        MPI_Bcast(name, temp1[2], MPI_CHAR, 0, mpi_info->comm);
    }
#endif

    // allocate mesh
    Mesh* mesh_p = new Mesh(name, numDim, mpi_info);

    // Each CPU will get at most chunkSize nodes so the message has to be
    // sufficiently large
    int chunkSize = numNodes / mpi_info->size + 1;
    int totalNodes = 0;
    int chunkNodes = 0;
    int nextCPU = 1;
    int *tempInts = new int[chunkSize*3+1]; // Stores the integer message data
    double *tempCoords = new double[chunkSize*numDim]; // Stores the double message data

    /*
    Read chunkSize nodes, send it in a chunk to worker CPU which copies chunk into its local mesh_p
    It doesn't matter that a CPU has the wrong nodes for its elements, this is sorted out later
    First chunk sent to CPU 1, second to CPU 2, ...
    Last chunk stays on CPU 0 (the master)
    The three columns of integers (Id, gDOF, Tag) are gathered into a single array tempInts and sent together in a single MPI message
    */

    if (mpi_info->rank == 0) {  /* Master */
        for (;;) {            /* Infinite loop */
#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*3+1; i0++) tempInts[i0] = -1;

#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*numDim; i0++) tempCoords[i0] = -1.0;

            chunkNodes = 0;
            for (i1=0; i1<chunkSize; i1++) {
                if (totalNodes >= numNodes) break;    /* End of inner loop */
                if (1 == numDim) {
                    scan_ret = fscanf(fileHandle_p, "%d %d %d %le\n",
                                        &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
                                        &tempCoords[i1*numDim+0]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                if (2 == numDim) {
                    scan_ret = fscanf(fileHandle_p, "%d %d %d %le %le\n",
                                        &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
                                        &tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                if (3 == numDim) {
                    scan_ret = fscanf(fileHandle_p, "%d %d %d %le %le %le\n",
                                        &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
                                        &tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1], &tempCoords[i1*numDim+2]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                totalNodes++; /* When do we quit the infinite loop? */
                chunkNodes++; /* How many nodes do we actually have in this chunk? It may be smaller than chunkSize. */
            }
            if (chunkNodes > chunkSize) {
                throw FinleyException("Mesh::read: error reading chunks of mesh, data too large for message size");
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunkSize nodes to each CPU numbered 1 ... mpi_info->size-1, here goes one of them
            if (nextCPU < mpi_info->size) {
                tempInts[chunkSize*3] = chunkNodes;   /* The message has one more int to send chunkNodes */
                MPI_Send(tempInts, chunkSize*3+1, MPI_INT, nextCPU, 81720, mpi_info->comm);
                MPI_Send(tempCoords, chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpi_info->comm);
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
            if (nextCPU > mpi_info->size) break; /* End infinite loop */
        } /* Infinite loop */
    }   /* End master */
    else { // Worker
#ifdef ESYS_MPI
        // Each worker receives two messages
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize*3+1, MPI_INT, 0, 81720, mpi_info->comm, &status);
        MPI_Recv(tempCoords, chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpi_info->comm, &status);
        // How many nodes are in this workers chunk?
        chunkNodes = tempInts[chunkSize*3];
#endif
    } // Worker

    /* Copy node data from tempMem to mesh_p */
    mesh_p->Nodes->allocTable(chunkNodes);

#pragma omp parallel for private (i0, i1) schedule(static)
    for (i0=0; i0<chunkNodes; i0++) {
        mesh_p->Nodes->Id[i0] = tempInts[0+i0];
        mesh_p->Nodes->globalDegreesOfFreedom[i0] = tempInts[chunkSize+i0];
        mesh_p->Nodes->Tag[i0] = tempInts[chunkSize*2+i0];
        for (i1=0; i1<numDim; i1++) {
            mesh_p->Nodes->Coordinates[INDEX2(i1,i0,numDim)] = tempCoords[i0*numDim+i1];
        }
    }
    delete[] tempInts;
    delete[] tempCoords;

    /* ************************  read elements ******************************/
    /* Read the element typeID */
    if (mpi_info->rank == 0) {
        scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");
        typeID=ReferenceElement::getTypeId(element_type);
    }
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        int temp1[2], mpi_error;
        temp1[0] = (int) typeID;
        temp1[1] = numEle;
        mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
        if (mpi_error != MPI_SUCCESS) {
            throw FinleyException("Mesh::read: broadcast of Element typeID failed");
        }
        typeID = (ElementTypeId) temp1[0];
        numEle = temp1[1];
    }
#endif
    if (typeID==NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << element_type;
        const std::string msg(ss.str());
        throw escript::IOError(msg);
    }

    // Allocate the ElementFile
    refElements.reset(new ReferenceElementSet(typeID, order, reduced_order));
    mesh_p->Elements=new ElementFile(refElements, mpi_info);
    // new meaning for numNodes: num nodes per element
    numNodes = mesh_p->Elements->numNodes;

    /********************** Read the element data ***************************/
    chunkSize = numEle / mpi_info->size + 1;
    int totalEle=0;
    int chunkEle=0;
    nextCPU=1;
    tempInts = new int[chunkSize*(2+numNodes)+1]; /* Store Id + Tag + node list (+ one int at end for chunkEle) */
    /* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
    if (mpi_info->rank == 0) {  /* Master */
        for (;;) {            /* Infinite loop */
#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;

            chunkEle = 0;
            for (i0=0; i0<chunkSize; i0++) {
                if (totalEle >= numEle) break; /* End inner loop */
                scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                for (i1 = 0; i1 < numNodes; i1++) {
                    scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                scan_ret = fscanf(fileHandle_p, "\n");
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunk of elements to each CPU except 0
            // itself, here goes one of them
            if (nextCPU < mpi_info->size) {
                tempInts[chunkSize*(2+numNodes)] = chunkEle;
                MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81722, mpi_info->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of the
            // worker nodes plus one more chunk for the master
            if (nextCPU > mpi_info->size)
                break; // End infinite loop
        } // Infinite loop
    } // End master
    else { // Worker
#ifdef ESYS_MPI
        // Each worker receives one message
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81722, mpi_info->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
    } // Worker

    mesh_p->Elements->allocTable(chunkEle);

    // Copy Element data from tempInts to mesh_p
    mesh_p->Elements->minColor=0;
    mesh_p->Elements->maxColor=chunkEle-1;
#pragma omp parallel for private (i0, i1) schedule(static)
    for (i0=0; i0<chunkEle; i0++) {
        mesh_p->Elements->Id[i0] = tempInts[i0*(2+numNodes)+0];
        mesh_p->Elements->Tag[i0] = tempInts[i0*(2+numNodes)+1];
        mesh_p->Elements->Owner[i0] = mpi_info->rank;
        mesh_p->Elements->Color[i0] = i0;
        for (i1 = 0; i1 < numNodes; i1++) {
            mesh_p->Elements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    delete[] tempInts;

    /******************** end of Read the element data **********************/

    /********************** read face elements ******************************/

    // Read the element typeID
    if (mpi_info->rank == 0) {
        scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");
        typeID=ReferenceElement::getTypeId(element_type);
    }
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        int temp1[2];
        temp1[0] = (int) typeID;
        temp1[1] = numEle;
        MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
        typeID = (ElementTypeId) temp1[0];
        numEle = temp1[1];
    }
#endif
    if (typeID==NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << element_type;
        const std::string msg(ss.str());
        throw escript::IOError(msg);
    }
    // Allocate the ElementFile
    refFaceElements.reset(new ReferenceElementSet(typeID, order, reduced_order));
    mesh_p->FaceElements=new ElementFile(refFaceElements, mpi_info);
    numNodes = mesh_p->FaceElements->numNodes; // new meaning for numNodes: num nodes per element

    /*********************** Read the face element data *********************/
    chunkSize = numEle / mpi_info->size + 1;
    totalEle=0;
    nextCPU=1;
    chunkEle=0;
    tempInts = new int[chunkSize*(2+numNodes)+1];
    // Store Id + Tag + node list (+ one int at end for chunkEle)
    // Elements are specified as a list of integers...only need one message
    // instead of two as with the nodes
    if (mpi_info->rank == 0) { // Master
        for (;;) { // Infinite loop
#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;

            chunkEle = 0;
            for (i0=0; i0<chunkSize; i0++) {
                if (totalEle >= numEle) break; /* End inner loop */
                scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                for (i1 = 0; i1 < numNodes; i1++) {
                    scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                scan_ret = fscanf(fileHandle_p, "\n");
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            /* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
            if (nextCPU < mpi_info->size) {
                tempInts[chunkSize*(2+numNodes)] = chunkEle;
                MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81723, mpi_info->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of the
            // worker nodes plus one more chunk for the master
            if (nextCPU > mpi_info->size)
                break; // End infinite loop
        } // Infinite loop
    } // End master
    else {  /* Worker */
#ifdef ESYS_MPI
        /* Each worker receives one message */
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81723, mpi_info->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
    } // Worker

    mesh_p->FaceElements->allocTable(chunkEle);

    // Copy Element data from tempInts to mesh_p
    mesh_p->FaceElements->minColor=0;
    mesh_p->FaceElements->maxColor=chunkEle-1;
#pragma omp parallel for private (i0, i1)
    for (i0=0; i0<chunkEle; i0++) {
        mesh_p->FaceElements->Id[i0] = tempInts[i0*(2+numNodes)+0];
        mesh_p->FaceElements->Tag[i0] = tempInts[i0*(2+numNodes)+1];
        mesh_p->FaceElements->Owner[i0] = mpi_info->rank;
        mesh_p->FaceElements->Color[i0] = i0;
        for (i1 = 0; i1 < numNodes; i1++) {
            mesh_p->FaceElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    delete[] tempInts;

    /******************* end of Read the face element data ******************/

    /************************* read contact elements ************************/

    // Read the element typeID
    if (mpi_info->rank == 0) {
        scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");
        typeID=ReferenceElement::getTypeId(element_type);
    }
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        int temp1[2];
        temp1[0] = (int) typeID;
        temp1[1] = numEle;
        MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
        typeID = (ElementTypeId) temp1[0];
        numEle = temp1[1];
    }
#endif
    if (typeID==NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << element_type;
        const std::string msg(ss.str());
        throw escript::IOError(msg);
    }

    // Allocate the ElementFile
    refContactElements.reset(new ReferenceElementSet(typeID, order, reduced_order));
    mesh_p->ContactElements=new ElementFile(refContactElements, mpi_info);
    numNodes = mesh_p->ContactElements->numNodes; // new meaning for numNodes: num nodes per element

    /******************* Read the contact element data **********************/
    chunkSize = numEle / mpi_info->size + 1;
    totalEle=0;
    nextCPU=1;
    chunkEle=0;
    tempInts = new int[chunkSize*(2+numNodes)+1];
    // Store Id + Tag + node list (+ one int at end for chunkEle)
    // Elements are specified as a list of integers...only need one message instead of two as with the nodes
    if (mpi_info->rank == 0) { // Master
        for (;;) { // Infinite loop
#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;

            chunkEle = 0;
            for (i0=0; i0<chunkSize; i0++) {
                if (totalEle >= numEle) break; /* End inner loop */
                scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                for (i1 = 0; i1 < numNodes; i1++) {
                    scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                scan_ret = fscanf(fileHandle_p, "\n");
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunk of elements to each CPU except
            // 0 itself, here goes one of them
            if (nextCPU < mpi_info->size) {
                tempInts[chunkSize*(2+numNodes)] = chunkEle;
                MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81724, mpi_info->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of the
            // worker nodes plus one more chunk for the master
            if (nextCPU > mpi_info->size)
                break; // End infinite loop
        } // Infinite loop
    } // End master
    else { // Worker
#ifdef ESYS_MPI
        // Each worker receives one message
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81724, mpi_info->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)] ;
#endif
    } // Worker

    // Copy Element data from tempInts to mesh_p
    mesh_p->ContactElements->allocTable(chunkEle);
    mesh_p->ContactElements->minColor=0;
    mesh_p->ContactElements->maxColor=chunkEle-1;
#pragma omp parallel for private (i0, i1)
    for (i0=0; i0<chunkEle; i0++) {
        mesh_p->ContactElements->Id[i0] = tempInts[i0*(2+numNodes)+0];
        mesh_p->ContactElements->Tag[i0] = tempInts[i0*(2+numNodes)+1];
        mesh_p->ContactElements->Owner[i0] = mpi_info->rank;
        mesh_p->ContactElements->Color[i0] = i0;
        for (i1 = 0; i1 < numNodes; i1++) {
            mesh_p->ContactElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    delete[] tempInts;

    // ****************** read nodal elements ******************

    // ***************  Read the element typeID ***********
    if (mpi_info->rank == 0) {
        scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
        if (scan_ret==EOF)
            throw escript::IOError("Mesh::read: scan error while reading file");
        typeID=ReferenceElement::getTypeId(element_type);
    }
#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        int temp1[2];
        temp1[0] = (int) typeID;
        temp1[1] = numEle;
        MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
        typeID = (ElementTypeId) temp1[0];
        numEle = temp1[1];
    }
#endif
    if (typeID==NoRef) {
        std::stringstream ss;
        ss << "Mesh::read: Unidentified element type " << element_type;
        const std::string msg(ss.str());
        throw escript::IOError(msg);
    }

    // Allocate the ElementFile
    refPoints.reset(new ReferenceElementSet(typeID, order, reduced_order));
    mesh_p->Points=new ElementFile(refPoints, mpi_info);
    // New meaning for numNodes: num nodes per element
    numNodes = mesh_p->Points->numNodes;

    // ******************* Read the nodal element data ****************
    chunkSize = numEle / mpi_info->size + 1;
    totalEle=0;
    nextCPU=1;
    chunkEle=0;
    // Store Id + Tag + node list (+ one int at end for chunkEle)
    tempInts = new int[chunkSize*(2+numNodes)+1];
    // Elements are specified as a list of integers...only need one
    // message instead of two as with the nodes
    if (mpi_info->rank == 0) {  // Master
        for (;;) { // Infinite loop
#pragma omp parallel for private (i0) schedule(static)
            for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;

            chunkEle = 0;
            for (i0=0; i0<chunkSize; i0++) {
                if (totalEle >= numEle) break; /* End inner loop */
                scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                for (i1 = 0; i1 < numNodes; i1++) {
                    scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (scan_ret==EOF)
                        throw escript::IOError("Mesh::read: scan error while reading file");
                }
                scan_ret = fscanf(fileHandle_p, "\n");
                if (scan_ret==EOF)
                    throw escript::IOError("Mesh::read: scan error while reading file");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            // Eventually we'll send chunk of elements to each CPU
            // except 0 itself, here goes one of them
            if (nextCPU < mpi_info->size) {
                tempInts[chunkSize*(2+numNodes)] = chunkEle;
                MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81725, mpi_info->comm);
            }
#endif
            nextCPU++;
            // Infinite loop ends when I've read a chunk for each of
            // the worker nodes plus one more chunk for the master
            if (nextCPU > mpi_info->size)
                break; // End infinite loop
        } // Infinite loop
    }   // End master
    else {  // Worker
#ifdef ESYS_MPI
        // Each worker receives one message
        MPI_Status status;
        MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81725, mpi_info->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
    } // Worker

    // Copy Element data from tempInts to mesh_p
    mesh_p->Points->allocTable(chunkEle);
    mesh_p->Points->minColor=0;
    mesh_p->Points->maxColor=chunkEle-1;
#pragma omp parallel for private (i0, i1) schedule(static)
    for (i0=0; i0<chunkEle; i0++) {
        mesh_p->Points->Id[i0] = tempInts[i0*(2+numNodes)+0];
        mesh_p->Points->Tag[i0] = tempInts[i0*(2+numNodes)+1];
        mesh_p->Points->Owner[i0] = mpi_info->rank;
        mesh_p->Points->Color[i0] = i0;
        for (i1 = 0; i1 < numNodes; i1++) {
            mesh_p->Points->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    delete[] tempInts;

    // ************** end of Read the nodal element data ***************

    // *****************  get the name tags ********************************
    char *remainder=NULL, *ptr;
    size_t len=0;
    int tag_key;
    if (mpi_info->rank == 0) {  // Master
        // Read the word 'Tag'
        if (!feof(fileHandle_p)) {
            scan_ret = fscanf(fileHandle_p, "%s\n", name);
            if (scan_ret==EOF)
                throw escript::IOError("Mesh::read: scan error while reading file");
        }

        // Read rest of file in one chunk, after using seek to find length
        long cur_pos = ftell(fileHandle_p);
        fseek(fileHandle_p, 0L, SEEK_END);
        long end_pos = ftell(fileHandle_p);
        fseek(fileHandle_p, (long)cur_pos, SEEK_SET);
        remainder = new char[end_pos-cur_pos+1];
        if (!feof(fileHandle_p)) {
            scan_ret = fread(remainder, (size_t) end_pos-cur_pos,
                             sizeof(char), fileHandle_p);

            if (scan_ret==EOF)
                throw escript::IOError("Mesh::read: scan error while reading file");
            remainder[end_pos-cur_pos] = 0;
        }
        len = strlen(remainder);
        // trim the string
        while (len>1 && isspace(remainder[--len])) {
            remainder[len]=0;
        }
        len = strlen(remainder);
    } // Master
#ifdef ESYS_MPI
    int len_i=static_cast<int>(len);
    MPI_Bcast(&len_i, 1, MPI_INT,  0, mpi_info->comm);
    len = static_cast<size_t>(len_i);
    if (mpi_info->rank != 0) {
        remainder = new char[len+1];
        remainder[0] = 0;
    }
    if (MPI_Bcast (remainder, len+1, MPI_CHAR,  0, mpi_info->comm) !=
            MPI_SUCCESS)
        throw FinleyException("Mesh::read: broadcast of remainder failed");
#endif

    if (remainder[0]) {
        ptr = remainder;
        do {
            sscanf(ptr, "%s %d\n", name, &tag_key);
            if (*name)
                mesh_p->addTagMap(name, tag_key);
            ptr++;
        } while(NULL != (ptr = strchr(ptr, '\n')) && *ptr);
    }
    delete[] remainder;

    // close file
    if (mpi_info->rank == 0)
        fclose(fileHandle_p);

    // resolve id's and rearrange elements
    mesh_p->resolveNodeIds();
    mesh_p->prepare(optimize);
    return mesh_p;
}

} // namespace finley

