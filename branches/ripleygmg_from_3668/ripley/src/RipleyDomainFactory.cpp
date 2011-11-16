
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

#include <ripley/RipleyDomainFactory.h>
#include <ripley/RipleyError.h>
extern "C" {
#include <esysUtils/blocktimer.h>
}

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/python/extract.hpp>
#include <sstream>

using namespace std;

namespace ripley {

#ifdef USE_NETCDF
// A convenience method to retrieve an integer attribute from a netCDF file
int ncGetIntAttribute(NcFile *dataFile, const string &fn, const string &an)
{
    NcAtt *attr;
    stringstream msg;
    if (! (attr=dataFile->get_att(an.c_str())) ) {
        msg << "loadMesh: Error retrieving integer attribute '"
                << an << "' from netCDF file '" << fn << "'";
        throw RipleyException(msg.str());
    }
    int val = attr->as_int(0);
    delete attr;
    return val;
}
#endif

escript::Domain_ptr loadMesh(const string& fileName)
{
#ifdef USE_NETCDF
    Esys_MPIInfo *mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    checkPasoError();
    char *fn = Esys_MPI_appendRankToFileName(fileName.c_str(), mpiInfo->size,
                                             mpiInfo->rank);
    string fName(fn);
    TMPMEMFREE(fn);
    double blocktimer_start = blocktimer_time();

    // set netCDF error handler
    NcError err(NcError::silent_nonfatal);
    NcFile dataFile(fName.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
        Esys_MPIInfo_free(mpiInfo);
        stringstream msg;
        msg << "loadMesh: Opening netCDF file '" << fName
            << "' for reading failed.";
        throw RipleyException(msg.str());
    }

    NcAtt *attr;
    NcVar *nc_var_temp;
    int mpiSize, mpiRank, numDim, numNodes, numElements, numFaceElements,
        numPoints, elementId, faceElementId, pointsId, numTags;

    // read integer attributes
    try {
        mpiSize = ncGetIntAttribute(&dataFile, fName, "mpi_size");
        mpiRank = ncGetIntAttribute(&dataFile, fName, "mpi_rank");
        numDim = ncGetIntAttribute(&dataFile, fName, "numDim");
        numNodes = ncGetIntAttribute(&dataFile, fName, "numNodes");
        numElements = ncGetIntAttribute(&dataFile, fName, "num_Elements");
        numFaceElements = ncGetIntAttribute(&dataFile, fName, "num_FaceElements");
        numPoints = ncGetIntAttribute(&dataFile, fName, "num_Points");
        elementId = ncGetIntAttribute(&dataFile, fName, "Elements_TypeId");
        faceElementId = ncGetIntAttribute(&dataFile, fName, "FaceElements_TypeId");
        pointsId = ncGetIntAttribute(&dataFile, fName, "Points_TypeId");
        numTags = ncGetIntAttribute(&dataFile, fName, "num_Tags");
    } catch (RipleyException(e)) {
        Esys_MPIInfo_free(mpiInfo);
        throw e;
    }

    // Verify size and rank
    if (mpiInfo->size != mpiSize) {
        Esys_MPIInfo_free(mpiInfo);
        stringstream msg;
        msg << "loadMesh: The netCDF file '" << fName
            << "' can only be read on " << mpiSize << " CPUs instead of "
            << mpiInfo->size;
        throw RipleyException(msg.str());
    }
    if (mpiInfo->rank != mpiRank) {
        Esys_MPIInfo_free(mpiInfo);
        stringstream msg;
        msg << "loadMesh: The netCDF file '" << fName
            << "' should be read on CPU #" << mpiRank << " instead of "
            << mpiInfo->rank;
        throw RipleyException(msg.str());
    }

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
        Esys_MPIInfo_free(mpiInfo);
        stringstream msg;
        msg << "loadMesh: Error retrieving mesh name from netCDF file '"
            << fName << "'";
        throw RipleyException(msg.str());
    }
    char *name = attr->as_string(0);
    delete attr;

    /* allocate mesh */
    RipleyDomain *dom=new RipleyDomain(name, numDim, mpiInfo);

    try {
        const string msgPrefix("loadMesh: netCDF operation failed - ");

        // read nodes
        NodeFile_ptr nodes(new NodeFile(numDim, mpiInfo));
        nodes->readFromNetCDF(dataFile, numNodes);
        dom->setNodes(nodes);

        // read elements
        ElementFile_ptr elements(new ElementFile(
                static_cast<ElementTypeId>(elementId), mpiInfo));
        elements->readFromNetCDF(dataFile, numElements, "Elements");
        dom->setElements(elements);

        // read face elements
        elements.reset(new ElementFile(
                static_cast<ElementTypeId>(faceElementId), mpiInfo));
        elements->readFromNetCDF(dataFile, numFaceElements, "FaceElements");
        dom->setFaceElements(elements);

        // read point elements
        elements.reset(new ElementFile(
                static_cast<ElementTypeId>(pointsId), mpiInfo));
        elements->readFromNetCDF(dataFile, numPoints, "Points");
        dom->setPoints(elements);

        // read tags
        if (numTags > 0) {
            // Temp storage to gather node IDs
            vector<int> tagKeys(numTags);
            char name_temp[4096];
            int i;

            // tagKeys
            if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) )
                throw RipleyException(msgPrefix+"get_var(Tags_keys)");

            if (! nc_var_temp->get(&tagKeys[0], numTags) )
                throw RipleyException(msgPrefix+"get(Tags_keys)");

            for (i=0; i<numTags; i++) {
                // Retrieve tag name
                sprintf(name_temp, "Tags_name_%d", i);
                if (! (attr=dataFile.get_att(name_temp)) ) {
                    stringstream msg;
                    msg << "get_att(" << name_temp << ")";
                    throw RipleyException(msgPrefix+msg.str());
                }
                char *name = attr->as_string(0);
                delete attr;
                dom->setTagMap(name, tagKeys[i]);
            }
        }

        // Nodes_DofDistribution
        IndexVector first_DofComponent(mpiSize+1);
        if (! (nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) )
            throw RipleyException(msgPrefix+"get_var(Nodes_DofDistribution)");
        if (!nc_var_temp->get(&first_DofComponent[0], mpiSize+1))
            throw RipleyException(msgPrefix+"get(Nodes_DofDistribution)");

        // Nodes_NodeDistribution
        IndexVector first_NodeComponent(mpiSize+1);
        if (! (nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) )
            throw RipleyException(msgPrefix+"get_var(Nodes_NodeDistribution)");
        if (!nc_var_temp->get(&first_NodeComponent[0], mpiSize+1))
            throw RipleyException(msgPrefix+"get(Nodes_NodeDistribution)");
        dom->createMappings(first_DofComponent, first_NodeComponent);

    } catch (RipleyException(e)) {
        Esys_MPIInfo_free(mpiInfo);
        throw e;
    }

    blocktimer_increment("LoadMesh()", blocktimer_start);
    return dom->getPtr();
#else
    throw RipleyException("loadMesh: not compiled with netCDF. Please contact your installation manager.");
#endif /* USE_NETCDF */
}

escript::Domain_ptr readMesh(const string& filename, bool optimize)
{
    if (filename.size() == 0)
        throw RipleyException("readMesh: Empty file name!");

    double blocktimer_start = blocktimer_time();
    dim_t i0, numNodes, numDim, numEle;
    char name[LenString_MAX], elementType[LenString_MAX], frm[20];
    FILE *handle = NULL;

    Esys_MPIInfo *mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    checkPasoError();

    if (mpiInfo->rank == 0) {
        /* get file handle */
        handle = fopen(filename.c_str(), "r");
        if (handle == NULL) {
            Esys_MPIInfo_free(mpiInfo);
            stringstream msg;
            msg << "readMesh: Opening file " << filename << " for reading failed.";
            throw RipleyException(msg.str());
        }

        /* read header */
        sprintf(frm, "%%%d[^\n]", LenString_MAX - 1);
        int ret = fscanf(handle, frm, name);
        if (ret!=1) {
            Esys_MPIInfo_free(mpiInfo);
            throw RipleyException("readMesh: Invalid file!");
        }
        /* get the number of nodes */
        ret = fscanf(handle, "%1d%*s %d\n", &numDim, &numNodes);
        if (ret!=2) {
            Esys_MPIInfo_free(mpiInfo);
            throw RipleyException("readMesh: Invalid file!");
        }
    }
#ifdef ESYS_MPI
    /* MPI Broadcast numDim, numNodes, name if there are multiple MPI procs */
    if (mpiInfo->size > 1) {
        int values[3];
        if (mpiInfo->rank == 0) {
            values[0] = numDim;
            values[1] = numNodes;
            values[2] = strlen(name) + 1;
        } else {
            values[0] = 0;
            values[1] = 0;
            values[2] = 1;
        }
        MPI_Bcast(values, 3, MPI_INT, 0, mpiInfo->comm);
        numDim = values[0];
        numNodes = values[1];
        MPI_Bcast(name, values[2], MPI_CHAR, 0, mpiInfo->comm);
    }
#endif

    /* allocate mesh */
    RipleyDomain *dom=new RipleyDomain(name, numDim, mpiInfo);

    /* Each CPU will get at most chunkSize nodes so the message has to
     * be sufficiently large */
    int chunkSize = numNodes / mpiInfo->size + 1, totalNodes = 0, chunkNodes = 0, nextCPU = 1;
    /* stores the integer message data */
    vector<int> tempInts(chunkSize*3+1, -1);
    /* stores the double message data */
    vector<double> tempCoords(chunkSize*numDim, -1.0);

    /*
       Read chunkSize nodes, send it in a chunk to worker CPU which copies
       chunk into its local mesh. It doesn't matter that a CPU has the wrong
       nodes for its elements, this is sorted out later. First chunk sent to
       CPU 1, second to CPU 2, ... Last chunk stays on CPU 0 (the master).
       The three columns of integers (Id, gDOF, Tag) are gathered into a
       single array tempInts and sent together in a single MPI message.
     */

    if (mpiInfo->rank == 0) { /* Master */
        for (;;) { /* Infinite loop */
            chunkNodes = 0;
            for (dim_t i = 0; i < chunkSize; i++) {
                if (totalNodes >= numNodes)
                    break;  /* End of inner loop */
                if (1 == numDim) {
                    int ret = fscanf(handle, "%d %d %d %le\n",
                            &tempInts[i], &tempInts[chunkSize+i],
                            &tempInts[chunkSize*2+i], &tempCoords[i*numDim]);
                    if (ret!=4)
                        throw RipleyException("readMesh: Invalid file!");
                }
                if (2 == numDim) {
                    int ret = fscanf(handle, "%d %d %d %le %le\n",
                            &tempInts[i], &tempInts[chunkSize+i],
                            &tempInts[chunkSize*2+i], &tempCoords[i*numDim],
                            &tempCoords[i*numDim+1]);
                    if (ret!=5)
                        throw RipleyException("readMesh: Invalid file!");
                }
                if (3 == numDim) {
                    int ret = fscanf(handle, "%d %d %d %le %le %le\n",
                            &tempInts[i], &tempInts[chunkSize+i],
                            &tempInts[chunkSize*2+i], &tempCoords[i*numDim],
                            &tempCoords[i*numDim+1], &tempCoords[i*numDim+2]);
                    if (ret!=6)
                        throw RipleyException("readMesh: Invalid file!");
                }
                totalNodes++; /* When do we quit the infinite loop? */
                /* Count number of nodes we actually have - it may be less
                 * than chunkSize. */
                chunkNodes++;
            }
            if (chunkNodes > chunkSize)
                throw RipleyException("readMesh: error reading chunks of mesh, data too large for message size");

#ifdef ESYS_MPI
            /* Eventually we'll send chunkSize nodes to each CPU numbered
             * 1 ... mpiInfo->size-1, here goes one of them */
            if (nextCPU < mpiInfo->size) {
                // the message has one more int to send chunkNodes
                tempInts[chunkSize * 3] = chunkNodes;
                MPI_Send(&tempInts[0], chunkSize*3+1, MPI_INT, nextCPU, 81720, mpiInfo->comm);
                MPI_Send(&tempCoords[0], chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpiInfo->comm);
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the
             * worker nodes plus one more chunk for the master */
            if (nextCPU > mpiInfo->size)
                break; // End infinite loop
        } // infinite loop

    } else { /* Worker */
#ifdef ESYS_MPI
        // Each worker receives two messages
        MPI_Status status;
        MPI_Recv(&tempInts[0], chunkSize*3+1, MPI_INT, 0, 81720, mpiInfo->comm, &status);
        MPI_Recv(&tempCoords[0], chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpiInfo->comm, &status);
        // how many nodes are in this workers chunk?
        chunkNodes = tempInts[chunkSize*3];
#endif
    } /* Worker */

    /* Copy data from temporary arrays to node file */
    NodeFile_ptr nodes(new NodeFile(numDim, mpiInfo));
    IndexVector id(&tempInts[0], &tempInts[chunkNodes]);
    IndexVector gDOF(&tempInts[chunkSize], &tempInts[chunkSize+chunkNodes]);
    IndexVector tag(&tempInts[chunkSize*2], &tempInts[chunkSize*2+chunkNodes]);
    vector<double> coordinates(&tempCoords[0], &tempCoords[chunkNodes*numDim]);
    nodes->swapEntries(id, tag, gDOF, coordinates);
    dom->setNodes(nodes);
    tempInts.clear();
    tempCoords.clear();

    /*************************** read elements ******************************/

    // read the element typeID
    ElementTypeId typeID = InvalidElementType;
    if (mpiInfo->rank == 0) {
        int ret = fscanf(handle, "%s %d\n", elementType, &numEle);
        if (ret!=2)
            throw RipleyException("readMesh: Invalid file!");
        typeID = ElementFile::ripleyTypeFromString(elementType);
    }
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        int array[2];
        array[0] = static_cast<int>(typeID);
        array[1] = numEle;
        if (MPI_Bcast(array, 2, MPI_INT, 0, mpiInfo->comm) != MPI_SUCCESS) {
            string msg("readMesh: Broadcast of element type id failed");
            throw RipleyException(msg);
        }
        typeID = static_cast<ElementTypeId>(array[0]);
        numEle = array[1];
    }
#endif
    if (typeID == InvalidElementType) {
        stringstream msg;
        msg << "readMesh: Unidentified element type " << elementType;
        throw RipleyException(msg.str());
    }

    // allocate the ElementFile
    ElementFile_ptr elements(new ElementFile(typeID, mpiInfo));
    // new meaning for numNodes: number of nodes per element
    numNodes = elements->getNumNodes();

    /*********************** Read the element data **************************/
    chunkSize = numEle / mpiInfo->size + 1;
    int totalEle = 0, chunkEle = 0;
    nextCPU = 1;
    /* Store Id + Tag + node list (+ one int at end for chunkEle) */
    tempInts.resize(chunkSize*(2+numNodes)+1, -1);
    /* Elements are specified as a list of integers.
     * Only need one message instead of two as with the nodes */
    if (mpiInfo->rank == 0) { /* Master */
        for (;;) { /* Infinite loop */
            chunkEle = 0;
            for (i0 = 0; i0 < chunkSize; i0++) {
                if (totalEle >= numEle)
                    break; /* End inner loop */
                int ret = fscanf(handle, "%d %d", &tempInts[i0*(2+numNodes)],
                        &tempInts[i0*(2+numNodes)+1]);
                if (ret!=2)
                    throw RipleyException("readMesh: Invalid file!");
                for (dim_t i1 = 0; i1 < numNodes; i1++) {
                    ret = fscanf(handle, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (ret!=1)
                        throw RipleyException("readMesh: Invalid file!");
                }
                ret = fscanf(handle, "\n");
                if (ret!=0)
                    throw RipleyException("readMesh: Invalid file!");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            /* Eventually we'll send chunk of elements to each CPU except
             * 0 itself, here goes one of them */
            if (nextCPU < mpiInfo->size) {
                tempInts[chunkSize*(2+numNodes)] = chunkEle;
                MPI_Send(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81722, mpiInfo->comm);
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the
             * worker nodes plus one more chunk for the master */
            if (nextCPU > mpiInfo->size)
                break; /* End infinite loop */
        } // Infinite loop

    } else { // Worker
#ifdef ESYS_MPI
        // each worker receives one message
        MPI_Status status;
        MPI_Recv(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, 0, 81722, mpiInfo->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
    } // Worker

    // copy element data from temporary array to mesh
    id.resize(chunkEle);
    tag.resize(chunkEle);
    IndexVector color(chunkEle, -1);
    RankVector owner(chunkEle);
    IndexVector nodesVec(chunkEle*numNodes);
#pragma omp parallel for private (i0) schedule(static)
    for (i0 = 0; i0 < chunkEle; i0++) {
        id[i0] = tempInts[i0*(2+numNodes)];
        tag[i0] = tempInts[i0*(2+numNodes)+1];
        owner[i0] = mpiInfo->rank;
        for (dim_t i1 = 0; i1 < numNodes; i1++) {
            nodesVec[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    elements->swapEntries(id, tag, color, owner, nodesVec);
    tempInts.clear();
    dom->setElements(elements);

    /************************* read face elements ***************************/

    // read the element typeID
    if (mpiInfo->rank == 0) {
        int ret = fscanf(handle, "%s %d\n", elementType, &numEle);
        if (ret!=2)
            throw RipleyException("readMesh: Invalid file!");
        typeID = ElementFile::ripleyTypeFromString(elementType);
    }
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        int array[2];
        array[0] = (int)typeID;
        array[1] = numEle;
        if (MPI_Bcast(array, 2, MPI_INT, 0, mpiInfo->comm) != MPI_SUCCESS) {
            string msg("readMesh: Broadcast of element type id failed");
            throw RipleyException(msg);
        }
        typeID = static_cast<ElementTypeId>(array[0]);
        numEle = array[1];
    }
#endif
    if (typeID == InvalidElementType) {
        stringstream msg;
        msg << "readMesh: Unidentified element type " << elementType;
        throw RipleyException(msg.str());
    }

    // allocate the ElementFile
    ElementFile_ptr faceElements(new ElementFile(typeID, mpiInfo));

    // new meaning for numNodes: number of nodes per face element */
    numNodes = faceElements->getNumNodes();

    /********************* Read the face element data ***********************/
    chunkSize = numEle / mpiInfo->size + 1;
    totalEle = 0;
    chunkEle = 0;
    nextCPU = 1;
    /* Store Id + Tag + node list (+ one int at end for chunkEle) */
    tempInts.resize(chunkSize*(2+numNodes)+1, -1);
    /* Elements are specified as a list of integers.
     * Only need one message instead of two as with the nodes */
    if (mpiInfo->rank == 0) { /* Master */
        for (;;) { /* Infinite loop */
            chunkEle = 0;
            for (i0 = 0; i0 < chunkSize; i0++) {
                if (totalEle >= numEle)
                    break;  /* End inner loop */
                int ret = fscanf(handle, "%d %d", &tempInts[i0*(2+numNodes)],
                        &tempInts[i0*(2+numNodes)+1]);
                if (ret!=2)
                    throw RipleyException("readMesh: Invalid file!");
                for (dim_t i1 = 0; i1 < numNodes; i1++) {
                    ret = fscanf(handle, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (ret!=1)
                        throw RipleyException("readMesh: Invalid file!");
                }
                ret = fscanf(handle, "\n");
                if (ret!=0)
                    throw RipleyException("readMesh: Invalid file!");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            /* Eventually we'll send chunk of elements to each CPU except
             * 0 itself, here goes one of them */
            if (nextCPU < mpiInfo->size) {
                tempInts[chunkSize * (2 + numNodes)] = chunkEle;
                MPI_Send(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81723, mpiInfo->comm);
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the
             * worker nodes plus one more chunk for the master */
            if (nextCPU > mpiInfo->size)
                break; /* End infinite loop */
        } // Infinite loop
    } else { // Worker
#ifdef ESYS_MPI
        // each worker receives one message
        MPI_Status status;
        MPI_Recv(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, 0, 81723, mpiInfo->comm, &status);
        chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
    } // Worker

    // copy data from temporary array to mesh
    id.resize(chunkEle);
    tag.resize(chunkEle);
    color.assign(chunkEle, -1);
    owner.resize(chunkEle);
    nodesVec.resize(chunkEle*numNodes);
#pragma omp parallel for private (i0)
    for (i0 = 0; i0 < chunkEle; i0++) {
        id[i0] = tempInts[i0*(2+numNodes)];
        tag[i0] = tempInts[i0*(2+numNodes)+1];
        owner[i0] = mpiInfo->rank;
        for (dim_t i1 = 0; i1 < numNodes; i1++) {
            nodesVec[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    faceElements->swapEntries(id, tag, color, owner, nodesVec);
    tempInts.clear();
    dom->setFaceElements(faceElements);

    /************************* read nodal elements **************************/

    // read the element typeID
    if (mpiInfo->rank == 0) {
        int ret = fscanf(handle, "%s %d\n", elementType, &numEle);
        if (ret!=2)
            throw RipleyException("readMesh: Invalid file!");
        typeID = ElementFile::ripleyTypeFromString(elementType);
    }
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        int array[2];
        array[0] = (int)typeID;
        array[1] = numEle;
        if (MPI_Bcast(array, 2, MPI_INT, 0, mpiInfo->comm) != MPI_SUCCESS) {
            string msg("readMesh: Broadcast of element type id failed");
            throw RipleyException(msg);
        }
        typeID = static_cast<ElementTypeId>(array[0]);
        numEle = array[1];
    }
#endif
    if (typeID == InvalidElementType) {
        stringstream msg;
        msg << "readMesh: Unidentified element type " << elementType;
        throw RipleyException(msg.str());
    }

    // allocate the ElementFile
    ElementFile_ptr points(new ElementFile(typeID, mpiInfo));
    // new meaning for numNodes: number of nodes per point element
    numNodes = points->getNumNodes();

    /********************** Read the nodal element data *********************/
    chunkSize = numEle / mpiInfo->size + 1;
    totalEle = 0;
    nextCPU = 1;
    chunkEle = 0;
    /* Store Id + Tag + node list (+ one int at end for chunkEle) */
    tempInts.resize(chunkSize*(2+numNodes)+1, -1);
    /* Elements are specified as a list of integers.
     * Only need one message instead of two as with the nodes */
    if (mpiInfo->rank == 0) { /* Master */
        for (;;) { /* Infinite loop */
            chunkEle = 0;
            for (i0 = 0; i0 < chunkSize; i0++) {
                if (totalEle >= numEle)
                    break;  /* End inner loop */
                int ret = fscanf(handle, "%d %d", &tempInts[i0*(2+numNodes)],
                        &tempInts[i0*(2+numNodes)+1]);
                if (ret!=2)
                    throw RipleyException("readMesh: Invalid file!");
                for (dim_t i1 = 0; i1 < numNodes; i1++) {
                    ret = fscanf(handle, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
                    if (ret!=1)
                        throw RipleyException("readMesh: Invalid file!");
                }
                ret = fscanf(handle, "\n");
                if (ret!=0)
                    throw RipleyException("readMesh: Invalid file!");
                totalEle++;
                chunkEle++;
            }
#ifdef ESYS_MPI
            /* Eventually we'll send chunk of elements to each CPU except
             * 0 itself, here goes one of them */
            if (nextCPU < mpiInfo->size) {
                tempInts[chunkSize * (2 + numNodes)] = chunkEle;
                MPI_Send(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81725, mpiInfo->comm);
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the
             * worker nodes plus one more chunk for the master */
            if (nextCPU > mpiInfo->size)
                break; /* End infinite loop */
        } // Infinite loop

    } else  { // Worker
#ifdef ESYS_MPI
        // each worker receives one message
        MPI_Status status;
        MPI_Recv(&tempInts[0], chunkSize*(2+numNodes)+1, MPI_INT, 0, 81725, mpiInfo->comm, &status);
        chunkEle = tempInts[chunkSize * (2 + numNodes)];
#endif
    } // Worker

    // copy element data from temporary array to mesh
    id.resize(chunkEle);
    tag.resize(chunkEle);
    color.assign(chunkEle, -1);
    owner.resize(chunkEle);
    nodesVec.resize(chunkEle*numNodes);
#pragma omp parallel for private (i0) schedule(static)
    for (i0 = 0; i0 < chunkEle; i0++) {
        id[i0] = tempInts[i0*(2+numNodes)];
        tag[i0] = tempInts[i0*(2+numNodes)+1];
        owner[i0] = mpiInfo->rank;
        for (dim_t i1 = 0; i1 < numNodes; i1++) {
            nodesVec[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
        }
    }
    points->swapEntries(id, tag, color, owner, nodesVec);
    tempInts.clear();
    dom->setPoints(points);

    /*************************  get the tag names ***************************/
    char *remainder = NULL, *ptr;
    size_t len = 0;
    int tag_key;
    if (mpiInfo->rank == 0) { /* Master */
        // Read the word 'Tag'
        if (!feof(handle)) {
            int ret = fscanf(handle, "%s\n", name);
            if (ret!=1)
                throw RipleyException("readMesh: Invalid file!");
        }
#ifdef _WIN32
        /* windows ftell lies on unix formatted text files */
        len = 0;
        while (1) {
            size_t malloc_chunk = 1024;
            size_t buff_size = 0;
            int ch;
            ch = fgetc(handle);
            if (ch == '\r')
                continue;
            if (len + 1 > buff_size) {
                TMPMEMREALLOC(remainder, remainder, buff_size + malloc_chunk, char);
            }
            if (ch == EOF) {
                /* hit EOF */
                remainder[len] = '\0'
                break;
            }
            remainder[len] = (char)ch;
            len++;
        }
#else
        /* Read rest of file in one chunk, after using seek to find length */
        {
            long cur_pos, end_pos;
            cur_pos = ftell(handle);
            fseek(handle, 0L, SEEK_END);
            end_pos = ftell(handle);
            fseek(handle, (long)cur_pos, SEEK_SET);
            remainder = TMPMEMALLOC(end_pos - cur_pos + 1, char);
            if (!feof(handle)) {
                fread(remainder, (size_t) end_pos-cur_pos, sizeof(char), handle);
                remainder[end_pos - cur_pos] = 0;
            }
        }
#endif
        len = strlen(remainder);
        while ((len > 1) && isspace(remainder[--len]))
            remainder[len] = '\0';
        len = strlen(remainder);
        TMPMEMREALLOC(remainder, remainder, len + 1, char);
    } /* Master */

#ifdef ESYS_MPI
    int len_i = (int)len;
    MPI_Bcast(&len_i, 1, MPI_INT, 0, mpiInfo->comm);
    len = (size_t) len_i;
    if (mpiInfo->rank != 0) {
        remainder = TMPMEMALLOC(len + 1, char);
        remainder[0] = '\0';
    }
    if (MPI_Bcast(remainder, len + 1, MPI_CHAR, 0, mpiInfo->comm) != MPI_SUCCESS)
        throw RipleyException("readMesh: broadcast of remainder failed");
#endif

    if (remainder[0]) {
        ptr = remainder;
        do {
            sscanf(ptr, "%s %d\n", name, &tag_key);
            if (*name)
                dom->setTagMap(name, tag_key);
            ptr++;
        } while (NULL != (ptr = strchr(ptr, '\n')) && *ptr);
    }
    if (remainder)
        TMPMEMFREE(remainder);

    /* close file */
    if (mpiInfo->rank == 0)
        fclose(handle);

    dom->prepare(optimize);

    Esys_MPIInfo_free(mpiInfo);
    blocktimer_increment("ReadMesh()", blocktimer_start);
    return dom->getPtr();
}

escript::Domain_ptr readGmsh(const string& filename, int numDim, bool optimize)
{
    if (filename.size() == 0)
        throw RipleyException("readGmsh: Empty file name!");

    double blocktimer_start = blocktimer_time();
    double version;
    const int MAX_NN=20;
    dim_t numNodes, totalNumElements = 0, numNodesPerElement = 0;
    index_t numElements = 0, numFaceElements = 0;
    char line[LenString_MAX + 1];

    Esys_MPIInfo *mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    checkPasoError();

    if (mpiInfo->size > 1) {
        Esys_MPIInfo_free(mpiInfo);
        throw RipleyException("reading gmsh with MPI is not supported yet.");
    }

    /* open file */
    FILE *handle = fopen(filename.c_str(), "r");
    if (handle == NULL) {
        Esys_MPIInfo_free(mpiInfo);
        stringstream msg;
        msg << "Opening gmsh file " << filename << " for reading failed.";
        throw RipleyException(msg.str());
    }

    /* allocate mesh */
    RipleyDomain *dom=new RipleyDomain(filename, numDim, mpiInfo);

    /* start reading */
    try {
        while (1) {
            /* find line starting with $ */
            do {
                if (!fgets(line, sizeof(line), handle))
                    break;
            } while (line[0] != '$');

            if (feof(handle))
                break;

            /* format */
            if (!strncmp(&line[1], "MeshFormat", 10)) {
                int format, size;
                int ret = fscanf(handle, "%lf %d %d\n", &version, &format, &size);
                if (ret!=3)
                    throw RipleyException("readGmsh: Invalid gmsh file!");
            }
            /* read nodes */
            if (!strncmp(&line[1], "NOD", 3) || !strncmp(&line[1], "NOE", 3) || !strncmp(&line[1], "Nodes", 5)) {

                int ret = fscanf(handle, "%d", &numNodes);
                if (ret!=1)
                    throw RipleyException("readGmsh: Invalid gmsh file!");
                NodeFile_ptr nodes(new NodeFile(numDim, mpiInfo));
                IndexVector id(numNodes);
                IndexVector gDOF(numNodes);
                IndexVector tag(numNodes);
                vector<double> coordinates(numNodes*numDim);
                for (index_t i = 0; i < numNodes; i++) {
                    if (1 == numDim) {
                        double rtmp0, rtmp1;
                        ret = fscanf(handle, "%d %le %le %le\n",
                                &id[i],
                                &coordinates[INDEX2(0, i, numDim)],
                                &rtmp0, &rtmp1);
                        if (ret!=4)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                    } else if (2 == numDim) {
                        double rtmp;
                        ret = fscanf(handle, "%d %le %le %le\n",
                                &id[i],
                                &coordinates[INDEX2(0, i, numDim)],
                                &coordinates[INDEX2(1, i, numDim)],
                                &rtmp);
                        if (ret!=4)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                    } else if (3 == numDim) {
                        ret = fscanf(handle, "%d %le %le %le\n",
                                &id[i],
                                &coordinates[INDEX2(0, i, numDim)],
                                &coordinates[INDEX2(1, i, numDim)],
                                &coordinates[INDEX2(2, i, numDim)]);
                        if (ret!=4)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                    }
                    gDOF[i] = id[i];
                    tag[i] = 0;
                }
                nodes->swapEntries(id, tag, gDOF, coordinates);
                dom->setNodes(nodes);
            }
            /* read elements */
            else if (!strncmp(&line[1], "ELM", 3) || !strncmp(&line[1], "Elements", 8)) {
                ElementTypeId final_element_type = InvalidElementType;
                ElementTypeId final_face_element_type = InvalidElementType;
                numElements = 0;
                numFaceElements = 0;
                int ret = fscanf(handle, "%d", &totalNumElements);
                if (ret!=1)
                    throw RipleyException("readGmsh: Invalid gmsh file!");

                IndexVector id(totalNumElements);
                IndexVector tag(totalNumElements);
                vector<ElementTypeId> element_type(totalNumElements);
                IndexVector vertices(totalNumElements*MAX_NN);

                /* read all elements */
                for (index_t e = 0; e < totalNumElements; e++) {
                    int gmsh_type;
                    dim_t element_dim = 0;
                    ret = fscanf(handle, "%d %d", &id[e], &gmsh_type);
                    if (ret!=1)
                        throw RipleyException("readGmsh: Invalid gmsh file!");
                    switch (gmsh_type) {
                        case 1: /* line order 1 */
                            element_type[e] = Line2;
                            element_dim = 1;
                            numNodesPerElement = 2;
                            break;
                        case 3: /* quad order 1 */
                            element_type[e] = Rec4;
                            numNodesPerElement = 4;
                            element_dim = 2;
                            break;
                        case 5: /* hexahedron order 1 */
                            element_type[e] = Hex8;
                            numNodesPerElement = 8;
                            element_dim = 3;
                            break;
                        case 15: /* point */
                            element_type[e] = Point1;
                            numNodesPerElement = 1;
                            element_dim = 0;
                            break;
                        default:
                            element_type[e] = InvalidElementType;
                            stringstream msg;
                            msg << "Unexpected/Unsupported gmsh element type "
                                << gmsh_type << " in mesh file " << filename;
                            throw RipleyException(msg.str());
                    }
                    if (element_dim == numDim) {
                        if (final_element_type == InvalidElementType) {
                            final_element_type = element_type[e];
                        } else if (final_element_type != element_type[e]) {
                            throw RipleyException(
                                    "Ripley can handle a single type of internal elements only.");
                        }
                        numElements++;
                    } else if (element_dim == numDim - 1) {
                        if (final_face_element_type == InvalidElementType) {
                            final_face_element_type = element_type[e];
                        } else if (final_face_element_type != element_type[e]) {
                            throw RipleyException(
                                    "Ripley can handle a single type of face elements only.");
                        }
                        numFaceElements++;
                    }

                    if (version <= 1.0) {
                        dim_t nodesPerEl;
                        index_t elId;
                        ret = fscanf(handle, "%d %d %d", &tag[e], &elId, &nodesPerEl);
                        if (ret!=3)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                        if (nodesPerEl != numNodesPerElement) {
                            stringstream msg;
                            msg << "Illegal number of nodes for element " << id[e]
                                << " in mesh file " << filename;
                            throw RipleyException(msg.str());
                        }
                    } else {
                        dim_t numTags;
                        ret = fscanf(handle, "%d", &numTags);
                        if (ret!=1)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                        tag[e] = 1;
                        for (index_t j = 0; j < numTags; j++) {
                            index_t itmp;
                            ret = fscanf(handle, "%d", &itmp);
                            if (ret!=1)
                                throw RipleyException("readGmsh: Invalid gmsh file!");
                            if (j == 0)
                                tag[e] = itmp;
                            /* ignore any other tags */
                        }
                    }
                    for (index_t j = 0; j < numNodesPerElement; j++) {
                        ret = fscanf(handle, "%d", &vertices[INDEX2(j, e, MAX_NN)]);
                        if (ret!=1)
                            throw RipleyException("readGmsh: Invalid gmsh file!");
                    }
                }

                /* all elements have been read, now we have to identify the
                 * elements for ripley */
                if (final_element_type == InvalidElementType) {
                    if (numDim == 1) {
                        final_element_type = Line2;
                    } else if (numDim == 2) {
                        final_element_type = Rec4;
                    } else if (numDim == 3) {
                        final_element_type = Hex8;
                    }
                }
                if (final_face_element_type == InvalidElementType) {
                    if (numDim == 1) {
                        final_face_element_type = Point1;
                    } else if (numDim == 2) {
                        final_face_element_type = Line2;
                    } else if (numDim == 3) {
                        final_face_element_type = Rec4;
                    }
                }
                ElementFile_ptr elements(new ElementFile(final_element_type, mpiInfo));
                ElementFile_ptr faceElements(new ElementFile(final_face_element_type, mpiInfo));
                ElementFile_ptr points(new ElementFile(Point1, mpiInfo));
                IndexVector eId, eTag, eNodes;
                IndexVector fId, fTag, fNodes;
                numElements = 0;
                numFaceElements = 0;
                for (index_t e = 0; e < totalNumElements; e++) {
                    if (element_type[e] == final_element_type) {
                        eId.push_back(id[e]);
                        eTag.push_back(tag[e]);
                        for (index_t j = 0; j < elements->getNumNodes(); ++j) {
                            eNodes.push_back(vertices[INDEX2(j, e, MAX_NN)]);
                        }
                        numElements++;
                    } else if (element_type[e] == final_face_element_type) {
                        fId.push_back(id[e]);
                        fTag.push_back(tag[e]);
                        for (index_t j = 0; j < faceElements->getNumNodes(); ++j) {
                            fNodes.push_back(vertices[INDEX2(j, e, MAX_NN)]);
                        }
                        numFaceElements++;
                    }
                }
                RankVector eOwner(numElements, 0);
                RankVector fOwner(numFaceElements, 0);
                IndexVector eColor(numElements, -1);
                IndexVector fColor(numFaceElements, -1);
                elements->swapEntries(eId, eTag, eColor, eOwner, eNodes);
                faceElements->swapEntries(fId, fTag, fColor, fOwner, fNodes);
                dom->setElements(elements);
                dom->setFaceElements(faceElements);
                dom->setPoints(points);
            }
            /* search for end of data block */
            do {
                if (!fgets(line, sizeof(line), handle)) {
                    stringstream msg;
                    msg << "Unexpected end of file in " << filename;
                    throw RipleyException(msg.str());
                }
                if (feof(handle)) {
                    stringstream msg;
                    msg << "Unexpected end of file in " << filename;
                    throw RipleyException(msg.str());
                }
            } while (line[0] != '$');

        } // endless loop

        dom->prepare(optimize);

    } catch (RipleyException(e)) {
        fclose(handle);
        Esys_MPIInfo_free(mpiInfo);
        throw e;
    }

    fclose(handle);
    Esys_MPIInfo_free(mpiInfo);
    blocktimer_increment("ReadGmsh()", blocktimer_start);
    return dom->getPtr();
}

escript::Domain_ptr brick(int n0, int n1, int n2, double l0, double l1, double l2, bool optimize)
{
    Esys_MPIInfo *mpiInfo = NULL;
    dim_t N0, N1, N2, NE0, NE1, NE2, Nstride0 = 0, Nstride1 = 0, Nstride2 = 0,
          local_NE0, local_NE1, local_NE2, local_N0 = 0, local_N1 = 0,
          local_N2 = 0, totalNECount, faceNECount, NDOF0 = 0, NDOF1 = 0,
          NDOF2 = 0;
    index_t myRank, offset0 = 0, offset1 = 0, offset2 = 0;
    const int DIM = 3;
    const int LEFTTAG = 1;     /* boundary x1=0 */
    const int RIGHTTAG = 2;    /* boundary x1=1 */
    const int BOTTOMTAG = 100; /* boundary x3=1 */
    const int TOPTAG = 200;    /* boundary x3=0 */
    const int FRONTTAG = 10;   /* boundary x2=0 */
    const int BACKTAG = 20;    /* boundary x2=1 */

    mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    checkPasoError();

    myRank=mpiInfo->rank;

    /* set up the global dimensions of the mesh */
    NE0=MAX(1, n0);
    NE1=MAX(1, n1);
    NE2=MAX(1, n2);
    N0 = NE0 + 1;
    N1 = NE1 + 1;
    N2 = NE2 + 1;

    /* allocate mesh */
    stringstream name;
    name << "Rectangular " << N0 << " x " << N1 << " x " << N2 << " mesh";
    RipleyDomain *dom=new RipleyDomain(name.str(), DIM, mpiInfo);

    try {
        dim_t NFaceElements = 0;

        // work out the largest dimension
        if (N2 == MAX3(N0, N1, N2)) {
            Nstride0 = 1;
            Nstride1 = N0;
            Nstride2 = N0*N1;
            local_NE0 = NE0;
            local_NE1 = NE1;
            Esys_MPIInfo_Split(mpiInfo, NE2, &local_NE2, &offset2);
        } else if (N1 == MAX3(N0,N1,N2)) {
            Nstride0 = N2;
            Nstride1 = N0*N2;
            Nstride2 = 1;
            local_NE0 = NE0;
            Esys_MPIInfo_Split(mpiInfo, NE1, &local_NE1, &offset1);
            local_NE2 = NE2;
        } else {
            Nstride0 = N1*N2;
            Nstride1 = 1;
            Nstride2 = N1;
            Esys_MPIInfo_Split(mpiInfo, NE0, &local_NE0, &offset0);
            local_NE1 = NE1;
            local_NE2 = NE2;
        }
        local_N0 = local_NE0>0 ? local_NE0+1 : 0;
        local_N1 = local_NE1>0 ? local_NE1+1 : 0;
        local_N2 = local_NE2>0 ? local_NE2+1 : 0;
        dim_t local_NN = local_N0*local_N1*local_N2;

        /* get the number of surface elements */
        if (local_NE2 > 0) {
            NDOF2=N2;
            if (offset2 == 0)
                NFaceElements+=local_NE1*local_NE0;
            if (local_NE2+offset2 == NE2)
                NFaceElements+=local_NE1*local_NE0;
        } else {
            NDOF2=N2-1;
        }

        if (local_NE0 > 0) {
            NDOF0 = N0;
            if (offset0 == 0)
                NFaceElements+=local_NE1*local_NE2;
            if (local_NE0+offset0 == NE0)
                NFaceElements+=local_NE1*local_NE2;
        } else {
            NDOF0=N0-1;
        }
        if (local_NE1 > 0) {
            NDOF1 = N1;
            if (offset1 == 0)
                NFaceElements+=local_NE0*local_NE2;
            if (local_NE1+offset1 == NE1)
                NFaceElements+=local_NE0*local_NE2;
        } else {
            NDOF1=N1-1;
        }

        /*  allocate tables */
        NodeFile_ptr nodes(new NodeFile(DIM, mpiInfo));
        ElementFile_ptr points(new ElementFile(Point1, mpiInfo));
        ElementFile_ptr faceElements(new ElementFile(Rec4, mpiInfo));
        ElementFile_ptr elements(new ElementFile(Hex8, mpiInfo));
        IndexVector id(local_NN);
        IndexVector gDOF(local_NN);
        IndexVector tag(local_NN);
        vector<double> coordinates(local_NN*DIM);

        dim_t i0, i1, i2;
        /* create nodes */
#pragma omp parallel for private(i0,i1,i2)
        for (i2=0; i2<local_N2; i2++) {
            for (i1=0; i1<local_N1; i1++) {
                for (i0=0; i0<local_N0; i0++) {
                    const dim_t k=i0+local_N0*i1+local_N0*local_N1*i2;
                    const dim_t global_i0=i0+offset0;
                    const dim_t global_i1=i1+offset1;
                    const dim_t global_i2=i2+offset2;
                    coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*l0;
                    coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*l1;
                    coordinates[INDEX2(2,k,DIM)]=DBLE(global_i2)/DBLE(N2-1)*l2;
                    id[k] = Nstride0*global_i0 + Nstride1*global_i1 + Nstride2*global_i2;
                    tag[k] = 0;
                    gDOF[k] = Nstride0*(global_i0 % NDOF0) +
                              Nstride1*(global_i1 % NDOF1) +
                              Nstride2*(global_i2 % NDOF2);
                }
            }
        }
        nodes->swapEntries(id, tag, gDOF, coordinates);
        dom->setNodes(nodes);

        // create elements
        const dim_t localNE = local_NE0 * local_NE1 * local_NE2;
        dim_t NN = elements->getNumNodes(); // =8
        id.resize(localNE);
        tag.assign(localNE, 0);
        IndexVector color(localNE, -1);
        IndexVector nodesVec(localNE*NN);
        RankVector owner(localNE, myRank);

#pragma omp parallel for private(i0, i1, i2)
        for (i2=0; i2<local_NE2; i2++) {
            for (i1=0; i1<local_NE1; i1++) {
                for (i0=0; i0<local_NE0; i0++) {
                    const dim_t k=i0 + local_NE0*i1 + local_NE0*local_NE1*i2;
                    const dim_t node0=Nstride0*(i0+offset0)+Nstride1*(i1+offset1)+Nstride2*(i2+offset2);
                    id[k]=(i0+offset0)+NE0*(i1+offset1)+NE0*NE1*(i2+offset2);
                    nodesVec[INDEX2(0,k,NN)]=node0;
                    nodesVec[INDEX2(1,k,NN)]=node0+Nstride0;
                    nodesVec[INDEX2(2,k,NN)]=node0+Nstride1+Nstride0;
                    nodesVec[INDEX2(3,k,NN)]=node0+Nstride1;
                    nodesVec[INDEX2(4,k,NN)]=node0+Nstride2;
                    nodesVec[INDEX2(5,k,NN)]=node0+Nstride2+Nstride0;
                    nodesVec[INDEX2(6,k,NN)]=node0+Nstride2+Nstride1+Nstride0;
                    nodesVec[INDEX2(7,k,NN)]=node0+Nstride2+Nstride1;
                }
            }
        }
        elements->swapEntries(id, tag, color, owner, nodesVec);
        dom->setElements(elements);

        // create face elements
        NN = faceElements->getNumNodes(); // =4
        id.resize(NFaceElements);
        tag.resize(NFaceElements);
        color.assign(NFaceElements, -1);
        nodesVec.resize(NFaceElements*NN);
        owner.assign(NFaceElements, myRank);

        totalNECount = NE0*NE1*NE2;
        faceNECount = 0;
        // quadrilateral elements on boundary 1 (x3=0):
        if (local_NE2 > 0) {
            // boundary 100 (x3=0):
            if (offset2==0) {
#pragma omp parallel for private(i0, i1)
                for (i1=0; i1<local_NE1; i1++) {
                    for (i0=0; i0<local_NE0; i0++) {
                        const dim_t k=i0+local_NE0*i1+faceNECount;
                        const dim_t node0=Nstride0*(i0+offset0)+Nstride1*(i1+offset1);
                        id[k]=(i0+offset0)+NE0*(i1+offset1)+totalNECount;
                        tag[k]=BOTTOMTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride1;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride1+Nstride0;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride0;
                    }
                }
                faceNECount+=local_NE1*local_NE0;
            }
            totalNECount+=NE1*NE0;
            /* boundary 200 (x3=1): */
            if (local_NE2+offset2 == NE2) {
#pragma omp parallel for private(i0, i1)
                for (i1=0; i1<local_NE1; i1++) {
                    for (i0=0; i0<local_NE0; i0++) {
                        const dim_t k=i0+local_NE0*i1+faceNECount;
                        const dim_t node0=Nstride0*(i0+offset0)+Nstride1*(i1+offset1)+Nstride2*(NE2-1);
                        id[k]=(i0+offset0)+NE0*(i1+offset1)+totalNECount;
                        tag[k]=TOPTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0+Nstride2;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride2+Nstride0;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride2+Nstride1+Nstride0;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride2+Nstride1;
                    }
                }
                faceNECount+=local_NE1*local_NE0;
            }
            totalNECount+=NE1*NE0;
        }
        if (local_NE0>0) {
            /* boundary 001 (x1=0): */
            if (offset0 == 0) {
#pragma omp parallel for private(i1, i2)
                for (i2=0; i2<local_NE2; i2++) {
                    for (i1=0; i1<local_NE1; i1++) {
                        const dim_t k=i1+local_NE1*i2+faceNECount;
                        const dim_t node0=Nstride1*(i1+offset1)+Nstride2*(i2+offset2);
                        id[k]=(i1+offset1)+NE1*(i2+offset2)+totalNECount;
                        tag[k]=LEFTTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride2;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride2+Nstride1;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride1;
                    }
                }
                faceNECount+=local_NE1*local_NE2;
            }
            totalNECount+=NE1*NE2;
            /* boundary 002 (x1=1): */
            if (local_NE0+offset0 == NE0) {
#pragma omp parallel for private(i1, i2)
                for (i2=0; i2<local_NE2; i2++) {
                    for (i1=0; i1<local_NE1; i1++) {
                        const dim_t k=i1+local_NE1*i2+faceNECount;
                        const dim_t node0=Nstride0*(NE0-1)+Nstride1*(i1+offset1)+Nstride2*(i2+offset2);
                        id[k]=(i1+offset1)+NE1*(i2+offset2)+totalNECount;
                        tag[k]=RIGHTTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0+Nstride0;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride1+Nstride0;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride2+Nstride1+Nstride0;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride2+Nstride0;
                    }
                }
                faceNECount+=local_NE1*local_NE2;
            }
            totalNECount+=NE1*NE2;
        }
        if (local_NE1>0) {
            /* boundary 010 (x2=0): */
            if (offset1 == 0) {
#pragma omp parallel for private(i0, i2)
                for (i2=0; i2<local_NE2; i2++) {
                    for (i0=0; i0<local_NE0; i0++) {
                        const dim_t k=i0+local_NE0*i2+faceNECount;
                        const dim_t node0=Nstride0*(i0+offset0)+Nstride2*(i2+offset2);
                        id[k]=(i2+offset2)+NE2*(offset0+i0)+totalNECount;
                        tag[k]=FRONTTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride0;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride2+Nstride0;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride2;
                    }
                }
                faceNECount+=local_NE0*local_NE2;
            }
            totalNECount+=NE0*NE2;
            /* boundary 020 (x2=1): */
            if (local_NE1+offset1 == NE1) {
#pragma omp parallel for private(i0, i2)
                for (i2=0; i2<local_NE2; i2++) {
                    for (i0=0; i0<local_NE0; i0++) {
                        const dim_t k=i0+local_NE0*i2+faceNECount;
                        const dim_t node0=Nstride0*(i0+offset0)+Nstride1*(NE1-1)+Nstride2*(i2+offset2);
                        id[k]=(i2+offset2)+NE2*(i0+offset0)+totalNECount;
                        tag[k]=BACKTAG;
                        nodesVec[INDEX2(0,k,NN)]=node0+Nstride1;
                        nodesVec[INDEX2(1,k,NN)]=node0+Nstride2+Nstride1;
                        nodesVec[INDEX2(2,k,NN)]=node0+Nstride2+Nstride1+Nstride0;
                        nodesVec[INDEX2(3,k,NN)]=node0+Nstride1+Nstride0;
                    }
                }
                faceNECount+=local_NE0*local_NE2;
            }
            //totalNECount+=NE0*NE2;
        }
        faceElements->swapEntries(id, tag, color, owner, nodesVec);
        dom->setFaceElements(faceElements);
        dom->setPoints(points);

        // add tag names
        dom->setTagMap("top", TOPTAG);
        dom->setTagMap("bottom", BOTTOMTAG);
        dom->setTagMap("left", LEFTTAG);
        dom->setTagMap("right", RIGHTTAG);
        dom->setTagMap("front", FRONTTAG);
        dom->setTagMap("back", BACKTAG);

        // prepare mesh for further calculations
        dom->prepare(optimize);

    } catch (RipleyException(e)) {
        Esys_MPIInfo_free(mpiInfo);
        throw e;
    }

    Esys_MPIInfo_free(mpiInfo);
    return dom->getPtr();
}

escript::Domain_ptr rectangle(int n0, int n1, double l0, double l1, bool optimize)
{
    const int DIM = 2;
    const int LEFTTAG = 1;    // boundary x1=0
    const int RIGHTTAG = 2;   // boundary x1=1
    const int BOTTOMTAG = 10; // boundary x2=0
    const int TOPTAG = 20;    // boundary x2=1

    Esys_MPIInfo *mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
    checkPasoError();

    // determine global number of elements and nodes for the mesh
    const dim_t NE0 = MAX(1, n0);
    const dim_t NE1 = MAX(1, n1);
    const dim_t N0 = NE0 + 1;
    const dim_t N1 = NE1 + 1;

    // allocate domain
    stringstream name;
    name << "Rectangular " << N0 << " x " << N1 << " mesh";
    RipleyDomain *dom = new RipleyDomain(name.str(), DIM, mpiInfo);

    try {
        dim_t NFaceElements = 0;
        dim_t local_NE0, local_NE1;
        index_t offset0, offset1;

        // work out the largest dimension
        if (N1 == MAX(N0,N1)) {
            Esys_MPIInfo_Split(mpiInfo, NE1, &local_NE1, &offset1);
            local_NE0 = NE0;
            offset0 = 0;
        } else {
            Esys_MPIInfo_Split(mpiInfo, NE0, &local_NE0, &offset0);
            local_NE1 = NE1;
            offset1 = 0;
        }
        dim_t local_N0 = local_NE0>0 ? local_NE0+1 : 0;
        dim_t local_N1 = local_NE1>0 ? local_NE1+1 : 0;
        dim_t local_NN = local_N0*local_N1;
        dim_t NDOF0, NDOF1;

        // get the number of surface elements
        if (local_NE0 > 0) {
            NDOF0 = N0;
            if (offset0 == 0)
                NFaceElements += local_NE1;
            if (local_NE0+offset0 == NE0)
                NFaceElements += local_NE1;
        } else {
            NDOF0 = N0-1;
        }
        if (local_NE1 > 0) {
            NDOF1 = N1;
            if (offset1 == 0)
                NFaceElements+=local_NE0;
            if (local_NE1+offset1 == NE1)
                NFaceElements+=local_NE0;
        } else {
            NDOF1 = N1-1;
        }

        // allocate tables
        NodeFile_ptr nodes(new NodeFile(DIM, mpiInfo));
        ElementFile_ptr elements(new ElementFile(Rec4, mpiInfo));
        ElementFile_ptr faceElements(new ElementFile(Line2, mpiInfo));
        ElementFile_ptr points(new ElementFile(Point1, mpiInfo));
        IndexVector id(local_NN);
        IndexVector gDOF(local_NN);
        IndexVector tag(local_NN);
        vector<double> coordinates(local_NN*DIM);

        dim_t i0, i1;

        /* create nodes */
#pragma omp parallel for private(i0,i1)
        for (i1=0; i1<local_N1; i1++) {
            for (i0=0; i0<local_N0; i0++) {
                const dim_t k=i0+local_N0*i1;
                const dim_t global_i0=i0+offset0;
                const dim_t global_i1=i1+offset1;
                coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*l0;
                coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*l1;
                id[k] = global_i0 + N0*global_i1;
                tag[k] = 0;
                gDOF[k] = (global_i0 % NDOF0) + N0*(global_i1 % NDOF1);
            }
        }
        nodes->swapEntries(id, tag, gDOF, coordinates);
        dom->setNodes(nodes);

        // create elements
        const dim_t localNE = local_NE0 * local_NE1;
        dim_t NN = elements->getNumNodes(); // =4
        id.resize(localNE);
        tag.assign(localNE, 0);
        IndexVector color(localNE, -1);
        IndexVector nodesVec(localNE*NN);
        RankVector owner(localNE, mpiInfo->rank);
#pragma omp parallel for private(i0,i1)
        for (i1=0; i1<local_NE1; i1++) {
            for (i0=0; i0<local_NE0; i0++) {
                const dim_t k=i0+local_NE0*i1;
                const index_t node0=i0+offset0 + N0*(i1+offset1);

                id[k]=(i0+offset0)+NE0*(i1+offset1);
                nodesVec[INDEX2(0,k,NN)]=node0;
                nodesVec[INDEX2(1,k,NN)]=node0+1;
                nodesVec[INDEX2(2,k,NN)]=node0+N0+1;
                nodesVec[INDEX2(3,k,NN)]=node0+N0;
            }
        }
        elements->swapEntries(id, tag, color, owner, nodesVec);
        dom->setElements(elements);

        // create face elements
        NN = faceElements->getNumNodes(); // =2
        id.resize(NFaceElements);
        tag.resize(NFaceElements);
        color.assign(NFaceElements, -1);
        nodesVec.resize(NFaceElements*NN);
        owner.assign(NFaceElements, mpiInfo->rank);

        dim_t totalNECount = NE0*NE1;
        dim_t faceNECount = 0;
        if (local_NE0 > 0) {
            // boundary 001 (x1=0):
            if (offset0 == 0) {
#pragma omp parallel for private(i1)
                for (i1=0; i1<local_NE1; i1++) {
                    const dim_t k=i1+faceNECount;
                    const dim_t node0=N0*(i1+offset1);

                    id[k]=i1+offset1+totalNECount;
                    tag[k]=LEFTTAG;
                    nodesVec[INDEX2(0,k,NN)]=node0+N0;
                    nodesVec[INDEX2(1,k,NN)]=node0;
                }
                faceNECount+=local_NE1;
            }
            totalNECount+=NE1;

            // boundary 002 (x1=1):
            if (local_NE0+offset0 == NE0) {
#pragma omp parallel for private(i1)
                for (i1=0; i1<local_NE1; i1++) {
                    const dim_t k=i1+faceNECount;
                    const dim_t node0=NE0-1+N0*(i1+offset1);

                    id[k]=(i1+offset1)+totalNECount;
                    tag[k]=RIGHTTAG;
                    nodesVec[INDEX2(0,k,NN)]=node0+1;
                    nodesVec[INDEX2(1,k,NN)]=node0+N0+1;
                }
                faceNECount+=local_NE1;
            }
            totalNECount+=NE1;
        }
        if (local_NE1 > 0) {
            // boundary 010 (x2=0):
            if (offset1 == 0) {
#pragma omp parallel for private(i0)
                for (i0=0; i0<local_NE0; i0++) {
                    const dim_t k=i0+faceNECount;
                    const dim_t node0=i0+offset0;

                    id[k]=offset0+i0+totalNECount;
                    tag[k]=BOTTOMTAG;
                    nodesVec[INDEX2(0,k,NN)]=node0;
                    nodesVec[INDEX2(1,k,NN)]=node0+1;
                }
                faceNECount+=local_NE0;
            }
            totalNECount+=NE0;

            // boundary 020 (x2=1):
            if (local_NE1+offset1 == NE1) {
#pragma omp parallel for private(i0)
                for (i0=0; i0<local_NE0; i0++) {
                    const dim_t k=i0+faceNECount;
                    const dim_t node0=i0+offset0+N0*(NE1-1);

                    id[k]=i0+offset0+totalNECount;
                    tag[k]=TOPTAG;
                    nodesVec[INDEX2(0,k,NN)]=node0+N0+1;
                    nodesVec[INDEX2(1,k,NN)]=node0+N0;
                }
                faceNECount+=local_NE0;
            }
            totalNECount+=NE0;
        }
        faceElements->swapEntries(id, tag, color, owner, nodesVec);
        dom->setFaceElements(faceElements);
        dom->setPoints(points);

        // add tag names
        dom->setTagMap("top", TOPTAG);
        dom->setTagMap("bottom", BOTTOMTAG);
        dom->setTagMap("left", LEFTTAG);
        dom->setTagMap("right", RIGHTTAG);

        // prepare mesh for further calculations
        dom->prepare(optimize);

    } catch (RipleyException(e)) {
        Esys_MPIInfo_free(mpiInfo);
        throw e;
    }

    Esys_MPIInfo_free(mpiInfo);
    return dom->getPtr();
}

} // end of namespace ripley

