
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
#include "FinleyException.h"

#include <escript/index.h>

#include <cstdio>

//can't return because the flag need to be shared across all nodes
#define SSCANF_CHECK(scan_ret) { if (scan_ret == EOF) { errorFlag = 1;} }
#define MAX_numNodes_gmsh 20

#define EARLY_EOF 1
#define MISSING_NODES 2
#define MISSING_ELEMENTS 3
#define THROW_ERROR 4
#define SUCCESS 5
#define ERROR 6
/*
    error flags include:

        0 - all ok
        1 - early EOF
        2 - EOF before nodes section found
        3 - EOF before elements section found
        4 - throw errorMsg
        5 - EOF at apropriate time.
        6 - !noError
*/

namespace {

using namespace finley;
using escript::IOError;

struct ElementInfo {
    ElementTypeId type;
    int id;
    int dim;
    int* vertex;
    int tag;
};

struct Msh4Entities {
    std::map<int,int> pointTags;
    std::map<int,int> curveTags;
    std::map<int,int> surfaceTags;
    std::map<int,int> volumeTags;
};

inline bool get_line(std::vector<char>& line, FILE* file)
{
    int capacity = 1024;
    line.clear();
    line.resize(capacity);
    char *tmp = &line[0];
    char *res = NULL;
    //not terribly efficient, but any line longer than 1024
    //is probably already bad
    while ((res = fgets(tmp, 1023, file)) == tmp
            && strchr(tmp, '\n') == NULL) {
        capacity += 1024;
        line.resize(capacity);
        tmp = strchr(&line[0], '\0'); //this bit is awful, O(n) instead of O(1)
        if (capacity > 1024) { //madness
            res = NULL;
            break;
        }
    }
    return res == tmp; //true if line read, false if EOF without or without \n
}

inline char* next_space(char** position, int count)
{
    for (int i = 0; i < count; i++) {
        *position = strchr(*position, ' ');
        if ((*position)++ == NULL) //move off the space
            return NULL;
    }
    return *position;
}

int getSingleElement(FILE* f, int dim, double version, struct ElementInfo& e,
        std::string& errorMsg, const std::string& filename,
        bool useMacroElements)
{
    int gmsh_type = -1;

    std::vector<char> line;
    if (!get_line(line, f))
        return EARLY_EOF;
    char *position = &line[0];
    if (sscanf(position, "%d %d", &e.id, &gmsh_type) != 2) {
        errorMsg = "malformed mesh file";
        return THROW_ERROR;
    }
    if (next_space(&position, 2) == NULL)
        return EARLY_EOF;

    int numNodesPerElement = 0;
    switch (gmsh_type) {
        case 1:  /* line order 1 */
            e.type = Line2;
            e.dim = 1;
            numNodesPerElement = 2;
            break;
        case 2:  /* triangle order 1 */
            e.type=Tri3;
            numNodesPerElement= 3;
            e.dim=2;
            break;
        case 3:  /* quadrilateral order 1 */
            e.type=Rec4;
            numNodesPerElement= 4;
            e.dim=2;
            break;
        case 4:  /* tetrahedron order 1 */
            e.type=Tet4;
            numNodesPerElement= 4;
            e.dim=3;
            break;
        case 5:  /* hexahedron order 1 */
            e.type=Hex8;
            numNodesPerElement= 8;
            e.dim=3;
            break;
        case 8:  /* line order 2 */
            if (useMacroElements) {
                e.type=Line3Macro;
            } else {
                e.type=Line3;
            }
            numNodesPerElement= 3;
            e.dim=1;
            break;
        case 9:  /* triangle order 2 */
            if (useMacroElements) {
                 e.type=Tri6Macro;
            } else {
                 e.type=Tri6;
            }
            numNodesPerElement= 6;
            e.dim=2;
            break;
        case 10:  /* quadrilateral order 2 */
            if (useMacroElements) {
                e.type=Rec9Macro;
            } else {
                e.type=Rec9;
            }
            numNodesPerElement= 9;
            e.dim=2;
            break;
        case 11:  /* tetrahedron order 2 */
            if (useMacroElements) {
                e.type=Tet10Macro;
            } else {
                e.type=Tet10;
            }
            numNodesPerElement= 10;
            e.dim=3;
            break;
        case 16:  /* rectangular order 2 */
            e.type=Rec8;
            numNodesPerElement= 8;
            e.dim=2;
            break;
        case 17:  /* hexahedron order 2 */
            e.type=Hex20;
            numNodesPerElement= 20;
            e.dim=3;
            break;
        case 15 :  /* point */
            e.type=Point1;
            numNodesPerElement= 1;
            e.dim=0;
            break;
        default:
            {
                e.type=NoRef;
                e.dim=-1;
                std::stringstream ss;
                ss << "readGmsh: Unexpected gmsh element type "
                    << gmsh_type << " in mesh file " << filename;
                errorMsg = ss.str();
                return THROW_ERROR;
            }
    }
    if (version <= 1.0) {
        int tmp = 0;
        if (sscanf(position, "%d %*d %d", &e.tag, &tmp) == 0
                || next_space(&position, 3) == NULL)
            return EARLY_EOF;
        if (tmp != numNodesPerElement) {
            std::stringstream ss;
            ss << "readGmsh: Illegal number of nodes for element " << e.id
                << " in mesh file " << filename;
            errorMsg = ss.str();
            return THROW_ERROR;
        }
    } else {
        e.tag = 1;
        int numTags = 0; //this is garbage and never used
        if (sscanf(position, "%d", &numTags) == 0
                || next_space(&position, 1) == NULL)
            return EARLY_EOF;
        if (sscanf(position, "%d", &e.tag) == 0
                || next_space(&position, numTags) == NULL)
            return EARLY_EOF;
        // ignore any other tags, second tag would be elementary id,
        // third tag would be partition id
    }

    for (int j = 0; j < numNodesPerElement; j++) {
        if (sscanf(position, "%d", e.vertex+j) == 0
                || next_space(&position, 1) == NULL)
            return EARLY_EOF;
    }
    return 0;
}

int getSingleElementMSH4(FILE* f, int dim, double version, struct ElementInfo& e,
        std::string& errorMsg, const std::string& filename,
        bool useMacroElements, int gmsh_type, char *position)
{
    int numNodesPerElement = 0;
    switch (gmsh_type) {
        case 1:  /* line order 1 */
            e.type = Line2;
            e.dim = 1;
            numNodesPerElement = 2;
            break;
        case 2:  /* triangle order 1 */
            e.type=Tri3;
            numNodesPerElement= 3;
            e.dim=2;
            break;
        case 3:  /* quadrilateral order 1 */
            e.type=Rec4;
            numNodesPerElement= 4;
            e.dim=2;
            break;
        case 4:  /* tetrahedron order 1 */
            e.type=Tet4;
            numNodesPerElement= 4;
            e.dim=3;
            break;
        case 5:  /* hexahedron order 1 */
            e.type=Hex8;
            numNodesPerElement= 8;
            e.dim=3;
            break;
        case 8:  /* line order 2 */
            if (useMacroElements) {
                e.type=Line3Macro;
            } else {
                e.type=Line3;
            }
            numNodesPerElement= 3;
            e.dim=1;
            break;
        case 9:  /* triangle order 2 */
            if (useMacroElements) {
                 e.type=Tri6Macro;
            } else {
                 e.type=Tri6;
            }
            numNodesPerElement= 6;
            e.dim=2;
            break;
        case 10:  /* quadrilateral order 2 */
            if (useMacroElements) {
                e.type=Rec9Macro;
            } else {
                e.type=Rec9;
            }
            numNodesPerElement= 9;
            e.dim=2;
            break;
        case 11:  /* tetrahedron order 2 */
            if (useMacroElements) {
                e.type=Tet10Macro;
            } else {
                e.type=Tet10;
            }
            numNodesPerElement= 10;
            e.dim=3;
            break;
        case 16:  /* rectangular order 2 */
            e.type=Rec8;
            numNodesPerElement= 8;
            e.dim=2;
            break;
        case 17:  /* hexahedron order 2 */
            e.type=Hex20;
            numNodesPerElement= 20;
            e.dim=3;
            break;
        case 15 :  /* point */
            e.type=Point1;
            numNodesPerElement= 1;
            e.dim=0;
            break;
        default:
            {
                e.type=NoRef;
                e.dim=-1;
                std::stringstream ss;
                ss << "readGmsh: Unexpected gmsh element type "
                    << gmsh_type << " in mesh file " << filename;
                errorMsg = ss.str();
                return THROW_ERROR;
            }
    }

    // char *position = &line[0];
    sscanf(position, "%d", &e.id);
    next_space(&position, 1);

    for (int j = 0; j < numNodesPerElement; j++) {
        if (sscanf(position, "%d", e.vertex+j) == 0
                || next_space(&position, 1) == NULL)
            return EARLY_EOF;
    }

    return 0;
}

int getElementsMaster(escript::JMPI& mpiInfo, FinleyDomain* dom,
                      FILE* fileHandle, std::string& errorMsg,
                      bool useMacroElements, const std::string& filename,
                      int numDim, double version, int order, int reducedOrder, Msh4Entities TagMap)
{
    /*
     *  This function should read in the elements and distribute
     *  them to the apropriate process.
     */
    int errorFlag = 0;
    ElementTypeId finalElementType = NoRef;
    ElementTypeId finalFaceElementType = NoRef;
    ElementTypeId contactElementType = NoRef;
    int numElements = 0, numFaceElements = 0, totalNumElements = 0;
    const_ReferenceElementSet_ptr refPoints, refContactElements;
    const_ReferenceElementSet_ptr refFaceElements, refElements;
    std::vector<char> line;
    if (!get_line(line, fileHandle))
        errorFlag = EARLY_EOF;
    int scan_ret, numEntityBlocks;
    if (version >= 4.1) {
        int minNodeTag, maxNodeTag;
        scan_ret = sscanf(&line[0], "%d %d %d %d", &numEntityBlocks, &totalNumElements, &minNodeTag, &maxNodeTag);
    } else if (version == 4.0){
        scan_ret = sscanf(&line[0], "%d %d", &numEntityBlocks, &totalNumElements);
    } else {
        scan_ret = sscanf(&line[0], "%d", &totalNumElements);
    }
    SSCANF_CHECK(scan_ret);

#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple MPI procs
    if (mpiInfo->size > 1) {
        int msg = totalNumElements;
        MPI_Bcast(&msg, 1, MPI_INT, 0, mpiInfo->comm);
    }
#endif

    int chunkSize = totalNumElements / mpiInfo->size;
    int rest = totalNumElements - (mpiInfo->size-1)*chunkSize;
    const size_t storage = std::max(chunkSize, rest);
    std::vector<int> id(storage+1);
    std::vector<int> tag(storage+1);
    std::vector<int> vertices(storage*MAX_numNodes_gmsh, -1);
    std::vector<ElementTypeId> elementType(storage+1);
    std::vector<int> elementIndices(storage, -1);
    std::vector<int> faceElementIndices(storage, -1);

#pragma omp parallel for schedule(static)
    for (int i=0; i<storage; i++) {
        id[i] = -1;
        tag[i] = -1;
        elementType[i] = NoRef;
    }

    int cpuId = 0;
    int chunkElements = 0;
    int chunkFaceElements = 0;
    int chunkOtherElements = 0;

    // read all in
    if(version >= 4.0){

        int count = 0;

        // Loop over the entity blocks
        for(index_t block_num = 0; block_num < numEntityBlocks; block_num++){

            if (!get_line(line, fileHandle))
                errorFlag = EARLY_EOF;

            int entityTag, entityDim, numElements, gmsh_type;

            if(version == 4.0)
                scan_ret = sscanf(&line[0], "%d %d %d %d", &entityTag, &entityDim, &gmsh_type, &numElements);
            else
                scan_ret = sscanf(&line[0], "%d %d %d %d", &entityDim, &entityTag, &gmsh_type, &numElements);            

            // Loop over elements in each entity block
            for (index_t e = 0; e < numElements; e++, count++) {

                if (cpuId >= mpiInfo->size-1) {
                    chunkSize = rest;
                }

                struct ElementInfo element = {NoRef, 0, 0, &vertices[count*MAX_numNodes_gmsh], 0};

                get_line(line, fileHandle);
                getSingleElementMSH4(fileHandle, numDim, version, element,
                                        errorMsg, filename, useMacroElements, gmsh_type, &line[0]);
                if(element.dim == 0)
                    element.tag = TagMap.pointTags.find(entityTag)->second;
                else if (element.dim == 1)
                    element.tag = TagMap.curveTags.find(entityTag)->second;
                else if (element.dim == 2)
                    element.tag = TagMap.surfaceTags.find(entityTag)->second;
                else if (element.dim == 3)
                    element.tag = TagMap.volumeTags.find(entityTag)->second;

                elementType[count] = element.type;
                id[count] = element.id;
                tag[count] = element.tag;

                // for tet10 the last two nodes need to be swapped
                if (element.type == Tet10 || element.type == Tet10Macro) {
                    int vertex = vertices[INDEX2(9, count, MAX_numNodes_gmsh)];
                    vertices[INDEX2(9, count, MAX_numNodes_gmsh)] = vertices[INDEX2(8, count, MAX_numNodes_gmsh)];
                    vertices[INDEX2(8, count, MAX_numNodes_gmsh)] = vertex;
                }

                if (element.dim == numDim) {
                    if (finalElementType == NoRef) {
                       finalElementType = element.type;
                    } else if (finalElementType != element.type) {
                        errorMsg = "Finley can only handle a single type of internal elements.";
                        errorFlag = THROW_ERROR;
                    }
                    elementIndices[chunkElements] = count;
                    chunkElements++;
                } else if (element.dim == numDim-1) {
                    if (finalFaceElementType == NoRef) {
                       finalFaceElementType = element.type;
                    } else if (finalFaceElementType != element.type) {
                       errorMsg = "Finley can only handle a single type of face elements.";
                       errorFlag = THROW_ERROR;
                    }
                    faceElementIndices[chunkFaceElements] = count;
                    numFaceElements++;
                    chunkFaceElements++;
                } else {
                    chunkOtherElements++;
                }
        #ifdef ESYS_MPI
                if (count < chunkSize - 1)
                    continue;

                // the last chunk is left for the master process
                if (++cpuId >= mpiInfo->size) {
                    continue;
                }

                if (errorFlag) {
                    for(; cpuId < mpiInfo->size; cpuId++) {
                        MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpiInfo->comm);
                    }
                    break;
                }
                int chunkInfo[2];
                chunkInfo[0] = chunkElements;
                chunkInfo[1] = chunkFaceElements;

                MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpiInfo->comm);
                MPI_Send(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, cpuId, 81720, mpiInfo->comm);
                MPI_Send(&id[0], chunkSize, MPI_INT, cpuId, 81721, mpiInfo->comm);
                MPI_Send(&tag[0], chunkSize, MPI_INT, cpuId, 81722, mpiInfo->comm);
                MPI_Send(&elementType[0], chunkSize, MPI_INT, cpuId, 81723, mpiInfo->comm);
                MPI_Send(chunkInfo, 2, MPI_INT, cpuId, 81724, mpiInfo->comm);
                MPI_Send(&elementIndices[0], chunkElements, MPI_INT, cpuId, 81725, mpiInfo->comm);
                MPI_Send(&faceElementIndices[0], chunkFaceElements, MPI_INT, cpuId, 81726, mpiInfo->comm);

                // reset arrays for next cpu
        #pragma omp parallel for schedule(static)
                for (index_t i = 0; i < chunkSize*MAX_numNodes_gmsh; i++)
                    vertices[i] = -1;
        #pragma omp parallel for schedule(static)
                for (index_t i = 0; i < chunkSize; i++) {
                    id[i] = -1;
                    tag[i] = -1;
                    elementType[i] = NoRef;
                }
                chunkElements = 0;
                chunkFaceElements = 0;
                chunkOtherElements = 0;
                count = -1;
        #endif
            }
        }

    } else { // Version < 4.0
        for (index_t e = 0, count = 0; e < totalNumElements; e++, count++) {
            if (cpuId >= mpiInfo->size-1) {
                chunkSize = rest;
            }

            struct ElementInfo element = {NoRef, 0, 0, &vertices[count*MAX_numNodes_gmsh], 0};
            getSingleElement(fileHandle, numDim, version, element, errorMsg,
                             filename, useMacroElements);
            elementType[count] = element.type;
            id[count] = element.id;
            tag[count] = element.tag;

            // for tet10 the last two nodes need to be swapped
            if (element.type == Tet10 || element.type == Tet10Macro) {
                int vertex = vertices[INDEX2(9, count, MAX_numNodes_gmsh)];
                vertices[INDEX2(9, count, MAX_numNodes_gmsh)] = vertices[INDEX2(8, count, MAX_numNodes_gmsh)];
                vertices[INDEX2(8, count, MAX_numNodes_gmsh)] = vertex;
            }

            if (element.dim == numDim) {
                if (finalElementType == NoRef) {
                   finalElementType = element.type;
                } else if (finalElementType != element.type) {
                    errorMsg = "Finley can only handle a single type of internal elements.";
                    errorFlag = THROW_ERROR;
                }
                elementIndices[chunkElements] = count;
                numElements++;
                chunkElements++;
            } else if (element.dim == numDim-1) {
                if (finalFaceElementType == NoRef) {
                   finalFaceElementType = element.type;
                } else if (finalFaceElementType != element.type) {
                   errorMsg = "Finley can only handle a single type of face elements.";
                   errorFlag = THROW_ERROR;
                }
                faceElementIndices[chunkFaceElements] = count;
                numFaceElements++;
                chunkFaceElements++;
            } else {
                chunkOtherElements++;
            }
    #ifdef ESYS_MPI
            if (count < chunkSize - 1)
                continue;

            // the last chunk is left for the master process
            if (++cpuId >= mpiInfo->size) {
                continue;
            }

            if (errorFlag) {
                for(; cpuId < mpiInfo->size; cpuId++) {
                    MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpiInfo->comm);
                }
                break;
            }
            int chunkInfo[2];
            chunkInfo[0] = chunkElements;
            chunkInfo[1] = chunkFaceElements;

            MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpiInfo->comm);
            MPI_Send(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, cpuId, 81720, mpiInfo->comm);
            MPI_Send(&id[0], chunkSize, MPI_INT, cpuId, 81721, mpiInfo->comm);
            MPI_Send(&tag[0], chunkSize, MPI_INT, cpuId, 81722, mpiInfo->comm);
            MPI_Send(&elementType[0], chunkSize, MPI_INT, cpuId, 81723, mpiInfo->comm);
            MPI_Send(chunkInfo, 2, MPI_INT, cpuId, 81724, mpiInfo->comm);
            MPI_Send(&elementIndices[0], chunkElements, MPI_INT, cpuId, 81725, mpiInfo->comm);
            MPI_Send(&faceElementIndices[0], chunkFaceElements, MPI_INT, cpuId, 81726, mpiInfo->comm);

            // reset arrays for next cpu
    #pragma omp parallel for schedule(static)
            for (index_t i = 0; i < chunkSize*MAX_numNodes_gmsh; i++)
                vertices[i] = -1;
    #pragma omp parallel for schedule(static)
            for (index_t i = 0; i < chunkSize; i++) {
                id[i] = -1;
                tag[i] = -1;
                elementType[i] = NoRef;
            }
            chunkElements = 0;
            chunkFaceElements = 0;
            chunkOtherElements = 0;
            count = -1;
    #endif
        }
    }

#ifdef ESYS_MPI
    if (mpiInfo->size > 1)
        MPI_Bcast(&errorFlag, 1, MPI_INT,  0, mpiInfo->comm);
#endif
    if(errorFlag)
        return errorFlag;

    // all elements have been read and shared, now we have to identify the
    // elements for finley
    if (finalElementType == NoRef) {
        if (numDim == 1) {
           finalElementType = Line2;
        } else if (numDim == 2) {
           finalElementType = Tri3;
        } else if (numDim == 3) {
           finalElementType = Tet4;
        }
    }
    if (finalFaceElementType == NoRef) {
        if (numDim == 1) {
           finalFaceElementType = Point1;
        } else if (numDim == 2) {
           finalFaceElementType = Line2;
        } else if (numDim == 3) {
           finalFaceElementType = Tri3;
        }
    }
    if (finalFaceElementType == Line2) {
        contactElementType = Line2_Contact;
    } else if (finalFaceElementType == Line3 || finalFaceElementType == Line3Macro) {
        contactElementType = Line3_Contact;
    } else if (finalFaceElementType == Tri3) {
        contactElementType = Tri3_Contact;
    } else if (finalFaceElementType == Tri6 || finalFaceElementType == Tri6Macro) {
        contactElementType = Tri6_Contact;
    } else {
        contactElementType = Point1_Contact;
    }

#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpiInfo->size > 1) {
        int msg[3] = {finalElementType, finalFaceElementType,
                contactElementType};
        MPI_Bcast(msg, 3, MPI_INT,  0, mpiInfo->comm);
    }
#endif

    refElements.reset(new ReferenceElementSet(finalElementType, order, reducedOrder));
    refFaceElements.reset(new ReferenceElementSet(finalFaceElementType, order, reducedOrder));
    refContactElements.reset(new ReferenceElementSet(contactElementType, order, reducedOrder));
    refPoints.reset(new ReferenceElementSet(Point1, order, reducedOrder));
    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    dom->setElements(elements);
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    dom->setFaceElements(faces);
    ElementFile* contacts = new ElementFile(refContactElements, mpiInfo);
    dom->setContactElements(contacts);
    ElementFile* points = new ElementFile(refPoints, mpiInfo);
    dom->setPoints(points);

    elements->allocTable(chunkElements);
    faces->allocTable(chunkFaceElements);
    contacts->allocTable(0);
    points->allocTable(0);
    elements->minColor = 0;
    elements->maxColor = chunkElements - 1;
    faces->minColor = 0;
    faces->maxColor = chunkFaceElements - 1;
    contacts->minColor = 0;
    contacts->maxColor = 0;
    points->minColor = 0;
    points->maxColor = 0;

#pragma omp parallel for schedule(static)
    for (index_t e = 0; e < chunkElements; e++) {
        elements->Id[e] = id[elementIndices[e]];
        elements->Tag[e] = tag[elementIndices[e]];
        elements->Color[e] = elementIndices[e];
        elements->Owner[e] = mpiInfo->rank;
        for (int j = 0; j < elements->numNodes; ++j)  {
            int vertex = vertices[INDEX2(j, elementIndices[e], MAX_numNodes_gmsh)];
            elements->Nodes[INDEX2(j, e, elements->numNodes)] = vertex;
        }
    }

#pragma omp parallel for schedule(static)
    for (index_t e = 0; e < chunkFaceElements; e++) {
        faces->Id[e] = id[faceElementIndices[e]];
        faces->Tag[e] = tag[faceElementIndices[e]];
        faces->Color[e] = e;
        faces->Owner[e] = mpiInfo->rank;
        for (int j = 0; j < faces->numNodes; ++j) {
            int faceVertex = vertices[INDEX2(j, faceElementIndices[e], MAX_numNodes_gmsh)];
            faces->Nodes[INDEX2(j, e, faces->numNodes)] = faceVertex;
        }
    }

    return errorFlag;
}

int gather_nodes(FILE* f, std::map<int,int>& tags, std::string& errorMsg,
                 int dim, double version, const std::string& filename, Msh4Entities TagMap)
{
    int numNodes=0;
    int numEntityBlocks;
    std::vector<char> line;
    if (!get_line(line, f))
        return EARLY_EOF;
    int scan_ret;
    if(version >= 4.1){
        int minNodes, maxNodes;
        scan_ret = sscanf(&line[0], "%d %d %d %d", &numEntityBlocks, &numNodes, &minNodes, &maxNodes);
    } else if (version == 4.0){
        int tmp;
        scan_ret = sscanf(&line[0], "%d %d", &tmp, &numNodes);
        numNodes += tmp;
    } else {
        scan_ret = sscanf(&line[0], "%d", &numNodes);
    }

    if (scan_ret == EOF)
        return EARLY_EOF;
    if(version >= 4.1){
        for (int entity = 0; entity < numEntityBlocks; entity++) {
            int entityDim, entityTag, parametric, numNodesInBlock;
            if (!get_line(line, f))
                return EARLY_EOF;
            scan_ret = sscanf(&line[0], "%d %d %d %d", &entityDim, &entityTag, &parametric, &numNodesInBlock);

            if (parametric == 1){
                errorMsg = "eScript does not supprot nodefiles with parametric coordinates.";
                return THROW_ERROR;
            }

            // Tag information
            for (int nodes = 0; nodes < numNodesInBlock; nodes++){
                int tag;
                if (!get_line(line, f))
                    return EARLY_EOF;
                scan_ret = sscanf(&line[0], "%d", &tag);
                if (scan_ret != 1){
                    errorMsg = "malformed meshfile (broken node section)!";
                    return THROW_ERROR;
                }
            }
            // Node coordinate information
            for (int nodes = 0; nodes < numNodesInBlock; nodes++){
                if (!get_line(line, f))
                    return EARLY_EOF;
                float x, y, z;
                scan_ret = sscanf(&line[0], "%f %f %f", &x, &y, &z);
                if (scan_ret != 3){
                    errorMsg = "malformed meshfile (broken node section)!";
                    return THROW_ERROR;
                }
            }
        }
    } else {
        for (int node = 0; node < numNodes; node++) {
            int tmp = 0;
            if (!get_line(line, f))
                return EARLY_EOF;
            scan_ret = sscanf(&line[0], "%d", &tmp);
            if (scan_ret == EOF) {
                return EARLY_EOF;
            } else if (scan_ret != 1) {
                errorMsg = "malformed meshfile (broken node section)!";
                return THROW_ERROR;
            }
            tags[tmp] = -1;
        }
    }
    if (!get_line(line, f))
        return EARLY_EOF;
    if(!(!strncmp(&line[0], "$ENDNOD", 7) || !strncmp(&line[0], "$ENDNOE", 7) || !strncmp(&line[0], "$EndNodes", 9))){
        std::stringstream ss;
        ss << "readGmsh: malformed mesh file. Expected '$EndNodes', got '"
            << &line[0] << "'";
        errorMsg = ss.str();
        return THROW_ERROR;
    }
    if (!get_line(line, f))
        return EARLY_EOF;
    if (strncmp(&line[0], "$ELM", 4) && strncmp(&line[0], "$Elements", 9)) {
        std::stringstream ss;
        ss << "readGmsh: malformed mesh file. Expected '$Elements', got '"
            << &line[0] << "'";
        errorMsg = ss.str();
        return THROW_ERROR;
    }
    int numElements = -1;
    int numBlocks = -1;
    if (!get_line(line, f))
        return EARLY_EOF;
    if (version >= 4.1){
        int minElements, maxElements;
        scan_ret = sscanf(&line[0], "%d %d %d %d\n", &numBlocks, &numElements, &minElements, &maxElements);
    } else if (version == 4.0){
        scan_ret = sscanf(&line[0], "%d %d\n", &numBlocks, &numElements);
    } else {
        scan_ret = sscanf(&line[0], "%d\n", &numElements);
    }
    if (scan_ret == EOF) {
        return EARLY_EOF;
    } else if ((scan_ret != 1 && version < 4.0) || (scan_ret != 2 && version == 4.0) || (scan_ret != 4 && version >= 4.1)) {
        errorMsg = "readGmsh: malformed mesh file ($Elements section contains incorrect header information)";
        return THROW_ERROR;
    }
    struct ElementInfo e;
    std::vector<int> v(MAX_numNodes_gmsh, -1);
    e.vertex = &v[0];

    if(version >= 4.1){
        int entityTag, entityDim, gmsh_type;
        for(int blocks = 0; blocks < numBlocks; blocks++){
            if (!get_line(line, f))
                return EARLY_EOF;
            scan_ret = sscanf(&line[0], "%d %d %d %d", &entityDim, &entityTag, &gmsh_type, &numElements); //Note that entityDim and entityTag are switched in version 4.0
            for(int elementNumber = 0; elementNumber < numElements; elementNumber++){
                if (!get_line(line, f))
                    return EARLY_EOF;
                getSingleElementMSH4(f, dim, version, e, errorMsg, filename, false, gmsh_type, &line[0]);
                
                // if(e.dim == 0)
                //     e.tag = TagMap.pointTags.find(entityTag)->second;
                // else if (e.dim == 1)
                //     e.tag = TagMap.curveTags.find(entityTag)->second;
                // else if (e.dim == 2)
                //     e.tag = TagMap.surfaceTags.find(entityTag)->second;
                // else if (e.dim == 3)
                //     e.tag = TagMap.volumeTags.find(entityTag)->second;
            }
        }
    } else if (version == 4.0){
        for(int x = 0; x < numBlocks; x++){
            int entityTag, entityDim, gmsh_type;
            if (!get_line(line, f))
                return EARLY_EOF;
            scan_ret = sscanf(&line[0], "%d %d %d %d", &entityTag, &entityDim, &gmsh_type, &numElements);
            for(int i = 0; i < numElements; i++){
                if (!get_line(line, f))
                    return EARLY_EOF;
                getSingleElementMSH4(f, dim, version, e, errorMsg, filename, false, gmsh_type, &line[0]);

                // if(e.dim == 0)
                //     e.tag = TagMap.pointTags.find(entityTag)->second;
                // else if (e.dim == 1)
                //     e.tag = TagMap.curveTags.find(entityTag)->second;
                // else if (e.dim == 2)
                //     e.tag = TagMap.surfaceTags.find(entityTag)->second;
                // else if (e.dim == 3)
                //     e.tag = TagMap.volumeTags.find(entityTag)->second;
            }
        }
    } else {
        for (int element = 0; element < numElements; element++) {
            getSingleElement(f, dim, version, e, errorMsg, filename, false);
            for (int i = 0; i < MAX_numNodes_gmsh && v[i] >= 0; i++) {
                std::map<int,int>::iterator it = tags.find(v[i]);
                if (it == tags.end()) {
                    std::stringstream ss;
                    ss << "readGmsh: element contains unknown node (node " << v[i]
                        << ")";
                    errorMsg = ss.str();
                    return THROW_ERROR;
                }
                // the first tagged element using a node tags that node too
                if (it->second == -1 && e.tag != 0)
                    tags[v[i]] = e.tag;
            }
        }
    }
    return 0;
}

int getNodesMaster(escript::JMPI& mpiInfo, FinleyDomain* dom, FILE* fileHandle,
                   int numDim, std::string& errorMsg, std::map<int, int>& tags, double version)
{
    int errorFlag = 0;
    std::vector<char> line;
    if (!get_line(line, fileHandle))
        errorFlag = EARLY_EOF;

    int numNodes = 0, numBlocks = 0;
    int scan_ret;
    if(version >= 4.1){
        int minNodeTag, maxNodeTag;
        scan_ret = sscanf(&line[0], "%d %d %d %d", &numBlocks, &numNodes, &minNodeTag, &maxNodeTag);
    }else if(version == 4.0){
        scan_ret = sscanf(&line[0], "%d %d", &numBlocks, &numNodes);
    } else {
        scan_ret = sscanf(&line[0], "%d", &numNodes);
    }

    SSCANF_CHECK(scan_ret);
#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpiInfo->size > 1)
        MPI_Bcast(&numNodes, 1, MPI_INT,  0, mpiInfo->comm);
#endif
    int chunkSize = numNodes / mpiInfo->size;
    const int rest = numNodes - (mpiInfo->size-1)*chunkSize;
    const size_t storage = std::max(chunkSize, rest);
    std::vector<int> tempInts(storage+1, -1);
    std::vector<double> tempCoords(storage*numDim, -1.);

    int totalNodes = 0;

    for (int nextCPU = mpiInfo->size-1; nextCPU >= 0; nextCPU--) {
        if (nextCPU == 0)
            chunkSize = rest;

        if (!errorFlag) {
            //read in chunksize nodes
            if(version == 4.1){
                int nodeCounter = -1;
                for (int x = 0; x < numBlocks; x++) {
                    std::vector<char> line;
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;

                    if (!strncmp(&line[0], "$ENDNOD", 7) || !strncmp(&line[0], "$ENDNOE", 7) || !strncmp(&line[0], "$EndNodes", 9)) {
                        errorMsg = "readGmsh: found end node string while still reading nodes!";
                        errorFlag = THROW_ERROR;
                        break;
                    }

                    int entityTag, entityDim, parametric, numDataPoints;
                    scan_ret = sscanf(&line[0], "%d %d %d %d\n", &entityDim, &entityTag, &parametric, &numDataPoints);

                    if (parametric == 1){
                        errorMsg = "eScript does not support MSH files with parametric coordinates.";
                        return THROW_ERROR;
                    }

                    // Tag information. At the moment this information is discarded.
                    int temp = nodeCounter;
                    for (int nodes = 0; nodes < numDataPoints; nodes++){
                        if (!get_line(line, fileHandle))
                            return EARLY_EOF;
                        temp++;
                        scan_ret = sscanf(&line[0], "%d", &tempInts[temp]);
                        if (scan_ret != 1){
                            errorMsg = "malformed meshfile (broken node section)!";
                            return THROW_ERROR;
                        }
                    }

                    for(int j = 0; j < numDataPoints; j++){
                        // Get the next line
                        if (!get_line(line, fileHandle))
                            return EARLY_EOF;

                        // Read the information
                        nodeCounter++;
                        if (1 == numDim) {
                            scan_ret = sscanf(&line[0], "%le\n", &tempCoords[0+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (2 == numDim) {
                            scan_ret = sscanf(&line[0], "%le %le\n", &tempCoords[0+nodeCounter*numDim], &tempCoords[1+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (3 == numDim) {
                            scan_ret = sscanf(&line[0], "%le %le %le\n", &tempCoords[0+nodeCounter*numDim], &tempCoords[1+nodeCounter*numDim], &tempCoords[2+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        }
                    }
                    totalNodes += numDataPoints;
                }

                // Possible error (if nodes file is malformed)
                if (totalNodes > numNodes) {
                    std::stringstream ss;
                    ss << "readGmsh: too many nodes (" << totalNodes << " < "
                        << numNodes << ")";
                    errorMsg = ss.str();
                    errorFlag = THROW_ERROR;
                    break;
                }

            } else if(version == 4.0){
                int nodeCounter = -1;
                for (int x = 0; x < numBlocks; x++) {
                    std::vector<char> line;
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;

                    if (!strncmp(&line[0], "$ENDNOD", 7) || !strncmp(&line[0], "$ENDNOE", 7) || !strncmp(&line[0], "$EndNodes", 9)) {
                        errorMsg = "readGmsh: found end node string while still reading nodes!";
                        errorFlag = THROW_ERROR;
                        break;
                    }

                    int entityTag, entityDim, parametric, numNodes;
                    scan_ret = sscanf(&line[0], "%d %d %d %d\n", &entityTag, &entityDim, &parametric, &numNodes);
                    
                    if (parametric == 1){
                        errorMsg = "eScript does not support nodefiles with parametric coordinates.";
                        return THROW_ERROR;
                    }

                    for(int j = 0; j < numNodes; j++){
                        // Get the next line
                        if (!get_line(line, fileHandle))
                            return EARLY_EOF;

                        // Read the information
                        nodeCounter++;
                        if (1 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le\n", &tempInts[nodeCounter], &tempCoords[0+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (2 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le %le\n", &tempInts[nodeCounter], &tempCoords[0+nodeCounter*numDim], &tempCoords[1+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (3 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le %le %le\n", &tempInts[nodeCounter], &tempCoords[0+nodeCounter*numDim], &tempCoords[1+nodeCounter*numDim], &tempCoords[2+nodeCounter*numDim]);
                            SSCANF_CHECK(scan_ret);
                        }
                    }
                    totalNodes += numNodes;
                }

                // Possible error (if nodes file is malformed)
                if (totalNodes > numNodes) {
                    std::stringstream ss;
                    ss << "readGmsh: too many nodes (" << totalNodes << " < "
                        << numNodes << ")";
                    errorMsg = ss.str();
                    errorFlag = THROW_ERROR;
                    break;
                }

            } else { // Not msh version 4.0 or higher
                for (int chunkNodes = 0; chunkNodes < chunkSize; chunkNodes++) {
                    if (totalNodes > numNodes) {
                        std::stringstream ss;
                        ss << "readGmsh: too many nodes (" << totalNodes << " < "
                            << numNodes << ")";
                        errorMsg = ss.str();
                        errorFlag = THROW_ERROR;
                        break;
                    }
                    std::vector<char> line;
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;

                    if (!strncmp(&line[0], "$ENDNOD", 7) || !strncmp(&line[0], "$ENDNOE", 7) || !strncmp(&line[0], "$EndNodes", 9)) {
                        errorMsg = "readGmsh: found end node string while still reading nodes!";
                        errorFlag = THROW_ERROR;
                        break;
                    } else {
                        if (1 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le\n", &tempInts[chunkNodes], &tempCoords[0+chunkNodes*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (2 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le %le\n", &tempInts[chunkNodes], &tempCoords[0+chunkNodes*numDim], &tempCoords[1+chunkNodes*numDim]);
                            SSCANF_CHECK(scan_ret);
                        } else if (3 == numDim) {
                            scan_ret = sscanf(&line[0], "%d %le %le %le\n", &tempInts[chunkNodes], &tempCoords[0+chunkNodes*numDim], &tempCoords[1+chunkNodes*numDim], &tempCoords[2+chunkNodes*numDim]);
                            SSCANF_CHECK(scan_ret);
                        }
                        totalNodes++;
                    }
                }
            }
        }
#ifdef ESYS_MPI
        if (nextCPU != 0) {
            // if there was an error, we have to stop them waiting for more
            MPI_Send(&errorFlag, 1, MPI_INT, nextCPU, 81719, mpiInfo->comm);
            // send out this chunk of mesh to the next waiting node
            if (!errorFlag) {
                // The message has one more int to send chunkSize
                tempInts[chunkSize] = chunkSize;
                MPI_Send(&tempInts[0], chunkSize+1, MPI_INT, nextCPU, 81720, mpiInfo->comm);
                if (chunkSize > 0)
                    MPI_Send(&tempCoords[0], chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpiInfo->comm);
            }
        }
#endif
    }

#ifdef ESYS_MPI
    if (mpiInfo->size > 1)
        MPI_Bcast(&errorFlag, 1, MPI_INT, 0, mpiInfo->comm);
#endif
    if (errorFlag)
        return errorFlag;

    NodeFile* nodes = dom->getNodes();
    nodes->allocTable(chunkSize);

//#pragma omp parallel for schedule(static)
    for(index_t i = 0; i < chunkSize; i++) {
        nodes->Id[i] = tempInts[i];
        nodes->globalDegreesOfFreedom[i] = tempInts[i];
        int tag = tags[tempInts[i]];
        if (tag == -1) {
            nodes->Tag[i] = tempInts[i]; //set tag to node label
        } else {
            nodes->Tag[i] = tag; //set tag of element
        }
        for (int j=0; j<numDim; j++) {
            nodes->Coordinates[INDEX2(j,i,numDim)] = tempCoords[i*numDim+j];
        }
    }

    return errorFlag;
}

int get_next_state(FILE *f, bool nodesRead, bool elementsRead, int *logicFlag) {
    std::vector<char> line;
    do {
        if (!get_line(line, f)) { //got no line
            //check to see we at least have some nodes and elements
            if (!nodesRead) {
                //EOF before nodes section found
                return MISSING_NODES;
            } else if(!elementsRead){
                //EOF before elements section found
                return MISSING_ELEMENTS;
            }
            *logicFlag = 4;
            return SUCCESS; //EOF as expected
        }
//        if (line[0] != '$')
//            fprintf(stderr, "consuming line: %s", &line[0]);
    } while(line[0] != '$');

    if (!strncmp(&line[1], "MeshFormat", 10)) {
        *logicFlag = 1;
    } else if (!strncmp(&line[0], "$NOD", 4) || !strncmp(&line[0], "$NOE", 4) || !strncmp(&line[0], "$Nodes", 6)) {
        *logicFlag = 2;
    } else if (!strncmp(&line[1], "ELM", 3) || !strncmp(&line[1], "Elements", 8)) {
        *logicFlag = 3;
    } else if (!strncmp(&line[1], "PhysicalNames", 13)) {
        *logicFlag = 4;
    } else if (!strncmp(&line[1], "Entities", 8)) {
        *logicFlag = 5;
    } else {
        *logicFlag = 0;
    }
    return 0;
}

void send_state(escript::JMPI& mpiInfo, int error, int logic)
{
#ifdef ESYS_MPI
    int flags[2] = {error, logic};
    // Broadcast line
    if (mpiInfo->size > 1) {
        MPI_Bcast(&flags, 2, MPI_INT,  0, mpiInfo->comm);
    }
#endif
}

int check_error(int error, FILE* f, const std::string& errorMsg)
{
    //handle errors
    switch(error) {
        case 0:
            break;
        case ERROR:
            throw FinleyException("ERROR set for unknown reason");
        case EARLY_EOF: //early EOF while scanning
            throw IOError("early EOF while scanning");
        case MISSING_NODES:  //EOF before nodes section found
            throw IOError("EOF before nodes section found");
        case MISSING_ELEMENTS:
            throw IOError("EOF before elements section found");
        case THROW_ERROR: // throw errorMsg
            throw IOError(errorMsg);
        case SUCCESS: // EOF at apropriate time.
            if (f)
                fclose(f);
            break;
        default:
            throw FinleyException("an unknown error has occured in readGmsh");

    }
    return error;
}

FinleyDomain* readGmshMaster(escript::JMPI& mpiInfo,
                             const std::string& filename, int numDim,
                             int order, int reducedOrder, bool optimize,
                             bool useMacroElements)
{
    double version = 1.0;
    bool nodesRead = false, elementsRead = false;
    int numPhysicalNames = 0;
    Msh4Entities TagMap;
    int format = 0, size = sizeof(double), scan_ret, errorFlag = 0, logicFlag = 0;
    std::vector<char> line;
    std::map<int,int> nodeTags;
    std::string errorMsg;

    size_t found = filename.find("\n");
    if (found != std::string::npos) {
        errorFlag = THROW_ERROR;
        send_state(mpiInfo, errorFlag, logicFlag);
        throw escript::ValueError("readGmsh: filename contains newline characters!");
    }

    // allocate mesh
    FinleyDomain* dom = new FinleyDomain(filename, numDim, mpiInfo);

    // get file handle
#ifdef _WIN32
    FILE* fileHandle = fopen(filename.c_str(), "rb"); // open in binary mode to allow seek
#else
    FILE* fileHandle = fopen(filename.c_str(), "r");
#endif
    if (!fileHandle) {
        std::stringstream ss;
        ss << "readGmsh: opening file " << filename << " for reading failed.";
        errorMsg = ss.str();
        errorFlag = THROW_ERROR;
        send_state(mpiInfo, errorFlag, logicFlag);
        throw IOError(errorMsg);
    }
    // start reading
    while (!errorFlag) {
        // find line starting with $
        logicFlag = 0;
        while(logicFlag == 0){
            errorFlag = get_next_state(fileHandle, nodesRead, elementsRead, &logicFlag);
        }
        send_state(mpiInfo, errorFlag, logicFlag);
        //pre-logic error check
        if (check_error(errorFlag, fileHandle, errorMsg) == SUCCESS)
            break;
        // format
        if (logicFlag == 1 && !errorFlag) {
            std::vector<char> fmt;
            if (!get_line(fmt, fileHandle))
                errorFlag = EARLY_EOF;
            scan_ret = sscanf(&fmt[0], "%lf %d %d\n", &version, &format, &size);
            SSCANF_CHECK(scan_ret);

            if(version != 2.0 && version != 2.1 && version != 2.2 && version != 4.0 && version != 4.1)
                throw FinleyException("Cannot understand this msh format version. Please use version 2.2 or 4.1");
            
        }
        // nodes are read
        else if (logicFlag == 2 && !errorFlag) {
            nodesRead = true;
            std::vector<int> sendable_map;
            long current = ftell(fileHandle);
            errorFlag = gather_nodes(fileHandle, nodeTags, errorMsg,
                    numDim, version, filename.c_str(), TagMap);
            if (!errorFlag && fseek(fileHandle, current, SEEK_SET) < 0) {
                errorMsg = "Error in file operation";
                errorFlag = THROW_ERROR;
            }
            send_state(mpiInfo, errorFlag, logicFlag);
            check_error(errorFlag, fileHandle, errorMsg);
#ifdef ESYS_MPI
            int mapsize = 2 * nodeTags.size();
            sendable_map.resize(mapsize);
            std::map<int,int>::iterator i = nodeTags.begin();
            for (int j = 0; i != nodeTags.end(); i++, j += 2) {
                sendable_map[j] = i->first;
                sendable_map[j + 1] = i->second;
            }
            if (mpiInfo->size > 1) {
                MPI_Bcast(&mapsize, 1, MPI_INT, 0, mpiInfo->comm);
                sendable_map.resize(mapsize);
                MPI_Bcast(&sendable_map[0], mapsize, MPI_INT, 0, mpiInfo->comm);
                for (int j = 0; j < mapsize; j += 2)
                    nodeTags[sendable_map[j]] = sendable_map[j + 1];
            }
#endif
            errorFlag = getNodesMaster(mpiInfo, dom, fileHandle, numDim,
                                       errorMsg, nodeTags, version);
        }
        // elements
        else if (logicFlag == 3 && !errorFlag) {
            elementsRead = true;
            errorFlag = getElementsMaster(mpiInfo, dom, fileHandle, errorMsg,
                            useMacroElements, filename, numDim, version, order,
                            reducedOrder, TagMap);
        }
        // name tags
        // (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca)
        else if (logicFlag == 4 && !errorFlag) {
            std::vector<char> names;
            if (!get_line(names, fileHandle))
                errorFlag = EARLY_EOF;
            scan_ret = sscanf(&names[0], "%d", &numPhysicalNames);
            SSCANF_CHECK(scan_ret);
#ifdef ESYS_MPI
            // Broadcast numNames if there are multiple mpi procs
            if (mpiInfo->size > 1)
                MPI_Bcast(&numPhysicalNames, 1, MPI_INT,  0, mpiInfo->comm);
#endif
            for (int i = 0; i < numPhysicalNames; i++) {
                std::vector<char> line;
                char name[1024] = {0};
                if (!get_line(line, fileHandle))
                    errorFlag = EARLY_EOF;
                int tag_info[2] = { 0 };
                char* position = &line[0];
                //skip the first int, it's the physical dimension
                if (next_space(&position, 1) == NULL
                        || sscanf(position, "%d", tag_info) != 1
                        || next_space(&position, 1) == NULL
                        || sscanf(position, "%s", name) != 1) {
                    errorFlag = ERROR;
                }
                name[strlen(name)-1]='\0'; //strip trailing "

#ifdef ESYS_MPI
                // broadcast the tag info
                if (mpiInfo->size > 1) {
                    tag_info[1] = strlen(name) + 1; //include \0
                    MPI_Bcast(tag_info, 2, MPI_INT,  0, mpiInfo->comm);
                    MPI_Bcast(&name, tag_info[1], MPI_CHAR,  0, mpiInfo->comm);
                }
#endif
                dom->setTagMap(name+1, tag_info[0]); //skip leading "
            }
        } else if (logicFlag == 5 && !errorFlag) {
            // If necessary, read in physical tag information from the Entities section
            if(version >= 4.0)
            {
                int numPoints, numCurves, numSurfaces, numVolumes;
                if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;
                scan_ret = sscanf(&line[0], "%d %d %d %d\n", &numPoints, &numCurves, &numSurfaces, &numVolumes);

                // Skip over the curve and surface information
                for(int i = 0; i < numPoints; i++)
                {
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;
                    int pointTag, numPhysicalTags, physicalTag;
                    if(version == 4.0)
                    {
                        float tmp[6]={-1};
                        sscanf(&line[0], "%d %f %f %f %f %f %f %d %d", &pointTag, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &numPhysicalTags, &physicalTag);
                    }
                    else
                    {
                        float tmp[3]={-1};
                        sscanf(&line[0], "%d %f %f %f %d %d", &pointTag, &tmp[0], &tmp[1], &tmp[2], &numPhysicalTags, &physicalTag);
                    }
                    if(numPhysicalTags > 0)
                        TagMap.pointTags.insert(std::pair<int,int>(pointTag,physicalTag));
                }

                for(int i = 0; i < numCurves; i++)
                {
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;
                    int pointTag, numPhysicalTags, physicalTag;
                    float tmp[6]={-1};
                    scan_ret = sscanf(&line[0], "%d %f %f %f %f %f %f %d %d", &pointTag, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &numPhysicalTags, &physicalTag);
                    if(numPhysicalTags > 0)
                        TagMap.curveTags.insert(std::pair<int,int>(pointTag,physicalTag));
                }

                for(int i = 0; i < numSurfaces; i++)
                {
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;
                    signed int pointTag, numPhysicalTags, physicalTag;
                    float tmp[6]={-1};
                    scan_ret = sscanf(&line[0], "%d %f %f %f %f %f %f %d %d", &pointTag, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &numPhysicalTags, &physicalTag);
                    if(numPhysicalTags > 0)
                        TagMap.surfaceTags.insert(std::pair<int,int>(pointTag,physicalTag));
                }

                for(int i = 0; i < numVolumes; i++)
                {
                    if (!get_line(line, fileHandle))
                        errorFlag = EARLY_EOF;
                    int pointTag, numPhysicalTags, physicalTag;
                    float tmp[6]={-1};
                    scan_ret = sscanf(&line[0], "%d %f %f %f %f %f %f %d %d", &pointTag, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &numPhysicalTags, &physicalTag);
                    if(numPhysicalTags > 0)
                        TagMap.volumeTags.insert(std::pair<int,int>(pointTag,physicalTag));
                }
            }
        }

        if (!get_line(line, fileHandle)) {
            errorFlag = EARLY_EOF;
        }
        if (line[0] != '$') {
            errorFlag = THROW_ERROR;
            std::stringstream ss;
            ss << "readGmsh: expected closing tag, got '"
                << &line[0] << "'...";
            errorMsg = ss.str();
        }
        send_state(mpiInfo, errorFlag, logicFlag);
        //post logic error check, throws if relevant
        check_error(errorFlag, fileHandle, errorMsg);
    }
    return dom;
}

// slave-only functions follow

#ifdef ESYS_MPI
int getNodesSlave(escript::JMPI& mpiInfo, FinleyDomain* dom, int numDim,
                  std::string& errorMsg, std::map<int, int>& tags)
{
    if (mpiInfo->size == 1)
        throw FinleyException("Slave function called without master!");

    int errorFlag = 0;
    int numNodes = 0;

    // get numNodes from the master
    MPI_Bcast(&numNodes, 1, MPI_INT,  0, mpiInfo->comm);
    int chunkSize = numNodes / mpiInfo->size, chunkNodes = 0;
    std::vector<int> tempInts(chunkSize+1); // Stores the integer message data
    std::vector<double> tempCoords(chunkSize*numDim); // Stores the double message data
    // Each worker receives two messages
    MPI_Status status;
    MPI_Recv(&errorFlag, 1, MPI_INT, 0, 81719, mpiInfo->comm, &status);
    if (!errorFlag) {
        MPI_Recv(&tempInts[0], chunkSize+1, MPI_INT, 0, 81720, mpiInfo->comm, &status);
        if (chunkSize > 0)
            MPI_Recv(&tempCoords[0], chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpiInfo->comm, &status);
        chunkNodes = tempInts[chunkSize]; // How many nodes are in this worker's chunk?
    }

    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, mpiInfo->comm);
    if (errorFlag)
        return errorFlag;

    NodeFile* nodes = dom->getNodes();
    nodes->allocTable(chunkNodes);

//#pragma omp parallel for schedule(static)
    for (index_t i = 0; i < chunkNodes; i++) {
        nodes->Id[i] = tempInts[i];
        nodes->globalDegreesOfFreedom[i] = tempInts[i];
        int tag = tags[tempInts[i]];
        if (tag == -1) {
            nodes->Tag[i] = tempInts[i]; //set tag to node label
        } else {
            nodes->Tag[i] = tag; //set tag of element
        }
        for (int j = 0; j < numDim; j++) {
            nodes->Coordinates[INDEX2(j,i,numDim)] = tempCoords[i*numDim+j];
        }
    }

    return errorFlag;
}

int getElementsSlave(escript::JMPI& mpiInfo, FinleyDomain* dom,
                     std::string& errorMsg, bool useMacroElements,
                     int numDim, double version, int order, int reducedOrder)
{
    /*
     *  This function should read in the elements and distribute
     *  them to the apropriate process.
     */
    if (mpiInfo->size == 1) {
        errorMsg = "Slave function called with no master";
        return THROW_ERROR; //again, sillyness
    }

    int errorFlag = 0;

    ElementTypeId finalElementType = NoRef;
    ElementTypeId finalFaceElementType = NoRef;
    ElementTypeId contactElementType = NoRef;
    const_ReferenceElementSet_ptr refPoints, refContactElements;
    const_ReferenceElementSet_ptr refFaceElements, refElements;

    int totalNumElements = 0;
    MPI_Bcast(&totalNumElements, 1, MPI_INT, 0, mpiInfo->comm);

    int chunkSize = totalNumElements / mpiInfo->size, chunkElements = 0;
    int chunkFaceElements = 0;
    std::vector<int> id(chunkSize+1);
    std::vector<int> tag(chunkSize+1);
    std::vector<int> vertices(chunkSize*MAX_numNodes_gmsh, -1);
    std::vector<ElementTypeId> elementType(chunkSize+1);
    std::vector<int> elementIndices(chunkSize, -1);
    std::vector<int> faceElementIndices (chunkSize, -1);

    //chunkInfo stores the number of elements and number of face elements
    int chunkInfo[2];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < chunkSize; i++) {
        id[i] = -1;
        tag[i] = -1;
        elementType[i] = NoRef;
    }

    // Each worker receives messages
    MPI_Status status;

    MPI_Recv(&errorFlag, 1, MPI_INT, 0, 81719, mpiInfo->comm, &status);
    if (errorFlag)
        return errorFlag;

    MPI_Recv(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, 0, 81720, mpiInfo->comm, &status);
    MPI_Recv(&id[0], chunkSize, MPI_INT, 0, 81721, mpiInfo->comm, &status);
    MPI_Recv(&tag[0], chunkSize, MPI_INT, 0, 81722, mpiInfo->comm, &status);
    MPI_Recv(&elementType[0], chunkSize, MPI_INT, 0, 81723, mpiInfo->comm, &status);
    MPI_Recv(chunkInfo, 2, MPI_INT, 0, 81724, mpiInfo->comm, &status);
    chunkElements = chunkInfo[0];
    chunkFaceElements = chunkInfo[1];
    MPI_Recv(&elementIndices[0], chunkElements, MPI_INT, 0, 81725, mpiInfo->comm,&status);
    MPI_Recv(&faceElementIndices[0], chunkFaceElements, MPI_INT, 0, 81726, mpiInfo->comm,&status);

    MPI_Bcast(&errorFlag, 1, MPI_INT,  0, mpiInfo->comm);
    if (errorFlag)
        return errorFlag;

    // all elements have been read and shared, now we have to identify the
    // elements for finley
    int numNodes[3] = {0,0,0};
    MPI_Bcast(numNodes, 3, MPI_INT,  0, mpiInfo->comm);
    finalElementType = static_cast<ElementTypeId>(numNodes[0]);
    finalFaceElementType = static_cast<ElementTypeId>(numNodes[1]);
    contactElementType = static_cast<ElementTypeId>(numNodes[2]);

    refElements.reset(new ReferenceElementSet(finalElementType, order, reducedOrder));
    refFaceElements.reset(new ReferenceElementSet(finalFaceElementType, order, reducedOrder));
    refContactElements.reset(new ReferenceElementSet(contactElementType, order, reducedOrder));
    refPoints.reset(new ReferenceElementSet(Point1, order, reducedOrder));
    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    dom->setElements(elements);
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    dom->setFaceElements(faces);
    ElementFile* contacts = new ElementFile(refContactElements, mpiInfo);
    dom->setContactElements(contacts);
    ElementFile* points = new ElementFile(refPoints, mpiInfo);
    dom->setPoints(points);
    elements->allocTable(chunkElements);
    faces->allocTable(chunkFaceElements);
    contacts->allocTable(0);
    points->allocTable(0);
    elements->minColor = 0;
    elements->maxColor = chunkElements-1;
    faces->minColor = 0;
    faces->maxColor = chunkFaceElements-1;
    contacts->minColor = 0;
    contacts->maxColor = 0;
    points->minColor = 0;
    points->maxColor = 0;

#pragma omp parallel for schedule(static)
    for (index_t e = 0; e < chunkElements; e++) {
        elements->Id[e] = id[elementIndices[e]];
        elements->Tag[e] = tag[elementIndices[e]];
        elements->Color[e] = elementIndices[e];
        elements->Owner[e] = mpiInfo->rank;
        for (int j = 0; j < elements->numNodes; ++j)  {
            int vertex = vertices[INDEX2(j, elementIndices[e], MAX_numNodes_gmsh)];
            elements->Nodes[INDEX2(j, e, elements->numNodes)] = vertex;
        }
    }

#pragma omp parallel for schedule(static)
    for (index_t e = 0; e < chunkFaceElements; e++) {
        faces->Id[e] = id[faceElementIndices[e]];
        faces->Tag[e] = tag[faceElementIndices[e]];
        faces->Color[e] = e;
        faces->Owner[e] = mpiInfo->rank;
        for (int j = 0; j < faces->numNodes; ++j) {
            int faceVertex = vertices[INDEX2(j, faceElementIndices[e], MAX_numNodes_gmsh)];
            faces->Nodes[INDEX2(j, e, faces->numNodes)] = faceVertex;
        }
    }

    return errorFlag;
}

void recv_state(escript::JMPI& mpiInfo, int* error, int* logic)
{
    int flags[2] = { 0 };
    // Broadcast line
    MPI_Bcast(&flags, 2, MPI_INT, 0, mpiInfo->comm);
    *error = flags[0];
    if (logic)
        *logic = flags[1];
}
#endif // ESYS_MPI

FinleyDomain* readGmshSlave(escript::JMPI& mpiInfo,
                            const std::string& filename, int numDim, int order,
                            int reducedOrder, bool optimize,
                            bool useMacroElements)
{
#ifndef ESYS_MPI
    throw FinleyException("slave function called in non-MPI build!");
#else
    if (mpiInfo->size == 1)
        throw FinleyException("slave function called but only one process");

    const double version = 1.0;
    int errorFlag = 0, logicFlag = 0;
    std::string errorMsg;

    // allocate mesh
    FinleyDomain* dom = new FinleyDomain(filename, numDim, mpiInfo);

    while (errorFlag != SUCCESS) {
        logicFlag = 0;
        // pre logic state fetch
        recv_state(mpiInfo, &errorFlag, &logicFlag);
        if (check_error(errorFlag, NULL, errorMsg) == SUCCESS)
            break;

        // format
        // nodes are read
        if (logicFlag == 2) {
            int mapsize = 0;
            std::vector<int> sendable_map;
            recv_state(mpiInfo, &errorFlag, &logicFlag);
            check_error(errorFlag, NULL, errorMsg);
            MPI_Bcast(&mapsize, 1, MPI_INT, 0, mpiInfo->comm);
            sendable_map.resize(mapsize);
            MPI_Bcast(&sendable_map[0], mapsize, MPI_INT, 0, mpiInfo->comm);
            std::map<int,int> nodeTags;
            for (int j = 0; j < mapsize; j += 2)
                nodeTags[sendable_map[j]] = sendable_map[j + 1];

            errorFlag = getNodesSlave(mpiInfo, dom, numDim, errorMsg, nodeTags);
        }
        // elements
        else if (logicFlag == 3) {
            errorFlag = getElementsSlave(mpiInfo, dom, errorMsg,
                                         useMacroElements, numDim, version,
                                         order, reducedOrder);
        }
        // name tags
        // (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca)
        else if (logicFlag == 4) {
            // Broadcast numNames
            int numNames = 0;
            MPI_Bcast(&numNames, 1, MPI_INT,  0, mpiInfo->comm);
            for (int i = 0; i < numNames; i++) {
                char name[1024];
                int tagInfo[2];
                // broadcast the tag info
                tagInfo[0] = 0;
                tagInfo[1] = 0;
                MPI_Bcast(tagInfo, 2, MPI_INT,  0, mpiInfo->comm);
                //strlen + 1 for null terminator
                MPI_Bcast(&name, tagInfo[1], MPI_CHAR, 0, mpiInfo->comm);
                dom->setTagMap(&name[1], tagInfo[0]);
            }
        }
        //post logic error check
        recv_state(mpiInfo, &errorFlag, &logicFlag);
        if (check_error(errorFlag, NULL, errorMsg) == SUCCESS)
            break;
    } //end while loop

    return dom;
#endif // ESYS_MPI
}

} // anonymous namespace


namespace finley {

escript::Domain_ptr FinleyDomain::readGmsh(escript::JMPI mpiInfo,
                        const std::string& filename, int numDim, int order,
                        int reducedOrder, bool optimize, bool useMacroElements)
{
    FinleyDomain* dom;

    if (mpiInfo->rank == 0) {
        dom = readGmshMaster(mpiInfo, filename, numDim, order, reducedOrder,
                             optimize, useMacroElements);
    } else {
        dom = readGmshSlave(mpiInfo, filename, numDim, order, reducedOrder,
                            optimize, useMacroElements);
    }

    // resolve id's
    dom->resolveNodeIds();
    // rearrange elements
    dom->prepare(optimize);
    return dom->getPtr();
}

} // namespace finley
