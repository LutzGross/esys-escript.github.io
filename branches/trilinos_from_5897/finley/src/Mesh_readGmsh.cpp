
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

  Finley: read mesh from gmsh file

*****************************************************************************/

#include "Mesh.h"
#include "FinleyException.h"

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

struct ElementInfo {
    finley::ElementTypeId type;
    int id;
    int dim;
    int *vertex;
    int tag;
};

inline bool is_node_string(const char *line)
{
    if (!line)
        return false;
    return !strncmp(line, "$NOD", 4) || !strncmp(line, "$NOE", 4)
            || !strncmp(line, "$Nodes", 6);
}

inline bool is_endnode_string(const char *line)
{
    if (!line)
        return false;
    return !strncmp(line, "$ENDNOD", 7) || !strncmp(line, "$ENDNOE", 7)
            || !strncmp(line, "$EndNodes", 9);
}

inline bool get_line(std::vector<char>& line, FILE *file)
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
        if (capacity > 1024) {//madness
            res = NULL;
            break;
        }
    }
    return res == tmp; //true if line read, false if EOF without or without \n
}

inline char *next_space(char **position, int count)
{
    for (int i = 0; i < count; i++) {
        *position = strchr(*position, ' ');
        if ((*position)++ == NULL) //move off the space
            return NULL;
    }
    return *position;
}

namespace finley {

int getSingleElement(FILE *f, int dim, double version, struct ElementInfo& e,
        std::string& errorMsg, const char *fname, bool useMacroElements)
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
                    << gmsh_type << " in mesh file " << fname;
                errorMsg = ss.str();
                return THROW_ERROR;
            }
    }
    if (version <= 1.0){
        int tmp = 0;
        if (sscanf(position, "%d %*d %d", &e.tag, &tmp) == 0
                || next_space(&position, 3) == NULL )
            return EARLY_EOF;
        if (tmp != numNodesPerElement) {
            std::stringstream ss;
            ss << "readGmsh: Illegal number of nodes for element " << e.id
                << " in mesh file " << fname;
            errorMsg = ss.str();
            return THROW_ERROR;
        }
    } else {
        e.tag = 1;
        int numTags=0; //this is garbage and never used
        if (sscanf(position, "%d", &numTags) == 0 
                || next_space(&position, 1) == NULL)
            return EARLY_EOF;
        if (sscanf(position, "%d", &e.tag) == 0 
                || next_space(&position, numTags) == NULL)
            return EARLY_EOF;
        /* ignore any other tags, second tag would be elementary id,
         third tag would be partition id */
    }

    for (int j = 0; j < numNodesPerElement; j++) {
        if (sscanf(position, "%d", e.vertex+j) == 0 
                || next_space(&position, 1) == NULL)
            return EARLY_EOF;
    }
    return 0;
}

int getElementsMaster(escript::JMPI& mpi_info, Mesh* mesh_p, FILE* fileHandle_p,
        std::string& errorMsg, bool useMacroElements, const std::string fname,
        int numDim, double version, int order, int reduced_order) {
    /*
     *  This function should read in the elements and distribute
     *  them to the apropriate process.
     */
    int errorFlag=0;
    ElementTypeId final_element_type = NoRef;
    ElementTypeId final_face_element_type = NoRef;
    ElementTypeId contact_element_type = NoRef;
    int numElements=0, numFaceElements=0, totalNumElements=0;
    const_ReferenceElementSet_ptr refPoints, refContactElements;
    const_ReferenceElementSet_ptr refFaceElements, refElements;
    int *id, *tag;
    ElementTypeId * element_type;
    std::vector<char> line;
    if (!get_line(line, fileHandle_p))
        errorFlag = EARLY_EOF;
    int scan_ret = sscanf(&line[0], "%d", &totalNumElements);
    SSCANF_CHECK(scan_ret);

#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpi_info->size > 1) {
        int msg = totalNumElements;
        MPI_Bcast(&msg, 1, MPI_INT, 0, mpi_info->comm);
    }
#endif

    int chunkSize = totalNumElements / mpi_info->size + 1, chunkElements=0;
    int chunkFaceElements=0, chunkOtherElements=0;
    id = new int[chunkSize+1];
    tag = new int[chunkSize+1];
    std::vector<int>vertices(chunkSize*MAX_numNodes_gmsh, -1);
    element_type = new ElementTypeId[chunkSize+1];
    std::vector<int> elementIndices (chunkSize, -1);
    std::vector<int> faceElementIndices (chunkSize, -1);

#ifdef ESYS_MPI
    //chunkInfo stores the number of elements and number of face elements
    int chunkInfo[2];
    int cpuId = 0;
#endif

#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkSize; i++) {
        id[i] = -1;
        tag[i] = -1;
        element_type[i] = NoRef;
    }

    // read all in
    for(int e = 0, count = 0; e < totalNumElements; e++, count++) {
        struct ElementInfo element = {NoRef, 0, 0, &vertices[count*MAX_numNodes_gmsh], 0};
        getSingleElement(fileHandle_p, numDim, version,
                element, errorMsg, fname.c_str(), useMacroElements);
        element_type[count] = element.type;
        id[count] = element.id;
        tag[count] = element.tag;

        /* for tet10 the last two nodes need to be swapped */
        if ((element.type==Tet10) || (element.type == Tet10Macro)) {
            int vertex = vertices[INDEX2(9,count,MAX_numNodes_gmsh)];
            vertices[INDEX2(9,count,MAX_numNodes_gmsh)] = vertices[INDEX2(8,count,MAX_numNodes_gmsh)];
            vertices[INDEX2(8,count,MAX_numNodes_gmsh)] = vertex;
        }

        if (element.dim == numDim) {
            if (final_element_type == NoRef) {
               final_element_type = element.type;
            } else if (final_element_type != element.type) {
                errorMsg = "Finley can only handle a single type of internal elements.";
                errorFlag = THROW_ERROR;
            }
            elementIndices[chunkElements]=count;
            numElements++;
            chunkElements++;
        } else if (element.dim == numDim-1) {
            if (final_face_element_type == NoRef) {
               final_face_element_type = element.type;
            } else if (final_face_element_type != element.type) {
               errorMsg = "Finley can only handle a single type of face elements.";
               errorFlag = THROW_ERROR;
            }
            faceElementIndices[chunkFaceElements]=count;
            numFaceElements++;
            chunkFaceElements++;
        } else {
            chunkOtherElements++;
        }
#ifdef ESYS_MPI
        if(count < chunkSize - 1)
            continue;
        chunkInfo[0]=chunkElements;
        chunkInfo[1]=chunkFaceElements;

        if(cpuId++ > mpi_info->size) {
            continue;
        }
        if(errorFlag){
            for(; cpuId<mpi_info->size; cpuId++) {
                MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
            }
            break;
        }
        MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
        MPI_Send(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, cpuId, 81720, mpi_info->comm);
        MPI_Send(id, chunkSize, MPI_INT, cpuId, 81721, mpi_info->comm);
        MPI_Send(tag, chunkSize, MPI_INT, cpuId, 81722, mpi_info->comm);
        MPI_Send(element_type, chunkSize, MPI_INT, cpuId, 81723, mpi_info->comm);
        MPI_Send(chunkInfo, 2, MPI_INT, cpuId, 81724, mpi_info->comm);
        MPI_Send(&(elementIndices[0]), chunkElements, MPI_INT, cpuId, 81725, mpi_info->comm);
        MPI_Send(&(faceElementIndices[0]), chunkFaceElements, MPI_INT, cpuId, 81726, mpi_info->comm);

        // reset arrays for next cpu
#pragma omp parallel for schedule(static)
        for (int i=0; i<chunkSize*MAX_numNodes_gmsh; i++)
            vertices[i] = -1;
#pragma omp parallel for schedule(static)
        for (int i=0; i<chunkSize; i++) {
            id[i] = -1;
            tag[i] = -1;
            element_type[i] = NoRef;
        }
        chunkElements=0;
        chunkFaceElements=0;
        chunkOtherElements=0;
        count = -1;
#endif
    }

#ifdef ESYS_MPI
    if(mpi_info->size>1){
        MPI_Bcast(&errorFlag, 1, MPI_INT,  0, mpi_info->comm);
    }
#endif
    if(errorFlag){
        return errorFlag;
    }

    // all elements have been read and shared, now we have to identify the
    // elements for finley
    if (final_element_type == NoRef) {
        if (numDim==1) {
           final_element_type=Line2;
        } else if (numDim==2) {
           final_element_type=Tri3;
        } else if (numDim==3) {
           final_element_type=Tet4;
        }
    }
    if (final_face_element_type == NoRef) {
        if (numDim==1) {
           final_face_element_type=Point1;
        } else if (numDim==2) {
           final_face_element_type=Line2;
        } else if (numDim==3) {
           final_face_element_type=Tri3;
        }
    }
    if (final_face_element_type == Line2) {
        contact_element_type=Line2_Contact;
    } else  if ( (final_face_element_type == Line3) || (final_face_element_type == Line3Macro) ) {
        contact_element_type=Line3_Contact;
    } else  if (final_face_element_type == Tri3) {
        contact_element_type=Tri3_Contact;
    } else  if ( (final_face_element_type == Tri6) || (final_face_element_type == Tri6Macro)) {
        contact_element_type=Tri6_Contact;
    } else {
        contact_element_type=Point1_Contact;
    }

#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpi_info->size > 1) {
        int msg[3] = {final_element_type, final_face_element_type, 
                contact_element_type};
        MPI_Bcast(msg, 3, MPI_INT,  0, mpi_info->comm);
    }
#endif



    refElements.reset(new ReferenceElementSet(final_element_type, order, reduced_order));
    refFaceElements.reset(new ReferenceElementSet(final_face_element_type, order, reduced_order));
    refContactElements.reset(new ReferenceElementSet(contact_element_type, order, reduced_order));
    refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));
    mesh_p->Elements=new ElementFile(refElements, mpi_info);
    mesh_p->FaceElements=new ElementFile(refFaceElements, mpi_info);
    mesh_p->ContactElements=new ElementFile(refContactElements, mpi_info);
    mesh_p->Points=new ElementFile(refPoints, mpi_info);
    mesh_p->Elements->allocTable(chunkElements);
    mesh_p->FaceElements->allocTable(chunkFaceElements);
    mesh_p->ContactElements->allocTable(0);
    mesh_p->Points->allocTable(0);
    mesh_p->Elements->minColor=0;
    mesh_p->Elements->maxColor=chunkElements-1;
    mesh_p->FaceElements->minColor=0;
    mesh_p->FaceElements->maxColor=chunkFaceElements-1;
    mesh_p->ContactElements->minColor=0;
    mesh_p->ContactElements->maxColor=0;
    mesh_p->Points->minColor=0;
    mesh_p->Points->maxColor=0;
#pragma omp parallel for schedule(static)
    for(int e = 0; e < chunkElements; e++) {
        mesh_p->Elements->Id[e]=id[elementIndices[e]];
        mesh_p->Elements->Tag[e]=tag[elementIndices[e]];
        mesh_p->Elements->Color[e]=elementIndices[e];
        mesh_p->Elements->Owner[e]=mpi_info->rank;
        for (int j = 0; j < mesh_p->Elements->numNodes; ++j)  {
            int vertex = vertices[INDEX2(j,elementIndices[e],MAX_numNodes_gmsh)];
            mesh_p->Elements->Nodes[INDEX2(j, e, mesh_p->Elements->numNodes)]=vertex;
        }
    }

#pragma omp parallel for schedule(static)
    for (int e = 0; e < chunkFaceElements; e++) {
        mesh_p->FaceElements->Id[e]=id[faceElementIndices[e]];
        mesh_p->FaceElements->Tag[e]=tag[faceElementIndices[e]];
        mesh_p->FaceElements->Color[e]=e;
        mesh_p->FaceElements->Owner[e]=mpi_info->rank;
        for (int j=0; j<mesh_p->FaceElements->numNodes; ++j) {
            int faceVertex = vertices[INDEX2(j,faceElementIndices[e],MAX_numNodes_gmsh)];
            mesh_p->FaceElements->Nodes[INDEX2(j, e, mesh_p->FaceElements->numNodes)]= faceVertex;
        }
    }

    /* and clean up */
    delete[] id;
    delete[] tag;
    delete[] element_type;
    return errorFlag;
}

int getElementsSlave(escript::JMPI& mpi_info, Mesh *mesh_p, FILE *fileHandle_p,
        std::string& errorMsg, bool useMacroElements, const std::string fname,
        int numDim, double version, int order, int reduced_order) {
    /*
     *  This function should read in the elements and distribute
     *  them to the apropriate process.
     */
#ifndef ESYS_MPI
    errorMsg = "Slave function called in non-MPI build";
    return THROW_ERROR; // calling the slave from a non-mpi process is an awful idea
#else
    if (mpi_info->size == 1) {
        errorMsg = "Slave function called with no master";
        return THROW_ERROR; //again, sillyness
    }

    int errorFlag=0;

    ElementTypeId final_element_type = NoRef;
    ElementTypeId final_face_element_type = NoRef;
    ElementTypeId contact_element_type = NoRef;
    int totalNumElements=0;
    const_ReferenceElementSet_ptr refPoints, refContactElements;
    const_ReferenceElementSet_ptr refFaceElements, refElements;
    ElementTypeId * element_type;

    int msg = 0;
    MPI_Bcast(&msg, 1, MPI_INT, 0, mpi_info->comm);
    totalNumElements = msg;

    int chunkSize = totalNumElements / mpi_info->size + 1, chunkElements=0;
    int chunkFaceElements=0;
    int *id = new int[chunkSize+1];
    int *tag = new int[chunkSize+1];
    std::vector<int>vertices(chunkSize*MAX_numNodes_gmsh, -1);
    element_type = new ElementTypeId[chunkSize+1];
    std::vector<int> elementIndices (chunkSize, -1);
    std::vector<int> faceElementIndices (chunkSize, -1);

    //chunkInfo stores the number of elements and number of face elements
    int chunkInfo[2];

#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkSize; i++) {
        id[i] = -1;
        tag[i] = -1;
        element_type[i] = NoRef;
    }

    /* Each worker receives messages */
    MPI_Status status;

    MPI_Recv(&errorFlag, 1, MPI_INT,0, 81719, mpi_info->comm, &status);
    if (errorFlag) {
        return errorFlag;
    }
    MPI_Recv(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, 0, 81720, mpi_info->comm, &status);
    MPI_Recv(id, chunkSize, MPI_INT, 0, 81721, mpi_info->comm, &status);
    MPI_Recv(tag, chunkSize, MPI_INT, 0, 81722, mpi_info->comm, &status);
    MPI_Recv(element_type, chunkSize, MPI_INT, 0, 81723, mpi_info->comm, &status);
    MPI_Recv(chunkInfo, 2, MPI_INT, 0, 81724, mpi_info->comm, &status);
    chunkElements = chunkInfo[0];
    chunkFaceElements = chunkInfo[1];
    MPI_Recv(&(elementIndices[0]), chunkElements, MPI_INT, 0, 81725, mpi_info->comm,&status);
    MPI_Recv(&(faceElementIndices[0]), chunkFaceElements, MPI_INT, 0, 81726, mpi_info->comm,&status);


    MPI_Bcast(&errorFlag, 1, MPI_INT,  0, mpi_info->comm);
    if (errorFlag) {
        return errorFlag;
    }

    // all elements have been read and shared, now we have to identify the
    // elements for finley
    int numNodes[3] = {0,0,0};
    MPI_Bcast(numNodes, 3, MPI_INT,  0, mpi_info->comm);
    final_element_type = static_cast<ElementTypeId>(numNodes[0]);
    final_face_element_type = static_cast<ElementTypeId>(numNodes[1]);
    contact_element_type = static_cast<ElementTypeId>(numNodes[2]);

    refElements.reset(new ReferenceElementSet(final_element_type, order, reduced_order));
    refFaceElements.reset(new ReferenceElementSet(final_face_element_type, order, reduced_order));
    refContactElements.reset(new ReferenceElementSet(contact_element_type, order, reduced_order));
    refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));
    mesh_p->Elements=new ElementFile(refElements, mpi_info);
    mesh_p->FaceElements=new ElementFile(refFaceElements, mpi_info);
    mesh_p->ContactElements=new ElementFile(refContactElements, mpi_info);
    mesh_p->Points=new ElementFile(refPoints, mpi_info);
    mesh_p->Elements->allocTable(chunkElements);
    mesh_p->FaceElements->allocTable(chunkFaceElements);
    mesh_p->ContactElements->allocTable(0);
    mesh_p->Points->allocTable(0);
    mesh_p->Elements->minColor=0;
    mesh_p->Elements->maxColor=chunkElements-1;
    mesh_p->FaceElements->minColor=0;
    mesh_p->FaceElements->maxColor=chunkFaceElements-1;
    mesh_p->ContactElements->minColor=0;
    mesh_p->ContactElements->maxColor=0;
    mesh_p->Points->minColor=0;
    mesh_p->Points->maxColor=0;

#pragma omp parallel for schedule(static)
    for(int e = 0; e < chunkElements; e++) {
        mesh_p->Elements->Id[e]=id[elementIndices[e]];
        mesh_p->Elements->Tag[e]=tag[elementIndices[e]];
        mesh_p->Elements->Color[e]=elementIndices[e];
        mesh_p->Elements->Owner[e]=mpi_info->rank;
        for (int j = 0; j < mesh_p->Elements->numNodes; ++j)  {
            int vertex = vertices[INDEX2(j,elementIndices[e],MAX_numNodes_gmsh)];
            mesh_p->Elements->Nodes[INDEX2(j, e, mesh_p->Elements->numNodes)]=vertex;
        }
    }

#pragma omp parallel for schedule(static)
    for (int e = 0; e < chunkFaceElements; e++) {
        mesh_p->FaceElements->Id[e]=id[faceElementIndices[e]];
        mesh_p->FaceElements->Tag[e]=tag[faceElementIndices[e]];
        mesh_p->FaceElements->Color[e]=e;
        mesh_p->FaceElements->Owner[e]=mpi_info->rank;
        for (int j=0; j<mesh_p->FaceElements->numNodes; ++j) {
            int faceVertex = vertices[INDEX2(j,faceElementIndices[e],MAX_numNodes_gmsh)];
            mesh_p->FaceElements->Nodes[INDEX2(j, e, mesh_p->FaceElements->numNodes)]= faceVertex;
        }
    }

    /* and clean up */
    delete[] id;
    delete[] tag;
    delete[] element_type;
    return errorFlag;
#endif // ESYS_MPI
}

int getElements(escript::JMPI& mpi_info, Mesh * mesh_p, FILE * fileHandle_p,
        std::string& errorMsg, bool useMacroElements, const std::string fname,
        int numDim, double version, int order, int reduced_order) {
    if (mpi_info->rank == 0) {
        return getElementsMaster(mpi_info, mesh_p, fileHandle_p,
                errorMsg, useMacroElements, fname,
                numDim, version, order, reduced_order);
    }
    return getElementsSlave(mpi_info, mesh_p, fileHandle_p,
                errorMsg, useMacroElements, fname,
                numDim, version, order, reduced_order);
}

int gather_nodes(FILE *f, std::map<int,int>& tags, std::string& errorMsg,
                 int dim, double version, const char *fname)
{
    int numNodes=0;
    std::vector<char> line;
    if (!get_line(line, f))
        return EARLY_EOF;
    int scan_ret = sscanf(&line[0], "%d", &numNodes);
    if (scan_ret == EOF)
        return EARLY_EOF;
    for (int node = 0; node < numNodes; node++) {
        int tmp = 0;
        if (!get_line(line, f))
            return EARLY_EOF;
        int scan_ret = sscanf(&line[0], "%d", &tmp);
        if (scan_ret == EOF) {
            return EARLY_EOF;
        } else if (scan_ret != 1) {
            errorMsg = "malformed meshfile";
            return THROW_ERROR;
        }
        tags[tmp] = -1;
    }
    if (!get_line(line, f))
        return EARLY_EOF;
    if (!is_endnode_string(&line[0])) {
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
    if (!get_line(line, f))
        return EARLY_EOF;
    scan_ret = sscanf(&line[0], "%d\n", &numElements);
    if (scan_ret == EOF) {
        return EARLY_EOF;
    } else if (scan_ret != 1) {
        errorMsg = "readGmsh: malformed mesh file";
        return THROW_ERROR;
    }
    struct ElementInfo e;
    std::vector<int> v(MAX_numNodes_gmsh, -1);
    e.vertex = &v[0];

    for (int element = 0; element < numElements; element++) {
        getSingleElement(f, dim, version, e, errorMsg, fname, false);
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
    return 0;
}

int getNodesMaster(escript::JMPI& mpi_info, Mesh *mesh_p, FILE *fileHandle_p,
        int numDim, std::string& errorMsg, std::map< int, int>& tags, int errorFlag)
{
    int numNodes=0;
    std::vector<char> line;
    if (!get_line(line, fileHandle_p))
        errorFlag = EARLY_EOF;
    int scan_ret = sscanf(&line[0], "%d", &numNodes);
    SSCANF_CHECK(scan_ret);
#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpi_info->size > 1) {
        MPI_Bcast(&numNodes, 1, MPI_INT,  0, mpi_info->comm);
    }
#endif
    int chunkSize;
    if(mpi_info->size > 1) {
        chunkSize = (numNodes / mpi_info->size) + 1;
    } else {
        chunkSize = (numNodes / mpi_info->size);
    }
    int totalNodes=0, chunkNodes = 0;
    int *tempInts = new int[chunkSize+1];        /* Stores the integer message data */
    double *tempCoords = new double[chunkSize*numDim]; /* Stores the double message data */

    for (int nextCPU = mpi_info->size-1; nextCPU >= 0; nextCPU--) {
//intialise arrays
#pragma omp parallel for schedule(static)
        for (int i=0; i<chunkSize+1; i++)
            tempInts[i] = -1;
#pragma omp parallel for schedule(static)
        for (int i=0; i<chunkSize*numDim; i++)
            tempCoords[i] = -1.0;
        if (!errorFlag) {
            if (nextCPU ==0) {
                chunkSize = numNodes-totalNodes;
            }
            //read in chunksize nodes
            for (chunkNodes=0; chunkNodes<chunkSize; chunkNodes++) {
                if(totalNodes > numNodes) {
                    std::stringstream ss;
                    ss << "readGmsh: too many nodes (" << totalNodes << " < "
                        << numNodes << ")";
                    errorMsg = ss.str();
                    errorFlag = THROW_ERROR;
                    break;
                }
                std::vector<char> line;
                if (!get_line(line, fileHandle_p))
                    errorFlag = EARLY_EOF;
                
                if (is_endnode_string(&line[0])) {
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
                }
                totalNodes++;
            }
        }
#ifdef ESYS_MPI
        if (nextCPU != 0) {
            /* if there was an error, we have to stop them waiting for more */
            MPI_Send(&errorFlag, 1, MPI_INT, nextCPU, 81719, mpi_info->comm);
            /* send out this chunk of mesh to the next waiting node */
            if (!errorFlag) {
                tempInts[chunkSize] = chunkNodes;   /* The message has one more int to send chunkNodes */
                MPI_Send(tempInts, chunkSize+1, MPI_INT, nextCPU, 81720, mpi_info->comm);
                MPI_Send(tempCoords, chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpi_info->comm);
            }
        }
#endif
    }


#ifdef ESYS_MPI
    if (mpi_info->size > 1) {
        MPI_Bcast(&errorFlag, 1, MPI_INT, 0, mpi_info->comm);
    }
#endif
    if (errorFlag) {
        return errorFlag;
    }

    mesh_p->Nodes->allocTable(chunkNodes);

#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkNodes; i++) {
        mesh_p->Nodes->Id[i] = tempInts[i];
        mesh_p->Nodes->globalDegreesOfFreedom[i] = tempInts[i];
        int tag = tags[tempInts[i]];
        if (tag == -1) {
            mesh_p->Nodes->Tag[i] = tempInts[i]; //set tag to node label
        } else {
            mesh_p->Nodes->Tag[i] = tag; //set tag of element
        }
        for (int j=0; j<numDim; j++) {
            mesh_p->Nodes->Coordinates[INDEX2(j,i,numDim)] = tempCoords[i*numDim+j];
        }

    }

    delete[] tempInts;
    delete[] tempCoords;
    return errorFlag;
}

int getNodesSlave(escript::JMPI& mpi_info, Mesh *mesh_p, FILE *fileHandle_p,
        int numDim, std::string& errorMsg, std::map< int, int>& tags, int errorFlag)
{
#ifndef ESYS_MPI
    throw FinleyException("slave function called in non-MPI build");
#else

    if (mpi_info->size == 1)
        throw FinleyException("slave function called without master");

    int numNodes=0;

    // Broadcast numNodes if there are multiple mpi procs
    MPI_Bcast(&numNodes, 1, MPI_INT,  0, mpi_info->comm);
    int chunkSize = (numNodes / mpi_info->size) + 1, chunkNodes=0;
    int *tempInts = new int[chunkSize+1];        /* Stores the integer message data */
    double *tempCoords = new double[chunkSize*numDim]; /* Stores the double message data */
    /* Each worker receives two messages */
    MPI_Status status;
    MPI_Recv(&errorFlag, 1, MPI_INT,0, 81719, mpi_info->comm, &status);
    if(!errorFlag){
        MPI_Recv(tempInts, chunkSize+1, MPI_INT, 0, 81720, mpi_info->comm, &status);
        MPI_Recv(tempCoords, chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpi_info->comm, &status);
        chunkNodes = tempInts[chunkSize];   /* How many nodes are in this workers chunk? */
    }

    MPI_Bcast(&errorFlag, 1, MPI_INT, 0, mpi_info->comm);
    if (errorFlag) {
        return errorFlag;
    }

    mesh_p->Nodes->allocTable(chunkNodes);

#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkNodes; i++) {
        mesh_p->Nodes->Id[i] = tempInts[i];
        mesh_p->Nodes->globalDegreesOfFreedom[i] = tempInts[i];
        int tag = tags[tempInts[i]];
        if (tag == -1) {
            mesh_p->Nodes->Tag[i] = tempInts[i]; //set tag to node label
        } else {
            mesh_p->Nodes->Tag[i] = tag; //set tag of element
        }
        for (int j=0; j<numDim; j++) {
            mesh_p->Nodes->Coordinates[INDEX2(j,i,numDim)] = tempCoords[i*numDim+j];
        }

    }

    delete[] tempInts;
    delete[] tempCoords;
    return errorFlag;
#endif //#ifndef ESYS_MPI -> #else
}

int getNodes(escript::JMPI& mpi_info, Mesh *mesh_p, FILE *fileHandle_p,
        int numDim, std::string& errorMsg, std::map< int, int>& tags, int errorFlag)
{
    if (mpi_info->rank == 0)
        return getNodesMaster(mpi_info, mesh_p, fileHandle_p, numDim, errorMsg,
                tags, errorFlag);

    return getNodesSlave(mpi_info, mesh_p, fileHandle_p, numDim, errorMsg,
                tags, errorFlag);
}

int get_next_state(FILE *f, bool nodesRead, bool elementsRead, int *logicFlag) {
    std::vector<char> line;
    do {
        if (!get_line(line, f)) { //got no line
            //check to see we atleast have some nodes and elements
            if (!nodesRead) {
                //EOF before nodes section found
                return MISSING_NODES;
            } else if(!elementsRead){
                //EOF before elements section found
                return MISSING_ELEMENTS;
            }
            return SUCCESS; //EOF as expected
        }
//        if (line[0] != '$')
//            fprintf(stderr, "consuming line: %s", &line[0]);
    } while(line[0] != '$');

    if (!strncmp(&line[1], "MeshFormat", 10)) {
        *logicFlag = 1;
    } else if (is_node_string(&line[0])) {
        *logicFlag = 2;
    } else if (!strncmp(&line[1], "ELM", 3) || !strncmp(&line[1], "Elements", 8)) {
        *logicFlag = 3;
    } else if (!strncmp(&line[1], "PhysicalNames", 13)) {
        *logicFlag = 4;
    }
    return 0;
}

void recv_state(escript::JMPI& mpi_info, int *error, int *logic) {
#ifdef ESYS_MPI
    int flags[2] = {0};
    // Broadcast line
    MPI_Bcast(&flags, 2, MPI_INT, 0, mpi_info->comm);
    *error = flags[0];
    if (logic)
        *logic = flags[1];
#endif
}

void send_state(escript::JMPI& mpi_info, int error, int logic) {
#ifdef ESYS_MPI
    int flags[2] = {error, logic};
    // Broadcast line
    if (mpi_info->size > 1) {
        MPI_Bcast(&flags, 2, MPI_INT,  0, mpi_info->comm);
    }
#endif
}

int check_error(int error, FILE *f, const std::string& errorMsg)
{
    //handle errors
    switch(error) {
        case 0:
            break;
        case ERROR:
            throw FinleyException("ERROR set for unknown reason");
        case EARLY_EOF: //early eof while scanning
            throw escript::IOError("early eof while scanning");
        case MISSING_NODES:  //EOF before nodes section found
            throw escript::IOError("EOF before nodes section found");
        case MISSING_ELEMENTS:
            throw escript::IOError("EOF before elements section found");
        case THROW_ERROR: // throw errorMsg
            throw escript::IOError(errorMsg);
        case SUCCESS: // eof at apropriate time.
            if (f)
                fclose(f);
            break;
        default:
            throw FinleyException("an unknown error has occured in readGmsh");

    }
    return error;
}

Mesh* Mesh::readGmshMaster(escript::JMPI& mpi_info, const std::string fname, int numDim, int order,
                     int reduced_order, bool optimize, bool useMacroElements)
{
    double version = 1.0;
    bool nodesRead=false, elementsRead=false;
    int format = 0, size = sizeof(double), scan_ret,  errorFlag=0, logicFlag=0;
    std::vector<char> line;
    std::map<int,int> nodeTags;
    FILE* fileHandle_p = NULL;
    std::string errorMsg;

    size_t found = fname.find("\n");
    if (found != std::string::npos) {
        errorFlag=THROW_ERROR;
        send_state(mpi_info, errorFlag, logicFlag);
        throw escript::ValueError("readGmsh: filename contains newline characters!");
    }

    // allocate mesh
    Mesh* mesh_p = new Mesh(fname, numDim, mpi_info);

    // get file handle
    fileHandle_p = fopen(fname.c_str(), "r");
    if (fileHandle_p==NULL) {
        std::stringstream ss;
        ss << "readGmsh: opening file " << fname << " for reading failed.";
        errorMsg = ss.str();
        errorFlag=THROW_ERROR;
        send_state(mpi_info, errorFlag, logicFlag);
        throw escript::IOError(errorMsg);
    }
    // start reading
    while (errorFlag==0) {
        // find line starting with $
        logicFlag=0;
        errorFlag = get_next_state(fileHandle_p, nodesRead, elementsRead, &logicFlag);
        send_state(mpi_info, errorFlag, logicFlag);
        //pre-logic error check
        if (check_error(errorFlag, fileHandle_p, errorMsg) == SUCCESS)
            break;
        // format
        if (logicFlag == 1 && errorFlag == 0) {
            std::vector<char> fmt;
            if (!get_line(fmt, fileHandle_p))
                errorFlag = EARLY_EOF;
            scan_ret = sscanf(&fmt[0], "%lf %d %d\n", &version, &format, &size);
            SSCANF_CHECK(scan_ret);
        }
        // nodes are read
        else if (logicFlag == 2 && errorFlag ==0) {
            nodesRead=true;
            std::vector<int> sendable_map;
            long current = ftell(fileHandle_p);
            errorFlag = gather_nodes(fileHandle_p, nodeTags, errorMsg,
                    numDim, version, fname.c_str());
            if (!errorFlag && fseek(fileHandle_p, current, SEEK_SET) < 0) {
                errorMsg = "Error in file operation";
                errorFlag = THROW_ERROR;
            }
            send_state(mpi_info, errorFlag, logicFlag);
            check_error(errorFlag, fileHandle_p, errorMsg);
#ifdef ESYS_MPI
            int mapsize = 2*nodeTags.size();
            sendable_map.resize(mapsize);
            std::map<int,int>::iterator i = nodeTags.begin();
            for (int j = 0; i != nodeTags.end(); i++, j += 2) {
                sendable_map[j] = i->first;
                sendable_map[j + 1] = i->second;
            }
            if (mpi_info->size > 1) {
                MPI_Bcast(&mapsize, 1, MPI_INT, 0, mpi_info->comm);
                sendable_map.resize(mapsize);
                MPI_Bcast(&sendable_map[0], mapsize, MPI_INT, 0, mpi_info->comm);
                for (int j = 0; j < mapsize; j += 2)
                    nodeTags[sendable_map[j]] = sendable_map[j + 1];
            }
#endif
            errorFlag = getNodes(mpi_info, mesh_p, fileHandle_p, numDim,
                    errorMsg, nodeTags, errorFlag);
        }
        // elements
        else if(logicFlag==3 && errorFlag ==0) {
            elementsRead=true;
            errorFlag=getElements(mpi_info, mesh_p, fileHandle_p, errorMsg,
                    useMacroElements, fname, numDim, version, order,
                    reduced_order);
        }
        // name tags (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca )
        else if (logicFlag==4 && errorFlag == 0) {
            std::vector<char> names;
            if (!get_line(names, fileHandle_p))
                errorFlag = EARLY_EOF;
            int numNames=0;
            scan_ret = sscanf(&names[0], "%d", &numNames);
            SSCANF_CHECK(scan_ret);
#ifdef ESYS_MPI
            // Broadcast numNames if there are multiple mpi procs
            if (mpi_info->size > 1) {
                MPI_Bcast(&numNames, 1, MPI_INT,  0, mpi_info->comm);
            }
#endif
            for (int i = 0; i < numNames; i++) {
                std::vector<char> line;
                char name[1024] = {0};
                if (!get_line(line, fileHandle_p))
                    errorFlag = EARLY_EOF;
                int tag_info[2] = {0};
                char *position = &line[0];
                //skip the first int, it's the physical dimension
                if (next_space(&position, 1) == NULL 
                        || sscanf(position, "%d", tag_info) != 1 
                        || next_space(&position, 1) == NULL
                        || sscanf(position, "%s", name) != 1) {
                    errorFlag = ERROR;
                }
                name[strlen(name)-1]='\0'; //strip trailing "
                //mpi broadcast the tag info

#ifdef ESYS_MPI
                if (mpi_info->size > 1) {
                    tag_info[1]=strlen(name) + 1; //include \0
                    MPI_Bcast(tag_info, 2, MPI_INT,  0, mpi_info->comm);
                    MPI_Bcast(&name, tag_info[1], MPI_CHAR,  0, mpi_info->comm);
                }
#endif
                mesh_p->addTagMap(name+1, tag_info[0]); //skip leading "

            }
        }

        if (!get_line(line, fileHandle_p)) {
            errorFlag = EARLY_EOF;
        }
        if (line[0] != '$') {
            errorFlag = THROW_ERROR;
            std::stringstream ss;
            ss << "readGmsh: expected closing tag, got '"
                << &line[0] << "'...";
            errorMsg = ss.str();
        }
        send_state(mpi_info, errorFlag, logicFlag);
        //post logic error check, throws if relevant
        check_error(errorFlag, fileHandle_p, errorMsg);
    }
    // resolve id's
    mesh_p->resolveNodeIds();
    // rearrange elements
    mesh_p->prepare(optimize);
    return mesh_p;
}

Mesh* Mesh::readGmshSlave(escript::JMPI& mpi_info, const std::string fname, int numDim, int order,
                     int reduced_order, bool optimize, bool useMacroElements)
{
#ifndef ESYS_MPI
    throw FinleyException("slave function called in non-MPI build");
#else
    if (mpi_info->size == 1)
        throw FinleyException("slave function called but only one process");

    double version = 1.0;
    int errorFlag=0, logicFlag=0;
    int numNames=0;
    int i, tag_info[2];
    char name[1024];
    std::string errorMsg;
    std::map<int,int> nodeTags;
    FILE * fileHandle_p = NULL;

    // allocate mesh
    Mesh* mesh_p = new Mesh(fname, numDim, mpi_info);

    // get file handle
    /* start reading */
    while (errorFlag != SUCCESS) {
        logicFlag = 0;
        //pre logic state fetch
        recv_state(mpi_info, &errorFlag, &logicFlag);
        if (check_error(errorFlag, NULL, errorMsg) == SUCCESS)
            break;
         
        /* format */
        /* nodes are read */
        if (logicFlag == 2) {
            int mapsize = 0;
            std::vector<int> sendable_map;
            recv_state(mpi_info, &errorFlag, &logicFlag);
            check_error(errorFlag, NULL, errorMsg);
            MPI_Bcast(&mapsize, 1, MPI_INT, 0, mpi_info->comm);
            sendable_map.resize(mapsize);
            MPI_Bcast(&sendable_map[0], mapsize, MPI_INT, 0, mpi_info->comm);
            for (int j = 0; j < mapsize; j += 2)
                nodeTags[sendable_map[j]] = sendable_map[j + 1];

            errorFlag = getNodes(mpi_info, mesh_p, fileHandle_p, numDim,
                    errorMsg, nodeTags, errorFlag);
        }

        /* elements */
        else if(logicFlag==3) {
            errorFlag=getElements(mpi_info, mesh_p, fileHandle_p, errorMsg, useMacroElements,
                    fname, numDim, version, order, reduced_order);
        }
         /* name tags (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca ) */
        else if (logicFlag==4) {
            // Broadcast numNames if there are multiple mpi procs
            MPI_Bcast(&numNames, 1, MPI_INT,  0, mpi_info->comm);
            for (i = 0; i < numNames; i++) {
                //mpi broadcast the tag info
                tag_info[0]=0;
                tag_info[1]=0;
                MPI_Bcast(tag_info, 2, MPI_INT,  0, mpi_info->comm);
                MPI_Bcast(&name, tag_info[1], MPI_CHAR,  0, mpi_info->comm); //strlen + 1 for null terminator
                mesh_p->addTagMap(&name[1], tag_info[0]);
            }
        }
        //post logic error check
        recv_state(mpi_info, &errorFlag, &logicFlag);
        if (check_error(errorFlag, NULL, errorMsg) == SUCCESS)
            break;
    //end while loop
    }

    // resolve id's
    mesh_p->resolveNodeIds();
    // rearrange elements
    mesh_p->prepare(optimize);
    return mesh_p;
#endif // ESYS_MPI
}


Mesh* Mesh::readGmsh(escript::JMPI& mpi_info, const std::string fname,
                     int numDim, int order, int reduced_order, bool optimize,
                     bool useMacroElements)
{
    if (mpi_info->rank == 0)
        return readGmshMaster(mpi_info, fname, numDim, order, reduced_order,
                optimize, useMacroElements);

    return readGmshSlave(mpi_info, fname, numDim, order, reduced_order,
            optimize, useMacroElements);
}

} // namespace finley

