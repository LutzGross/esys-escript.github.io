
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

  Finley: read mesh from gmsh file

*****************************************************************************/

#include "Mesh.h"
#include <cstdio>
#include "CPPAdapter/FinleyAdapterException.h"

//can't return because the flag need to be shared across all nodes
#define FSCANF_CHECK(scan_ret) { if (scan_ret == EOF) { errorFlag = 1;} }
#define MAX_numNodes_gmsh 20

#define EARLY_EOF 1
#define THROW_ERROR 4
#define ERROR 6
/*
    error flags include:

        0 - all ok
        1 - early eof
        2 - EOF before nodes section found
        3 - EOF before elements section found
        4 - throw error_msg
        5 - eof at apropriate time.
        6 - !noError
*/

struct ElementInfo {
    finley::ElementTypeId type;
    int id;
    int dim;
    int *vertex;
    int tag;
};

bool is_node_string(char *line) {
    if (line == NULL)
        return false;
    return !strncmp(line, "$NOD", 4) || !strncmp(line, "$NOE", 4) 
            || !strncmp(line, "$Nodes", 6);
}

bool is_endnode_string(char *line) {
    if (line == NULL)
        return false;
    return !strncmp(line, "$ENDNOD", 7) || !strncmp(line, "$ENDNOE", 7) 
            || !strncmp(line, "$EndNodes", 9);
}

namespace finley {

int getSingleElement(FILE *f, int dim, double version, struct ElementInfo& e,
        char *error_msg, const char *fname, bool useMacroElements)
{
    int gmsh_type = -1;
    int numTags=0;
    int scan_ret = fscanf(f, "%d %d", &e.id, &gmsh_type);
    if (scan_ret == EOF)
        return EARLY_EOF;
    else if (scan_ret != 2) {
        sprintf(error_msg, "malformed mesh file");
        return THROW_ERROR;
    }
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
            e.type=NoRef;
            e.dim=-1;
            sprintf(error_msg,"Unexpected gmsh element type %d in mesh file %s.", gmsh_type, fname);
            return THROW_ERROR;
    }
    if (version <= 1.0){
        int tmp = 0;
        scan_ret = fscanf(f, "%d %*d %d", &e.tag, &tmp);
        if (scan_ret == EOF)
            return EARLY_EOF;
        if (tmp != numNodesPerElement) {
            sprintf(error_msg,"Illegal number of nodes for element %d in mesh file %s.", e.id, fname);
            return THROW_ERROR;
        }
    } else {
        scan_ret = fscanf(f, "%d", &numTags);
        if (scan_ret == EOF)
            return EARLY_EOF;
        e.tag = 1;
        for (int j = 0; j < numTags; j++){
            int tmp = 0;
            scan_ret = fscanf(f, "%d", &tmp);
            if (scan_ret == EOF)
                return EARLY_EOF;
            if (j == 0) {
              e.tag = tmp;
            }
            /* ignore any other tags, second tag would be elementary id,
             third tag would be partition id */
        }
    }

    if (!noError()) {
        return ERROR;
    }
    for(int j = 0; j < numNodesPerElement; j++) {
        scan_ret = fscanf(f, "%d", e.vertex+j);
        if (scan_ret == EOF)
            return EARLY_EOF;
    }
    return 0;
}



int getElements(esysUtils::JMPI& mpi_info, Mesh * mesh_p, FILE * fileHandle_p, 
        char * error_msg, bool useMacroElements, const std::string fname,
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
    if (mpi_info->rank == 0) {
        int scan_ret = fscanf(fileHandle_p, "%d", &totalNumElements);
        FSCANF_CHECK(scan_ret);
    }
    
#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpi_info->size > 1) {
        int msg[1];
        if (mpi_info->rank == 0) {
            msg[0] = totalNumElements;
        } else {
            msg[0] = 0;
        }
        MPI_Bcast(msg, 1, MPI_INT, 0, mpi_info->comm);
        totalNumElements = msg[0];
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
    int chunkInfo[2];//chunkInfo stores the number of element and number of face elements
    int cpuId = 0;
#endif


#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkSize; i++) {
        id[i] = -1;
        tag[i] = -1;
        element_type[i] = NoRef;
    }
    
    if (mpi_info->rank == 0) {
      /* read all in */
        for(int e = 0, count = 0; e < totalNumElements; e++, count++) {
            struct ElementInfo element = {NoRef, 0, 0, &vertices[count*MAX_numNodes_gmsh], 0};
            getSingleElement(fileHandle_p, numDim, version,
                    element, error_msg, fname.c_str(), useMacroElements);
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
                   sprintf(error_msg,"Finley can handle a single type of internal elements only.");
                   errorFlag = THROW_ERROR;
               }
               elementIndices[chunkElements]=count;
               numElements++;
               chunkElements++;
            } else if (element.dim == numDim-1) {
               if (final_face_element_type == NoRef) {
                  final_face_element_type = element.type;
               } else if (final_face_element_type != element.type) {
                   sprintf(error_msg,"Finley can handle a single type of face elements only.");
                   errorFlag = THROW_ERROR;
               }
               faceElementIndices[chunkFaceElements]=count;
               numFaceElements++;
               chunkFaceElements++;
            } else{
                chunkOtherElements++;
            }
#ifdef ESYS_MPI
            if(count < chunkSize - 1)
                continue;
            cpuId++;
            chunkInfo[0]=chunkElements;
            chunkInfo[1]=chunkFaceElements;

            if(cpuId >= mpi_info->size) {
                continue;
            }
            if(!errorFlag){         
                MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
                MPI_Send(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, cpuId, 81720, mpi_info->comm);
                MPI_Send(id, chunkSize, MPI_INT, cpuId, 81721, mpi_info->comm);
                MPI_Send(tag, chunkSize, MPI_INT, cpuId, 81722, mpi_info->comm);
                MPI_Send(element_type, chunkSize, MPI_INT, cpuId, 81723, mpi_info->comm);
                MPI_Send(chunkInfo, 2, MPI_INT, cpuId, 81724, mpi_info->comm);
                MPI_Send(&(elementIndices[0]), chunkElements, MPI_INT, cpuId, 81725, mpi_info->comm);
                MPI_Send(&(faceElementIndices[0]), chunkFaceElements, MPI_INT, cpuId, 81726, mpi_info->comm);
            } else{
                //let the remaining processes know something has gone wrong
                for(; cpuId<mpi_info->size; cpuId++) {
                    MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
                }                    
                break;
            }
            
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

            //end elment reading for loop
            if (errorFlag)
            {
                break;
            }
        }
    } else {
#ifdef ESYS_MPI
        /* Each worker receives messages */
        if(mpi_info->size>1){
            MPI_Status status;

            MPI_Recv(&errorFlag, 1, MPI_INT,0, 81719, mpi_info->comm, &status);
            if(!errorFlag){
                MPI_Recv(&vertices[0], chunkSize*MAX_numNodes_gmsh, MPI_INT, 0, 81720, mpi_info->comm, &status);
                MPI_Recv(id, chunkSize, MPI_INT, 0, 81721, mpi_info->comm, &status);
                MPI_Recv(tag, chunkSize, MPI_INT, 0, 81722, mpi_info->comm, &status);
                MPI_Recv(element_type, chunkSize, MPI_INT, 0, 81723, mpi_info->comm, &status);
                MPI_Recv(chunkInfo, 2, MPI_INT, 0, 81724, mpi_info->comm, &status);
                chunkElements = chunkInfo[0];
                chunkFaceElements = chunkInfo[1];
                MPI_Recv(&(elementIndices[0]), chunkElements, MPI_INT, 0, 81725, mpi_info->comm,&status);
                MPI_Recv(&(faceElementIndices[0]), chunkFaceElements, MPI_INT, 0, 81726, mpi_info->comm,&status);
            } else {
                return errorFlag;
            }
        }
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
    if (!noError()) {
        return ERROR;
    }
    if (mpi_info->rank == 0) {
        /* first we have to identify the elements to define Elements and FaceElements */
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
    } 
    
#ifdef ESYS_MPI
        // Broadcast numNodes if there are multiple mpi procs
        if (mpi_info->size > 1) {
            int msg[3];
            if (mpi_info->rank == 0) {
                msg[0] = final_element_type;
                msg[1] = final_face_element_type;
                msg[2] = contact_element_type;
            } else {
                msg[0] = 0;
                msg[1] = 0;
                msg[2] = 0;
            }
            MPI_Bcast(msg, 3, MPI_INT,  0, mpi_info->comm);
            final_element_type = static_cast<ElementTypeId>(msg[0]);
            final_face_element_type = static_cast<ElementTypeId>(msg[1]);
            contact_element_type = static_cast<ElementTypeId>(msg[2]);
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
    

    if (!noError())
        return ERROR;
    
    mesh_p->Elements->allocTable(chunkElements);
    mesh_p->FaceElements->allocTable(chunkFaceElements);
    mesh_p->ContactElements->allocTable(0);
    mesh_p->Points->allocTable(0);

    if (!noError())
        return ERROR;

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

int gather_nodes(FILE *f, std::map<int,int>& tags, char *error_msg,
        int dim, double version, const char *fname)
{
    int numNodes=0;

    if (fscanf(f, "%d", &numNodes) == EOF)
        return EARLY_EOF;
    for (int node = 0; node < numNodes; node++) {
        int tmp = 0;
        int scan_ret = fscanf(f, "%d %*e %*e %*e\n", &tmp);
        if (scan_ret == EOF) {
            return EARLY_EOF;
        } else if (scan_ret != 1) {
            sprintf(error_msg, "malformed meshfilea");
            return THROW_ERROR;
        }
        tags[tmp] = -1;
    }
    char line[LenString_MAX] = {0};
    if(fgets(line,LenString_MAX,f) == NULL) {
        return 1;
    }
    if (!is_endnode_string(line)) {
        sprintf(error_msg, "malformed meshfile, expected '$EndNodes', got '%s'", line);
        return THROW_ERROR;
    }
    if(fgets(line,LenString_MAX,f) == NULL) {
        return 1;
    }
    if (strncmp(line, "$ELM", 4) && strncmp(line, "$Elements", 9)) {
        sprintf(error_msg, "malformed meshfile, expected '$Elements', got '%s'", line);
        return THROW_ERROR;
    }
    int numElements = -1;
    int scan_ret = fscanf(f, "%d\n", &numElements);
    if (scan_ret == EOF) {
        return EARLY_EOF;
    } else if (scan_ret != 1) {
        sprintf(error_msg, "malformed meshfile");
        return THROW_ERROR;
    }
    struct ElementInfo e;
    std::vector<int> v(MAX_numNodes_gmsh, -1);
    e.vertex = &v[0];
    
    for (int element = 0; element < numElements; element++) {
        getSingleElement(f, dim, version, e, error_msg, fname, false);
        for (int i = 0; i < MAX_numNodes_gmsh && v[i] >= 0; i++) {
            std::map<int,int>::iterator it = tags.find(v[i]);
            if (it == tags.end()) {
                sprintf(error_msg, "element contains unknown node");
                return THROW_ERROR;
            }
            // the last element containing a node is the tag that stays
            tags[v[i]] = e.tag;
        }
    }
    return 0;
}

int getNodes(esysUtils::JMPI& mpi_info, Mesh *mesh_p, FILE *fileHandle_p,
        int numDim, char *error_msg, std::map< int, int>& tags, int errorFlag)
{
    int scan_ret, numNodes=0;
    double rtmp0, rtmp1;

    if (mpi_info->rank == 0) {  /* Master */
        scan_ret = fscanf(fileHandle_p, "%d", &numNodes);
        FSCANF_CHECK(scan_ret);
    }
#ifdef ESYS_MPI
    // Broadcast numNodes if there are multiple mpi procs
    if (mpi_info->size > 1) {
        int msg[1];
        if (mpi_info->rank == 0) {
            msg[0] = numNodes;
        } else {
            msg[0] = 0;
        }
        MPI_Bcast(msg, 1, MPI_INT,  0, mpi_info->comm);
        numNodes = msg[0];
    }
#endif
    int chunkSize = (numNodes / mpi_info->size) + 1, totalNodes=0, chunkNodes=0,  nextCPU=1;
    int *tempInts = new int[chunkSize+1];        /* Stores the integer message data */
    double *tempCoords = new double[chunkSize*numDim]; /* Stores the double message data */
    if (mpi_info->rank == 0) {  /* Master */
        while(1){
#pragma omp parallel for schedule(static)
            for (int i=0; i<chunkSize+1; i++)
                tempInts[i] = -1;
#pragma omp parallel for schedule(static)
            for (int i=0; i<chunkSize*numDim; i++)
                tempCoords[i] = -1.0;

            if(!errorFlag){
                chunkNodes=0;
                for(int i=0;i<chunkSize;i++){
                    if(totalNodes>=numNodes) break;//sprintf(error_msg, "more "); errorFlag=4;
                    if (1 == numDim) {
                        scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &tempInts[i], &tempCoords[0+i*numDim], &rtmp0, &rtmp1);
                        FSCANF_CHECK(scan_ret);
                    } else if (2 == numDim) {
                        scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &tempInts[i], &tempCoords[0+i*numDim], &tempCoords[1+i*numDim], &rtmp1);
                        FSCANF_CHECK(scan_ret);
                    } else if (3 == numDim) {
                        scan_ret = fscanf(fileHandle_p, "%d %le %le %le\n", &tempInts[i], &tempCoords[0+i*numDim], &tempCoords[1+i*numDim], &tempCoords[2+i*numDim]);
                        FSCANF_CHECK(scan_ret);
                    }
                    totalNodes++; /* When do we quit the infinite loop? */
                    chunkNodes++; /* How many nodes do we actually have in this chunk? It may be smaller than chunkSize. */
                }
            }
#ifdef ESYS_MPI
            if(mpi_info->size>1 && (nextCPU < mpi_info->size)) { 
                /* Eventually we'll send chunkSize nodes to each CPU numbered 1 ... mpi_info->size-1, here goes one of them */
                MPI_Send(&errorFlag, 1, MPI_INT, nextCPU, 81719, mpi_info->comm);
                if(!errorFlag){
                        tempInts[chunkSize] = chunkNodes;   /* The message has one more int to send chunkNodes */
                        MPI_Send(tempInts, chunkSize+1, MPI_INT, nextCPU, 81720, mpi_info->comm);
                        MPI_Send(tempCoords, chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpi_info->comm);
                }
            }
#endif
            nextCPU++;
            /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
            if (nextCPU > mpi_info->size) break; /* End infinite loop */
        }
    } else {
#ifdef ESYS_MPI
        /* Each worker receives two messages */
        MPI_Status status;
        MPI_Recv(&errorFlag, 1, MPI_INT,0, 81719, mpi_info->comm, &status);
        if(!errorFlag){
            MPI_Recv(tempInts, chunkSize+1, MPI_INT, 0, 81720, mpi_info->comm, &status);
            MPI_Recv(tempCoords, chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpi_info->comm, &status);
            chunkNodes = tempInts[chunkSize];   /* How many nodes are in this workers chunk? */
        }
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

    if (!noError()) return ERROR;
    mesh_p->Nodes->allocTable(chunkNodes);
    if (!noError()) return ERROR;

#pragma omp parallel for schedule(static)
    for (int i=0; i<chunkNodes; i++) {
        mesh_p->Nodes->Id[i] = tempInts[i];
        mesh_p->Nodes->globalDegreesOfFreedom[i] = tempInts[i];
        mesh_p->Nodes->Tag[i] = tags[tempInts[i]]; //set tag of element
        for (int j=0; j<numDim; j++) {
            mesh_p->Nodes->Coordinates[INDEX2(j,i,numDim)] = tempCoords[i*numDim+j];
        }

    }

    delete[] tempInts;
    delete[] tempCoords; 
    return errorFlag;
}




Mesh* Mesh::readGmsh(esysUtils::JMPI& mpi_info, const std::string fname, int numDim, int order,
                     int reduced_order, bool optimize, bool useMacroElements)
{
    double version = 1.0;
    bool nodesRead=false, elementsRead=false;
    int format = 0, size = sizeof(double), scan_ret,  errorFlag=0, logicFlag=0;
    int numNames=0;
    int i, tag_info[2], itmp;
    char line[LenString_MAX+1], name[LenString_MAX+1];
    char error_msg[LenErrorMsg_MAX];
    std::map<int,int> nodeTags;
#ifdef Finley_TRACE
    double time0=timer();
#endif
    FILE * fileHandle_p = NULL;

    resetError();
    // allocate mesh
    Mesh* mesh_p = new Mesh(fname, numDim, mpi_info);

    // get file handle
    if (mpi_info->rank == 0) {
        fileHandle_p = fopen(fname.c_str(), "r");
        if (fileHandle_p==NULL) {
            sprintf(error_msg, "Opening Gmsh file %s for reading failed.", fname.c_str());
            throw FinleyAdapterException(error_msg);
        }
    }
    /* start reading */
    while(noError() && errorFlag==0) {
        /* find line starting with $ */
        if (mpi_info->rank == 0) {
            do {
                if(fgets(line,LenString_MAX,fileHandle_p) == NULL) {
                    //check to see we atleast have some nodes and elements, if we do end the outer loop otherwise throw early eof.
                    if (!nodesRead){
                        //EOF before nodes section found
                        errorFlag = 2; 
                    } else if(!elementsRead){
                        //EOF before elements section found
                        errorFlag = 3;
                    } else{
                        errorFlag=5;
                    }
                    break;
                }
            } while(line[0] != '$');

            if (!strncmp(&line[1], "MeshFormat", 10)) {
                logicFlag = 1;
            }

            else if (is_node_string(line)) {
                logicFlag = 2;
            }
            else if(!strncmp(&line[1], "ELM", 3) ||
               !strncmp(&line[1], "Elements", 8)  ) {
                logicFlag = 3;
            }
            else if (!strncmp(&line[1], "PhysicalNames", 13)) {
                logicFlag = 4;
            }
        }
#ifdef ESYS_MPI
        int flags[2];
        // Broadcast line
        if (mpi_info->size > 1) {
            if (mpi_info -> rank==0) {
                flags[0] = errorFlag;
                flags[1] = logicFlag;

            } else {
                flags[0] = 0;
                flags[1] = 0;
            }
            MPI_Bcast(&flags, 2, MPI_INT,  0, mpi_info->comm);
            errorFlag = flags[0];
            logicFlag = flags[1];
        }
#endif
        // MPI_Barrier(mpi_info->comm);
        /* format */
        if (mpi_info->rank==0 && logicFlag ==1) {
                scan_ret = fscanf(fileHandle_p, "%lf %d %d\n", &version, &format, &size);
                FSCANF_CHECK(scan_ret);
        }
        /* nodes are read */
        else if (logicFlag == 2) {

            nodesRead=true;
            int mapsize = 0;
            std::vector<int> sendable_map;
            if (mpi_info->rank == 0) {
                long current = ftell(fileHandle_p);
                errorFlag = gather_nodes(fileHandle_p, nodeTags, error_msg,
                        numDim, version, fname.c_str());
                if (!errorFlag && fseek(fileHandle_p, current, SEEK_SET) < 0) {
                        sprintf(error_msg, "Error in file operation");
                        errorFlag = THROW_ERROR;
                }
#ifdef ESYS_MPI
                mapsize = 2*nodeTags.size();
                sendable_map.resize(mapsize);
                std::map<int,int>::iterator i = nodeTags.begin();
                for (int j = 0; i != nodeTags.end(); i++, j += 2) {
                    sendable_map[j] = i->first;
                    sendable_map[j + 1] = i->second;
                }
            }
            MPI_Bcast(&mapsize, 1, MPI_INT, 0, mpi_info->comm);
            if (mpi_info->rank != 0) {
                sendable_map.resize(mapsize);
            }
            MPI_Bcast(&sendable_map[0], mapsize, MPI_INT, 0, mpi_info->comm);
            if (mpi_info->rank != 0) {
                for (int j = 0; j < mapsize; j += 2)
                    nodeTags[sendable_map[j]] = sendable_map[j + 1];
            }
#else
            }
#endif
            errorFlag = getNodes(mpi_info, mesh_p, fileHandle_p, numDim,
                    error_msg, nodeTags, errorFlag); 
            logicFlag=0;
        }
       
        /* elements */
        else if(logicFlag==3) {
            elementsRead=true;
            errorFlag=getElements(mpi_info, mesh_p, fileHandle_p, error_msg, useMacroElements, 
                    fname, numDim, version, order, reduced_order);
            logicFlag=0;
        }
         /* name tags (thanks to Antoine Lefebvre, antoine.lefebvre2@mail.mcgill.ca ) */
        else if (logicFlag==4) {
            
            if (! noError()) errorFlag=6;
            if(mpi_info->rank == 0) {
                scan_ret = fscanf(fileHandle_p, "%d", &numNames);
                FSCANF_CHECK(scan_ret);
            }
#ifdef ESYS_MPI
            // Broadcast numNodes if there are multiple mpi procs
            if (mpi_info->size > 1) {
                MPI_Bcast(&numNames, 1, MPI_INT,  0, mpi_info->comm);
            }
#endif

            for (i = 0; i < numNames; i++) {
                if(mpi_info->rank == 0) {
                    scan_ret = fscanf(fileHandle_p, "%d %d %s\n", &itmp, &(tag_info[0]), name);
                    FSCANF_CHECK(scan_ret);
                    //if (! (itmp == 2)) setError(IO_ERROR,"Mesh_readGmsh: expecting two entries per physical name.");
                    if ( strlen(name) < 3 ) setError(IO_ERROR,"Mesh_readGmsh: illegal tagname (\" missing?)");
                    if (! noError()) errorFlag=6 ;
                    name[strlen(name)-1]='\0';
                } 
                //mpi broadcast the tag info

#ifdef ESYS_MPI
                if (mpi_info->size > 1) {
                    if(mpi_info->rank==0){
                        tag_info[1]=strlen(name)+1;
                    } else{
                        tag_info[0]=0;
                        tag_info[1]=0;
                    }

                    MPI_Bcast(tag_info, 2, MPI_INT,  0, mpi_info->comm);
                    MPI_Bcast(&name, tag_info[1], MPI_CHAR,  0, mpi_info->comm); //strlen + 1 for null terminator
                }
#endif         
                
                mesh_p->addTagMap(&name[1], tag_info[0]);

            }
            logicFlag=0;
        }

        //read closing tag $EndTHING
        if(mpi_info->rank==0 && errorFlag!=5) {
            do {
                if(fgets(line,LenString_MAX,fileHandle_p) == NULL) {
                    errorFlag = 1;
                    break;
                }
            } while( line[0]!='$' );
        }

        //handle errors
        if(errorFlag){
            switch(errorFlag) {
                case 1: //early eof while scanning
                    throw FinleyAdapterException("early eof while scanning");
                    break;
                case 2:  //EOF before nodes section found
                    throw FinleyAdapterException("EOF before nodes section found");
                    break;
                case 3:
                    throw FinleyAdapterException("EOF before elements section found");
                    break;
                case 4: // throw error_msg
                    throw FinleyAdapterException(error_msg);
                    break;
                case 5: // eof at apropriate time.
                    if(mpi_info->rank == 0) {
                        fclose(fileHandle_p);
                    }                    
                    break;
                case 6:
                    // throw FinleyAdapterException("not noError");
                    return NULL; 
                default:
                    throw FinleyAdapterException("an unknown error has occured in readGmsh");
            }
            if(errorFlag == 5) {
                break;
            }
        }


    //end while loop
    }

    // clean up
    if (!noError()) {
        delete mesh_p;
        return NULL;
    }
    // resolve id's
    if (noError())
        mesh_p->resolveNodeIds();
    // rearrange elements
    if (noError()) 
        mesh_p->prepare(optimize);

    if (!noError()) {
        delete mesh_p;
        return NULL;
    }
    return mesh_p;
 }

} // namespace finley

