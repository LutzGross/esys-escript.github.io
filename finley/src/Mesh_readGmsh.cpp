
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

namespace finley {

int getElements(esysUtils::JMPI& mpi_info, Mesh * mesh_p, FILE * fileHandle_p, char * error_msg, 
            bool useMacroElements, const std::string fname, int numDim,double version, int order, int reduced_order) {
    /* This function should read in the elements and distribute them to the apropriate process.
     *
     */
    int errorFlag=0, scan_ret, i, j, e,gmsh_type, element_dim, partition_id, itmp, elementary_id, numTags=0;
    ElementTypeId final_element_type = NoRef;
    ElementTypeId final_face_element_type = NoRef;
    ElementTypeId contact_element_type = NoRef;
    int numElements=0, numFaceElements=0, totalNumElements=0, numNodesPerElement=0, numNodesPerElement2;
    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    int *id, *tag, * vertices;
    ElementTypeId * element_type;
    if (mpi_info->rank == 0) {
        scan_ret = fscanf(fileHandle_p, "%d", &totalNumElements);
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
        MPI_Bcast(msg, 1, MPI_INT,  0, mpi_info->comm);
        totalNumElements = msg[0];
    }
#endif
    
    int chunkSize = totalNumElements / mpi_info->size + 1, chunkElements=0, chunkFaceElements=0, count=0;
    id = new int[chunkSize+1];
    tag = new int[chunkSize+1];
    vertices = new int[chunkSize*MAX_numNodes_gmsh];
    element_type = new ElementTypeId[chunkSize+1];
    std::vector<int> elementIndecies (chunkSize);
    std::vector<int> faceElementIndecies (chunkSize);
    


#ifdef ESYS_MPI        
    int chunkInfo[2];//chunkInfo stores the number of element and number of face elements
#endif

#pragma omp parallel for private (i) schedule(static)
    for (i=0; i<chunkSize*MAX_numNodes_gmsh; i++) vertices[i] = -1;

#pragma omp parallel for private (i) schedule(static)
    for (i=0; i<chunkSize; i++) {
        id[i] = -1;
        tag[i] = -1;
        element_type[i] = NoRef;
        elementIndecies[i]=-1;
        faceElementIndecies[i]=-1;
    }
    
    if (mpi_info->rank == 0) {
      /* read all in */
#ifdef ESYS_MPI
    int cpuId = 0;
#endif
        for(e = 0; e < totalNumElements; e++) {
            count = (chunkElements + chunkFaceElements);
            scan_ret = fscanf(fileHandle_p, "%d %d", &id[count], &gmsh_type);
            FSCANF_CHECK(scan_ret);
            switch (gmsh_type) {
                case 1:  /* line order 1 */
                    element_type[count]=Line2;
                    element_dim=1;
                    numNodesPerElement=2;
                    break;
                case 2:  /* triangle order 1 */
                    element_type[count]=Tri3;
                    numNodesPerElement= 3;
                    element_dim=2;
                    break;
                case 3:  /* quadrilateral order 1 */
                    element_type[count]=Rec4;
                    numNodesPerElement= 4;
                    element_dim=2;
                    break;
                case 4:  /* tetrahedron order 1 */
                    element_type[count]=Tet4;
                    numNodesPerElement= 4;
                    element_dim=3;
                    break;
                case 5:  /* hexahedron order 1 */
                    element_type[count]=Hex8;
                    numNodesPerElement= 8;
                    element_dim=3;
                    break;
                case 8:  /* line order 2 */
                    if (useMacroElements) {
                        element_type[count]=Line3Macro;
                    } else {
                        element_type[count]=Line3;
                    }
                    numNodesPerElement= 3;
                    element_dim=1;
                    break;
                case 9:  /* triangle order 2 */
                    if (useMacroElements) {
                         element_type[count]=Tri6Macro;
                    } else {
                         element_type[count]=Tri6;
                    }
                    numNodesPerElement= 6;
                    element_dim=2;
                    break;
                case 10:  /* quadrilateral order 2 */
                    if (useMacroElements) {
                        element_type[count]=Rec9Macro;
                    } else {
                        element_type[count]=Rec9;
                    }
                    numNodesPerElement= 9;
                    element_dim=2;
                    break;
                case 11:  /* tetrahedron order 2 */
                    if (useMacroElements) {
                        element_type[count]=Tet10Macro;
                    } else {
                        element_type[count]=Tet10;
                    }
                    numNodesPerElement= 10;
                    element_dim=3;
                    break;
                case 16:  /* rectangular order 2 */
                    element_type[count]=Rec8;
                    numNodesPerElement= 8;
                    element_dim=2;
                    break;
                case 17:  /* hexahedron order 2 */
                    element_type[count]=Hex20;
                    numNodesPerElement= 20;
                    element_dim=3;
                    break;
                case 15 :  /* point */
                    element_type[count]=Point1;
                    numNodesPerElement= 1;
                    element_dim=0;
                    break;
                default:
                   element_type[count]=NoRef;
                   sprintf(error_msg,"Unexpected gmsh element type %d in mesh file %s.", gmsh_type, fname.c_str());
                   errorFlag = 4;
            }
            if (element_dim == numDim) {
               if (final_element_type == NoRef) {
                  final_element_type = element_type[count];
               } else if (final_element_type != element_type[count]) {
                   sprintf(error_msg,"Finley can handle a single type of internal elements only.");
                   errorFlag = 4;
                   
               }
               elementIndecies[chunkElements]=count;
               numElements++;
               chunkElements++;
            } else if (element_dim == numDim-1) {
               if (final_face_element_type == NoRef) {
                  final_face_element_type = element_type[count];
               } else if (final_face_element_type != element_type[count]) {
                   sprintf(error_msg,"Finley can handle a single type of face elements only.");
                   errorFlag = 4;
               }
               faceElementIndecies[chunkFaceElements]=count;
               numFaceElements++;
               chunkFaceElements++;
            }
            if(version <= 1.0){
              scan_ret = fscanf(fileHandle_p, "%d %d %d", &tag[count], &elementary_id, &numNodesPerElement2);
              FSCANF_CHECK(scan_ret);
              partition_id = 1;
              if (numNodesPerElement2 != numNodesPerElement) {
                   sprintf(error_msg,"Illegal number of nodes for element %d in mesh file %s.", id[count], fname.c_str());
                   errorFlag = 4;
              }
            } else {
              scan_ret = fscanf(fileHandle_p, "%d", &numTags);
              FSCANF_CHECK(scan_ret);
              elementary_id = tag[count] = partition_id = 1;
              numNodesPerElement2=-1;
              for(j = 0; j < numTags; j++){
                scan_ret = fscanf(fileHandle_p, "%d", &itmp);
                FSCANF_CHECK(scan_ret);
                if (j == 0) {
                  tag[count] = itmp;
                } else if (j == 1) {
                  elementary_id = itmp;
                } else if (j == 2) {
                  partition_id = itmp;
                }
                /* ignore any other tags */
              }
            }

            if (!noError()) {
                errorFlag = 6;
            }
            // fprintf(stderr,"num elements=%d, chunkElements = %d\n",(numElements+numFaceElements),(chunkFaceElements+chunkElements));
            for(j = 0; j < numNodesPerElement; j++) {
              scan_ret = fscanf(fileHandle_p, "%d", &vertices[INDEX2(j,count,MAX_numNodes_gmsh)]);
              FSCANF_CHECK(scan_ret);
            }
            /* for tet10 the last two nodes need to be swapped */
            if ((element_type[count]==Tet10) || (element_type[count]==Tet10Macro)) {
                 itmp=vertices[INDEX2(9,count,MAX_numNodes_gmsh)];
                 vertices[INDEX2(9,count,MAX_numNodes_gmsh)]=vertices[INDEX2(8,count,MAX_numNodes_gmsh)];
                 vertices[INDEX2(8,count,MAX_numNodes_gmsh)]=itmp;
            }
#ifdef ESYS_MPI
            if(chunkElements+chunkFaceElements==chunkSize) {
                cpuId++;
                chunkInfo[0]=chunkElements;
                chunkInfo[1]=chunkFaceElements;

                if(cpuId < mpi_info->size) {
                    if(!errorFlag){                
                        MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
                        MPI_Send(vertices, chunkSize*MAX_numNodes_gmsh, MPI_INT, cpuId, 81720, mpi_info->comm);
                        MPI_Send(id, chunkSize, MPI_INT, cpuId, 81721, mpi_info->comm);
                        MPI_Send(tag, chunkSize, MPI_INT, cpuId, 81722, mpi_info->comm);
                        MPI_Send(element_type, chunkSize, MPI_INT, cpuId, 81723, mpi_info->comm);
                        MPI_Send(chunkInfo, 2, MPI_INT, cpuId, 81724, mpi_info->comm);
                        MPI_Send(&(elementIndecies[0]), chunkElements, MPI_INT, cpuId, 81725, mpi_info->comm);
                        MPI_Send(&(faceElementIndecies[0]), chunkFaceElements, MPI_INT, cpuId, 81726, mpi_info->comm);
                    } else{
                        //let the remaining processes know something has gone wrong
                        for(; cpuId<mpi_info->size; cpuId++) {
                            MPI_Send(&errorFlag, 1, MPI_INT, cpuId, 81719, mpi_info->comm);
                        }                    
                        break;
                    }
                    
                    // reset arrays for next cpu
                    // for i in vertices vertices[i]=-1 etc also use openMp
#pragma omp parallel for private (i) schedule(static)
                    for (i=0; i<chunkSize*MAX_numNodes_gmsh; i++) vertices[i] = -1;
#pragma omp parallel for private (i) schedule(static)
                    for (i=0; i<chunkSize; i++) {
                        id[i] = -1;
                        tag[i] = -1;
                        element_type[i] = NoRef;
                    }
                    chunkElements=0;
                    chunkFaceElements=0;
                }
            }
#endif
            
            //end elment reading for loop
        }
    } else {
#ifdef ESYS_MPI
        /* Each worker receives messages */
        if(mpi_info->size>1){
            MPI_Status status;

            MPI_Recv(&errorFlag, 1, MPI_INT,0, 81719, mpi_info->comm, &status);
            if(!errorFlag){
                MPI_Recv(vertices, chunkSize*MAX_numNodes_gmsh, MPI_INT, 0, 81720, mpi_info->comm, &status);
                MPI_Recv(id, chunkSize, MPI_INT, 0, 81721, mpi_info->comm, &status);
                MPI_Recv(tag, chunkSize, MPI_INT, 0, 81722, mpi_info->comm, &status);
                MPI_Recv(element_type, chunkSize, MPI_INT, 0, 81723, mpi_info->comm, &status);
                MPI_Recv(chunkInfo, 2, MPI_INT, 0, 81724, mpi_info->comm, &status);
                chunkElements = chunkInfo[0];
                chunkFaceElements = chunkInfo[1];
                MPI_Recv(&(elementIndecies[0]), chunkElements, MPI_INT, 0, 81725, mpi_info->comm,&status);
                MPI_Recv(&(faceElementIndecies[0]), chunkFaceElements, MPI_INT, 0, 81726, mpi_info->comm,&status);
            } else {
                return errorFlag;
            }
        }
#endif

    }

    // fprintf(stderr,"in elements rank=,%d\n",mpi_info->rank);    

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
        return 6;
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
    

    if (noError()) {
        mesh_p->Elements->allocTable(chunkElements);
        mesh_p->FaceElements->allocTable(chunkFaceElements);
        mesh_p->ContactElements->allocTable(0);
        mesh_p->Points->allocTable(0);
        if (noError()) {
            mesh_p->Elements->minColor=0;
            mesh_p->Elements->maxColor=chunkElements-1;
            mesh_p->FaceElements->minColor=0;
            mesh_p->FaceElements->maxColor=chunkFaceElements-1;
            mesh_p->ContactElements->minColor=0;
            mesh_p->ContactElements->maxColor=0;
            mesh_p->Points->minColor=0;
            mesh_p->Points->maxColor=0;
            

            chunkElements=0;
            chunkFaceElements=0;
            for(e = 0; e < chunkSize; e++) {
               if (element_type[e] == final_element_type) {
                  mesh_p->Elements->Id[chunkElements]=id[e];
                  mesh_p->Elements->Tag[chunkElements]=tag[e];
                  mesh_p->Elements->Color[chunkElements]=e;
                  mesh_p->Elements->Owner[chunkElements]=mpi_info->rank;
                  // fprintf(stderr,"element id%d: ",mesh_p->Elements->Id[chunkElements]);
                  for (j = 0; j<  mesh_p->Elements->numNodes; ++j)  {
                        mesh_p->Elements->Nodes[INDEX2(j, chunkElements, mesh_p->Elements->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                        // fprintf(stderr,"nodes[%d]=%d ",INDEX2(j, chunkElements, mesh_p->Elements->numNodes),vertices[INDEX2(j,e,MAX_numNodes_gmsh)]);
                  }
                  // fprintf(stderr,"\n");
                  chunkElements++;
               } else if (element_type[e] == final_face_element_type) {
                  mesh_p->FaceElements->Id[chunkFaceElements]=id[e];
                  mesh_p->FaceElements->Tag[chunkFaceElements]=tag[e];
                  mesh_p->FaceElements->Color[chunkFaceElements]=chunkFaceElements;
                  mesh_p->FaceElements->Owner[chunkFaceElements]=mpi_info->rank;
                  for (j=0; j<mesh_p->FaceElements->numNodes; ++j) {
                           mesh_p->FaceElements->Nodes[INDEX2(j, chunkFaceElements, mesh_p->FaceElements->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                  }
                  chunkFaceElements++;
               }
            }
        } else {
            return 6;
        }

    } else {
        return 6;
    }
    // fprintf(stderr,"in elements rank=,%d\n",mpi_info->rank);    

    /* and clean up */
    delete[] id;
    delete[] tag;
    delete[] element_type;
    delete[] vertices;
    return errorFlag;
}


int getNodes(esysUtils::JMPI& mpi_info, Mesh * mesh_p, FILE * fileHandle_p, int numDim, char * error_msg) {
    int i, j, scan_ret, numNodes=0, errorFlag=0;
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
#pragma omp parallel for private (i) schedule(static)
            for (i=0; i<chunkSize+1; i++) tempInts[i] = -1;
#pragma omp parallel for private (i) schedule(static)
            for (i=0; i<chunkSize*numDim; i++) tempCoords[i] = -1.0;

            if(!errorFlag){
                chunkNodes=0;
                for(i=0;i<chunkSize;i++){
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

    // fprintf(stderr, "chunkNodes=%d on rank %d\n",chunkNodes, mpi_info->rank);
    // throw FinleyAdapterException("die");

#ifdef ESYS_MPI    
    if(mpi_info->size>1){
        MPI_Bcast(&errorFlag, 1, MPI_INT,  0, mpi_info->comm);
    }
#endif
    if(errorFlag){
        return errorFlag;
    }    

    if (!noError()) return 6;
    mesh_p->Nodes->allocTable(chunkNodes);
    if (!noError()) return 6;

#pragma omp parallel for private (i, j) schedule(static)
    for (i=0; i<chunkNodes; i++) {
        mesh_p->Nodes->Id[i]                     = tempInts[i];
        mesh_p->Nodes->globalDegreesOfFreedom[i] = tempInts[i];
        mesh_p->Nodes->Tag[i]=0;
        // fprintf(stderr,"node id%d %d: ",mesh_p->Nodes->Id[i], mesh_p->Nodes->globalDegreesOfFreedom[i]);
        
        // mesh_p->Nodes->Tag[i]              = tempInts[chunkSize*2+i];
        for (j=0; j<numDim; j++) {
            mesh_p->Nodes->Coordinates[INDEX2(j,i,numDim)]  = tempCoords[i*numDim+j];
            // fprintf(stderr, "%g ", tempCoords[i*numDim+j]);
        }
        // fprintf(stderr, "\n");

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
        /* find line staring with $ */
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
                }
            } while(line[0] != '$');
        }


        if (mpi_info->rank == 0 ) {
            if (!strncmp(&line[1], "MeshFormat", 10)) {
                logicFlag = 1;
            }

            else if ( !strncmp(&line[1], "NOD", 3)   ||
                      !strncmp(&line[1], "NOE", 3)   ||
                      !strncmp(&line[1], "Nodes", 5) ) {
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
        //fprintf(stderr,"broadcasted\n");
#endif
        // MPI_Barrier(mpi_info->comm);
        //fprintf(stderr,"logic flag:%d on rank %d  at line %d\n", logicFlag,mpi_info->rank,__LINE__);
        /* format */
        if (mpi_info->rank==0 && logicFlag ==1) {
                scan_ret = fscanf(fileHandle_p, "%lf %d %d\n", &version, &format, &size);
                FSCANF_CHECK(scan_ret);
                //fprintf(stderr,"mesh format errorFlag:%d on rank %d \n",errorFlag,mpi_info->rank);
        }
        /* nodes are read */
        else if (logicFlag == 2) {

            nodesRead=true;
            errorFlag = getNodes(mpi_info, mesh_p, fileHandle_p, numDim,error_msg); 
            //fprintf(stderr,"nodes errorFlag:%d on rank %d \n",errorFlag,mpi_info->rank);
            logicFlag=0;     
        }
       
        /* elements */
        else if(logicFlag==3) {
            elementsRead=true;
            errorFlag=getElements(mpi_info, mesh_p, fileHandle_p, error_msg, useMacroElements, fname, numDim, version, order, reduced_order);
            //fprintf(stderr,"elements errorFlag:%d on rank %d \n",errorFlag,mpi_info->rank);
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
                //fprintf(stderr,"elements errorFlag:%d on rank %d \n",errorFlag,mpi_info->rank);

            }
            //fprintf(stderr,"physical errorFlag:%d on rank %d \n",errorFlag,mpi_info->rank);
            logicFlag=0;
        }

        //read closing tag $EndElements
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
                    fprintf(stderr,"returning NULL\n");
                    throw FinleyAdapterException("not noError");
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
    if (noError()) mesh_p->resolveNodeIds();
    // rearrange elements
    if (noError()) mesh_p->prepare(optimize);

    if (!noError()) {
        delete mesh_p;
        return NULL;
    }
    return mesh_p;
 }

} // namespace finley

