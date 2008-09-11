
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: read mesh */

/**************************************************************/

#include <ctype.h>
#include "Mesh.h"

/**************************************************************/

/*  reads a mesh from a Finley file of name fname */

Finley_Mesh* Finley_Mesh_read(char* fname,index_t order, index_t reduced_order,  bool_t optimize) 

{

  Paso_MPIInfo *mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  dim_t numNodes, numDim, numEle, i0, i1;
  index_t tag_key;
  Finley_Mesh *mesh_p=NULL;
  char name[LenString_MAX],element_type[LenString_MAX],frm[20];
  char error_msg[LenErrorMsg_MAX];
  double time0=Finley_timer();
  FILE *fileHandle_p = NULL;
  ElementTypeId typeID, faceTypeID, contactTypeID, pointTypeID;
  
  Finley_resetError();

  if (mpi_info->size > 1) { 
     Finley_setError(SYSTEM_ERROR,"Finley_Mesh_read: MPI is not suporrted yet.");
  } else {
     /* get file handle */
     fileHandle_p = fopen(fname, "r");
     if (fileHandle_p==NULL) {
       sprintf(error_msg,"Finley_Mesh_read: Opening file %s for reading failed.",fname);
       Finley_setError(IO_ERROR,error_msg);
       Paso_MPIInfo_free( mpi_info );
       return NULL;
     }
   
     /* read header */
     sprintf(frm,"%%%d[^\n]",LenString_MAX-1);
     fscanf(fileHandle_p, frm, name);
   
     /* get the nodes */
   
     fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
     /* allocate mesh */
     mesh_p = Finley_Mesh_alloc(name,numDim,order,reduced_order,mpi_info);
     if (Finley_noError()) {
   
        /* read nodes */
        Finley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
        if (Finley_noError()) {
           if (1 == numDim) {
               for (i0 = 0; i0 < numNodes; i0++)
   	            fscanf(fileHandle_p, "%d %d %d %le\n", &mesh_p->Nodes->Id[i0],
   	                   &mesh_p->Nodes->globalDegreesOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
   	                   &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)]);
           } else if (2 == numDim) {
                    for (i0 = 0; i0 < numNodes; i0++)
   	                      fscanf(fileHandle_p, "%d %d %d %le %le\n", &mesh_p->Nodes->Id[i0],
   	                             &mesh_p->Nodes->globalDegreesOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
   	                             &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
   	                             &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)]);
           } else if (3 == numDim) {
                    for (i0 = 0; i0 < numNodes; i0++)
   	                      fscanf(fileHandle_p, "%d %d %d %le %le %le\n", &mesh_p->Nodes->Id[i0],
   	                             &mesh_p->Nodes->globalDegreesOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
   	                             &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
   	                             &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
   	                             &mesh_p->Nodes->Coordinates[INDEX2(2,i0,numDim)]);
           } /* if else else */
        }
        /* read elements */
        if (Finley_noError()) {
   
           fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
           typeID=Finley_RefElement_getTypeId(element_type);
           if (typeID==NoType) {
             sprintf(error_msg,"Finley_Mesh_read :Unidentified element type %s",element_type);
             Finley_setError(VALUE_ERROR,error_msg);
           } else {
             /* read the elements */
             mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
             if (Finley_noError()) {
                 Finley_ElementFile_allocTable(mesh_p->Elements, numEle);
                 mesh_p->Elements->minColor=0;
                 mesh_p->Elements->maxColor=numEle-1;
                 if (Finley_noError()) {
                    for (i0 = 0; i0 < numEle; i0++) {
                      fscanf(fileHandle_p, "%d %d", &mesh_p->Elements->Id[i0], &mesh_p->Elements->Tag[i0]);
                      mesh_p->Elements->Color[i0]=i0;
                      mesh_p->Elements->Owner[i0]=0;
                      for (i1 = 0; i1 < mesh_p->Elements->ReferenceElement->Type->numNodes; i1++) {
                           fscanf(fileHandle_p, " %d",
                              &mesh_p->Elements->Nodes[INDEX2(i1, i0, mesh_p->Elements->ReferenceElement->Type->numNodes)]);
                      }	/* for i1 */
                      fscanf(fileHandle_p, "\n");
                    } /* for i0 */
                 }
             }
          }
        }
        /* get the face elements */
        if (Finley_noError()) {
             fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
             faceTypeID=Finley_RefElement_getTypeId(element_type);
             if (faceTypeID==NoType) {
               sprintf(error_msg,"Finley_Mesh_read :Unidentified element type %s for face elements",element_type);
               Finley_setError(VALUE_ERROR,error_msg);
             } else {
                mesh_p->FaceElements=Finley_ElementFile_alloc(faceTypeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
                if (Finley_noError()) {
                   Finley_ElementFile_allocTable(mesh_p->FaceElements, numEle);
                   if (Finley_noError()) {
                      mesh_p->FaceElements->minColor=0;
                      mesh_p->FaceElements->maxColor=numEle-1;
                      for (i0 = 0; i0 < numEle; i0++) {
                        fscanf(fileHandle_p, "%d %d", &mesh_p->FaceElements->Id[i0], &mesh_p->FaceElements->Tag[i0]);
                        mesh_p->FaceElements->Color[i0]=i0;
                        mesh_p->FaceElements->Owner[i0]=0;
                        for (i1 = 0; i1 < mesh_p->FaceElements->ReferenceElement->Type->numNodes; i1++) {
                             fscanf(fileHandle_p, " %d",
                                &mesh_p->FaceElements->Nodes[INDEX2(i1, i0, mesh_p->FaceElements->ReferenceElement->Type->numNodes)]);
                        }	/* for i1 */
                        fscanf(fileHandle_p, "\n");
                      } /* for i0 */
                   }
                }
             }
        }
        /* get the Contact face element */
        if (Finley_noError()) {
             fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
             contactTypeID=Finley_RefElement_getTypeId(element_type);
             if (contactTypeID==NoType) {
               sprintf(error_msg,"Finley_Mesh_read: Unidentified element type %s for contact elements",element_type);
               Finley_setError(VALUE_ERROR,error_msg);
             } else {
               mesh_p->ContactElements=Finley_ElementFile_alloc(contactTypeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
               if (Finley_noError()) {
                   Finley_ElementFile_allocTable(mesh_p->ContactElements, numEle);
                   if (Finley_noError()) {
                      mesh_p->ContactElements->minColor=0;
                      mesh_p->ContactElements->maxColor=numEle-1;
                      for (i0 = 0; i0 < numEle; i0++) {
                        fscanf(fileHandle_p, "%d %d", &mesh_p->ContactElements->Id[i0], &mesh_p->ContactElements->Tag[i0]);
                        mesh_p->ContactElements->Color[i0]=i0;
                        mesh_p->ContactElements->Owner[i0]=0;
                        for (i1 = 0; i1 < mesh_p->ContactElements->ReferenceElement->Type->numNodes; i1++) {
                            fscanf(fileHandle_p, " %d",
                               &mesh_p->ContactElements->Nodes[INDEX2(i1, i0, mesh_p->ContactElements->ReferenceElement->Type->numNodes)]);
                        }	/* for i1 */
                        fscanf(fileHandle_p, "\n");
                      } /* for i0 */
                  }
               }
             }
        }  
        /* get the nodal element */
        if (Finley_noError()) {
             fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
             pointTypeID=Finley_RefElement_getTypeId(element_type);
             if (pointTypeID==NoType) {
               sprintf(error_msg,"Finley_Mesh_read: Unidentified element type %s for points",element_type);
               Finley_setError(VALUE_ERROR,error_msg);
             }
             mesh_p->Points=Finley_ElementFile_alloc(pointTypeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
             if (Finley_noError()) {
                Finley_ElementFile_allocTable(mesh_p->Points, numEle);
                if (Finley_noError()) {
                   mesh_p->Points->minColor=0;
                   mesh_p->Points->maxColor=numEle-1;
                   for (i0 = 0; i0 < numEle; i0++) {
                     fscanf(fileHandle_p, "%d %d", &mesh_p->Points->Id[i0], &mesh_p->Points->Tag[i0]);
                     mesh_p->Points->Color[i0]=i0;
                     mesh_p->Points->Owner[i0]=0;
                     for (i1 = 0; i1 < mesh_p->Points->ReferenceElement->Type->numNodes; i1++) {
                         fscanf(fileHandle_p, " %d",
                            &mesh_p->Points->Nodes[INDEX2(i1, i0, mesh_p->Points->ReferenceElement->Type->numNodes)]);
                     }	/* for i1 */
                     fscanf(fileHandle_p, "\n");
                   } /* for i0 */
                }
             }
        }
        /* get the name tags */
        if (Finley_noError()) {
           if (feof(fileHandle_p) == 0) {
              fscanf(fileHandle_p, "%s\n", name);
              while (feof(fileHandle_p) == 0) {
                   fscanf(fileHandle_p, "%s %d\n", name, &tag_key);
                   Finley_Mesh_addTagMap(mesh_p,name,tag_key);
              }
           }
        }
     }
     /* close file */
     fclose(fileHandle_p);
   
     /*   resolve id's : */
     /* rearrange elements: */
   
     if (Finley_noError()) Finley_Mesh_resolveNodeIds(mesh_p);
     if (Finley_noError()) Finley_Mesh_prepare(mesh_p, optimize);
   
     /* that's it */
     #ifdef Finley_TRACE
     printf("timing: reading mesh: %.4e sec\n",Finley_timer()-time0);
     #endif
  }

  /* that's it */
  if (! Finley_noError()) {
       Finley_Mesh_free(mesh_p);
  }
  Paso_MPIInfo_free( mpi_info );
  return mesh_p;
}

Finley_Mesh* Finley_Mesh_read_MPI(char* fname,index_t order, index_t reduced_order,  bool_t optimize)

{

  Paso_MPIInfo *mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  dim_t numNodes, numDim, numEle, i0, i1;
  Finley_Mesh *mesh_p=NULL;
  char name[LenString_MAX],element_type[LenString_MAX],frm[20];
  char error_msg[LenErrorMsg_MAX];
  double time0=Finley_timer();
  FILE *fileHandle_p = NULL;
  ElementTypeId typeID, faceTypeID, contactTypeID, pointTypeID;
  Finley_TagMap* tag_map;
  index_t tag_key;

  Finley_resetError();

  if (mpi_info->rank == 0) {
     /* get file handle */
     fileHandle_p = fopen(fname, "r");
     if (fileHandle_p==NULL) {
       sprintf(error_msg,"Finley_Mesh_read: Opening file %s for reading failed.",fname);
       Finley_setError(IO_ERROR,error_msg);
       Paso_MPIInfo_free( mpi_info );
       return NULL;
     }

     /* read header */
     sprintf(frm,"%%%d[^\n]",LenString_MAX-1);
     fscanf(fileHandle_p, frm, name);

     /* get the number of nodes */
     fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
  }

#ifdef PASO_MPI
  /* MPI Broadcast numDim, numNodes, name if there are multiple MPI procs*/
  if (mpi_info->size > 1) {
    int temp1[3], error_code;
    temp1[0] = numDim;
    temp1[1] = numNodes;
    temp1[2] = strlen(name) + 1;
    error_code = MPI_Bcast (temp1, 3, MPI_INT,  0, mpi_info->comm);
    if (error_code != MPI_SUCCESS) {
      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of temp1 failed");
      return NULL;
    }
    numDim = temp1[0];
    numNodes = temp1[1];
    error_code = MPI_Bcast (name, temp1[2], MPI_CHAR, 0, mpi_info->comm);
    if (error_code != MPI_SUCCESS) {
      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of name failed");
      return NULL;
    }
  }
#endif

     /* allocate mesh */
     mesh_p = Finley_Mesh_alloc(name,numDim,order,reduced_order,mpi_info);
     if (Finley_noError()) {
	/* Each CPU will get at most chunkSize nodes so the message has to be sufficiently large */
	int chunkSize = numNodes / mpi_info->size + 1, totalNodes=0, chunkNodes=0, chunkEle=0, nextCPU=1;
	int *tempInts = TMPMEMALLOC(chunkSize*3+1, index_t);		/* Stores the integer message data */
	double *tempCoords = TMPMEMALLOC(chunkSize*numDim, double);	/* Stores the double message data */

	/*
	  Read chunkSize nodes, send it in a chunk to worker CPU which copies chunk into its local mesh_p
	  It doesn't matter that a CPU has the wrong nodes for its elements, this is sorted out later
	  First chunk sent to CPU 1, second to CPU 2, ...
	  Last chunk stays on CPU 0 (the master)
	  The three columns of integers (Id, gDOF, Tag) are gathered into a single array tempInts and sent together in a single MPI message
	*/

	if (mpi_info->rank == 0) {	/* Master */
	  for (;;) {			/* Infinite loop */
	    for (i0=0; i0<chunkSize*3+1; i0++) tempInts[i0] = -1;
	    for (i0=0; i0<chunkSize*numDim; i0++) tempCoords[i0] = -1.0;
	    chunkNodes = 0;
	    for (i1=0; i1<chunkSize; i1++) {
	      if (totalNodes >= numNodes) break;	/* End of inner loop */
              if (1 == numDim)
		fscanf(fileHandle_p, "%d %d %d %le\n",
		  &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
		  &tempCoords[i1*numDim+0]);
              if (2 == numDim)
		fscanf(fileHandle_p, "%d %d %d %le %le\n",
		  &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
		  &tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1]);
              if (3 == numDim)
		fscanf(fileHandle_p, "%d %d %d %le %le %le\n",
		  &tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
		  &tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1], &tempCoords[i1*numDim+2]);
	      totalNodes++; /* When do we quit the infinite loop? */
	      chunkNodes++; /* How many nodes do we actually have in this chunk? It may be smaller than chunkSize. */
	    }
	    if (chunkNodes > chunkSize) {
              Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: error reading chunks of mesh, data too large for message size");
              return NULL;
	    }
#ifdef PASO_MPI
	    /* Eventually we'll send chunkSize nodes to each CPU numbered 1 ... mpi_info->size-1, here goes one of them */
	    if (nextCPU < mpi_info->size) {
              int mpi_error;
	      tempInts[chunkSize*3] = chunkNodes;	/* The message has one more int to send chunkNodes */
	      mpi_error = MPI_Send(tempInts, chunkSize*3+1, MPI_INT, nextCPU, 81720, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempInts failed");
                return NULL;
	      }
	      mpi_error = MPI_Send(tempCoords, chunkSize*numDim, MPI_DOUBLE, nextCPU, 81721, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempCoords failed");
                return NULL;
	      }
	    }
#endif
	    nextCPU++;
	    /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
	    if (nextCPU > mpi_info->size) break; /* End infinite loop */
	  }	/* Infinite loop */
	}	/* End master */
	else {	/* Worker */
#ifdef PASO_MPI
	  /* Each worker receives two messages */
	  MPI_Status status;
          int mpi_error;
	  mpi_error = MPI_Recv(tempInts, chunkSize*3+1, MPI_INT, 0, 81720, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempInts failed");
            return NULL;
	  }
	  mpi_error = MPI_Recv(tempCoords, chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempCoords failed");
            return NULL;
	  }
	  chunkNodes = tempInts[chunkSize*3];	/* How many nodes are in this workers chunk? */
#endif
	}	/* Worker */

	/* Copy node data from tempMem to mesh_p */
        Finley_NodeFile_allocTable(mesh_p->Nodes, chunkNodes);
        if (Finley_noError()) {
	  for (i0=0; i0<chunkNodes; i0++) {
	    mesh_p->Nodes->Id[i0]				= tempInts[0+i0];
	    mesh_p->Nodes->globalDegreesOfFreedom[i0]		= tempInts[chunkSize+i0];
	    mesh_p->Nodes->Tag[i0]				= tempInts[chunkSize*2+i0];
	    for (i1=0; i1<numDim; i1++) {
	      mesh_p->Nodes->Coordinates[INDEX2(i1,i0,numDim)]	= tempCoords[i0*numDim+i1];
	    }
          }
        }

	TMPMEMFREE(tempInts);
	TMPMEMFREE(tempCoords);

        /* read elements */

	/* Read the element typeID */
        if (Finley_noError()) {
	  if (mpi_info->rank == 0) {
            fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
            typeID=Finley_RefElement_getTypeId(element_type);
	  }
#ifdef PASO_MPI
	  if (mpi_info->size > 1) {
	    int temp1[2], mpi_error;
	    temp1[0] = (int) typeID;
	    temp1[1] = numEle;
	    mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
	    if (mpi_error != MPI_SUCCESS) {
	      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of Element typeID failed");
	      return NULL;
	    }
	    typeID = (ElementTypeId) temp1[0];
	    numEle = temp1[1];
	  }
#endif
          if (typeID==NoType) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
          }
	}

      /* Allocate the ElementFile */
      mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
      numNodes = mesh_p->Elements->ReferenceElement->Type->numNodes; /* New meaning for numNodes: num nodes per element */

      /* Read the element data */
      if (Finley_noError()) {
	int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1;
	int *tempInts = TMPMEMALLOC(chunkSize*(2+numNodes)+1, index_t); /* Store Id + Tag + node list (+ one int at end for chunkEle) */
	/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
	if (mpi_info->rank == 0) {	/* Master */
	  for (;;) {			/* Infinite loop */
	    for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
	    chunkEle = 0;
	    for (i0=0; i0<chunkSize; i0++) {
	      if (totalEle >= numEle) break; /* End inner loop */
	      fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
	      for (i1 = 0; i1 < numNodes; i1++) fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
	      fscanf(fileHandle_p, "\n");
	      totalEle++;
	      chunkEle++;
	    }
#ifdef PASO_MPI
	    /* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
	    if (nextCPU < mpi_info->size) {
              int mpi_error;
	      tempInts[chunkSize*(2+numNodes)] = chunkEle;
	      mpi_error = MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81722, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempInts for Elements failed");
                return NULL;
	      }
	    }
#endif
	    nextCPU++;
	    /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
	    if (nextCPU > mpi_info->size) break; /* End infinite loop */
	  }	/* Infinite loop */
	}	/* End master */
	else {	/* Worker */
#ifdef PASO_MPI
	  /* Each worker receives one message */
	  MPI_Status status;
	  int mpi_error;
	  mpi_error = MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81722, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempInts for Elements failed");
            return NULL;
	  }
	  chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
	}	/* Worker */

	/* Copy Element data from tempInts to mesh_p */
	Finley_ElementFile_allocTable(mesh_p->Elements, chunkEle);
        mesh_p->Elements->minColor=0;
        mesh_p->Elements->maxColor=chunkEle-1;
        if (Finley_noError()) {
          #pragma omp parallel for private (i0, i1)
	  for (i0=0; i0<chunkEle; i0++) {
	    mesh_p->Elements->Id[i0]	= tempInts[i0*(2+numNodes)+0];
	    mesh_p->Elements->Tag[i0]	= tempInts[i0*(2+numNodes)+1];
            mesh_p->Elements->Owner[i0]  =mpi_info->rank;
            mesh_p->Elements->Color[i0] = i0;
	    for (i1 = 0; i1 < numNodes; i1++) {
	      mesh_p->Elements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
	    }
          }
        }

	TMPMEMFREE(tempInts);
      } /* end of Read the element data */

        /* read face elements */

	/* Read the element typeID */
        if (Finley_noError()) {
	  if (mpi_info->rank == 0) {
            fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
            typeID=Finley_RefElement_getTypeId(element_type);
	  }
#ifdef PASO_MPI
	  if (mpi_info->size > 1) {
	    int temp1[2], mpi_error;
	    temp1[0] = (int) typeID;
	    temp1[1] = numEle;
	    mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
	    if (mpi_error != MPI_SUCCESS) {
	      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of Element typeID failed");
	      return NULL;
	    }
	    typeID = (ElementTypeId) temp1[0];
	    numEle = temp1[1];
	  }
#endif
          if (typeID==NoType) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
          }
	}

      /* Allocate the ElementFile */
      mesh_p->FaceElements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
      numNodes = mesh_p->FaceElements->ReferenceElement->Type->numNodes; /* New meaning for numNodes: num nodes per element */

      /* Read the face element data */
      if (Finley_noError()) {
	int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1;
	int *tempInts = TMPMEMALLOC(chunkSize*(2+numNodes)+1, index_t); /* Store Id + Tag + node list (+ one int at end for chunkEle) */
	/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
	if (mpi_info->rank == 0) {	/* Master */
	  for (;;) {			/* Infinite loop */
	    for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
	    chunkEle = 0;
	    for (i0=0; i0<chunkSize; i0++) {
	      if (totalEle >= numEle) break; /* End inner loop */
	      fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
	      for (i1 = 0; i1 < numNodes; i1++) fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
	      fscanf(fileHandle_p, "\n");
	      totalEle++;
	      chunkEle++;
	    }
#ifdef PASO_MPI
	    /* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
	    if (nextCPU < mpi_info->size) {
              int mpi_error;
	      tempInts[chunkSize*(2+numNodes)] = chunkEle;
	      mpi_error = MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81723, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempInts for FaceElements failed");
                return NULL;
	      }
	    }
#endif
	    nextCPU++;
	    /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
	    if (nextCPU > mpi_info->size) break; /* End infinite loop */
	  }	/* Infinite loop */
	}	/* End master */
	else {	/* Worker */
#ifdef PASO_MPI
	  /* Each worker receives one message */
	  MPI_Status status;
	  int mpi_error;
	  mpi_error = MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81723, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempInts for FaceElements failed");
            return NULL;
	  }
	  chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
	}	/* Worker */

	/* Copy Element data from tempInts to mesh_p */
	Finley_ElementFile_allocTable(mesh_p->FaceElements, chunkEle);
        mesh_p->FaceElements->minColor=0;
        mesh_p->FaceElements->maxColor=chunkEle-1;
        if (Finley_noError()) {
          #pragma omp parallel for private (i0, i1)
	  for (i0=0; i0<chunkEle; i0++) {
	    mesh_p->FaceElements->Id[i0]	= tempInts[i0*(2+numNodes)+0];
	    mesh_p->FaceElements->Tag[i0]	= tempInts[i0*(2+numNodes)+1];
            mesh_p->FaceElements->Owner[i0]  =mpi_info->rank;
            mesh_p->FaceElements->Color[i0] = i0;
	    for (i1 = 0; i1 < numNodes; i1++) {
	      mesh_p->FaceElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
	    }
          }
        }

	TMPMEMFREE(tempInts);
      } /* end of Read the face element data */

        /* read contact elements */

	/* Read the element typeID */
        if (Finley_noError()) {
	  if (mpi_info->rank == 0) {
            fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
            typeID=Finley_RefElement_getTypeId(element_type);
	  }
#ifdef PASO_MPI
	  if (mpi_info->size > 1) {
	    int temp1[2], mpi_error;
	    temp1[0] = (int) typeID;
	    temp1[1] = numEle;
	    mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
	    if (mpi_error != MPI_SUCCESS) {
	      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of Element typeID failed");
	      return NULL;
	    }
	    typeID = (ElementTypeId) temp1[0];
	    numEle = temp1[1];
	  }
#endif
          if (typeID==NoType) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
          }
	}

      /* Allocate the ElementFile */
      mesh_p->ContactElements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
      numNodes = mesh_p->ContactElements->ReferenceElement->Type->numNodes; /* New meaning for numNodes: num nodes per element */

      /* Read the contact element data */
      if (Finley_noError()) {
	int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1;
	int *tempInts = TMPMEMALLOC(chunkSize*(2+numNodes)+1, index_t); /* Store Id + Tag + node list (+ one int at end for chunkEle) */
	/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
	if (mpi_info->rank == 0) {	/* Master */
	  for (;;) {			/* Infinite loop */
	    for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
	    chunkEle = 0;
	    for (i0=0; i0<chunkSize; i0++) {
	      if (totalEle >= numEle) break; /* End inner loop */
	      fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
	      for (i1 = 0; i1 < numNodes; i1++) fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
	      fscanf(fileHandle_p, "\n");
	      totalEle++;
	      chunkEle++;
	    }
#ifdef PASO_MPI
	    /* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
	    if (nextCPU < mpi_info->size) {
              int mpi_error;
	      tempInts[chunkSize*(2+numNodes)] = chunkEle;
	      mpi_error = MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81724, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempInts for ContactElements failed");
                return NULL;
	      }
	    }
#endif
	    nextCPU++;
	    /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
	    if (nextCPU > mpi_info->size) break; /* End infinite loop */
	  }	/* Infinite loop */
	}	/* End master */
	else {	/* Worker */
#ifdef PASO_MPI
	  /* Each worker receives one message */
	  MPI_Status status;
	  int mpi_error;
	  mpi_error = MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81724, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempInts for ContactElements failed");
            return NULL;
	  }
	  chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
	}	/* Worker */

	/* Copy Element data from tempInts to mesh_p */
	Finley_ElementFile_allocTable(mesh_p->ContactElements, chunkEle);
        mesh_p->ContactElements->minColor=0;
        mesh_p->ContactElements->maxColor=chunkEle-1;
        if (Finley_noError()) {
          #pragma omp parallel for private (i0, i1)
	  for (i0=0; i0<chunkEle; i0++) {
	    mesh_p->ContactElements->Id[i0]	= tempInts[i0*(2+numNodes)+0];
	    mesh_p->ContactElements->Tag[i0]	= tempInts[i0*(2+numNodes)+1];
            mesh_p->ContactElements->Owner[i0]  =mpi_info->rank;
            mesh_p->ContactElements->Color[i0] = i0;
	    for (i1 = 0; i1 < numNodes; i1++) {
	      mesh_p->ContactElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
	    }
          }
        }

	TMPMEMFREE(tempInts);
      } /* end of Read the contact element data */

        /* read nodal elements */

	/* Read the element typeID */
        if (Finley_noError()) {
	  if (mpi_info->rank == 0) {
            fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
            typeID=Finley_RefElement_getTypeId(element_type);
	  }
#ifdef PASO_MPI
	  if (mpi_info->size > 1) {
	    int temp1[2], mpi_error;
	    temp1[0] = (int) typeID;
	    temp1[1] = numEle;
	    mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
	    if (mpi_error != MPI_SUCCESS) {
	      Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of Element typeID failed");
	      return NULL;
	    }
	    typeID = (ElementTypeId) temp1[0];
	    numEle = temp1[1];
	  }
#endif
          if (typeID==NoType) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
          }
	}

      /* Allocate the ElementFile */
      mesh_p->Points=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, mpi_info);
      numNodes = mesh_p->Points->ReferenceElement->Type->numNodes; /* New meaning for numNodes: num nodes per element */

      /* Read the nodal element data */
      if (Finley_noError()) {
	int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1;
	int *tempInts = TMPMEMALLOC(chunkSize*(2+numNodes)+1, index_t); /* Store Id + Tag + node list (+ one int at end for chunkEle) */
	/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
	if (mpi_info->rank == 0) {	/* Master */
	  for (;;) {			/* Infinite loop */
	    for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
	    chunkEle = 0;
	    for (i0=0; i0<chunkSize; i0++) {
	      if (totalEle >= numEle) break; /* End inner loop */
	      fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
	      for (i1 = 0; i1 < numNodes; i1++) fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
	      fscanf(fileHandle_p, "\n");
	      totalEle++;
	      chunkEle++;
	    }
#ifdef PASO_MPI
	    /* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
	    if (nextCPU < mpi_info->size) {
              int mpi_error;
	      tempInts[chunkSize*(2+numNodes)] = chunkEle;
	      mpi_error = MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81725, mpi_info->comm);
	      if ( mpi_error != MPI_SUCCESS ) {
                Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: send of tempInts for Points failed");
                return NULL;
	      }
	    }
#endif
	    nextCPU++;
	    /* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
	    if (nextCPU > mpi_info->size) break; /* End infinite loop */
	  }	/* Infinite loop */
	}	/* End master */
	else {	/* Worker */
#ifdef PASO_MPI
	  /* Each worker receives one message */
	  MPI_Status status;
	  int mpi_error;
	  mpi_error = MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81725, mpi_info->comm, &status);
	  if ( mpi_error != MPI_SUCCESS ) {
            Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: receive of tempInts for Points failed");
            return NULL;
	  }
	  chunkEle = tempInts[chunkSize*(2+numNodes)];
#endif
	}	/* Worker */

	/* Copy Element data from tempInts to mesh_p */
	Finley_ElementFile_allocTable(mesh_p->Points, chunkEle);
        mesh_p->Points->minColor=0;
        mesh_p->Points->maxColor=chunkEle-1;
        if (Finley_noError()) {
          #pragma omp parallel for private (i0, i1)
	  for (i0=0; i0<chunkEle; i0++) {
	    mesh_p->Points->Id[i0]	= tempInts[i0*(2+numNodes)+0];
	    mesh_p->Points->Tag[i0]	= tempInts[i0*(2+numNodes)+1];
            mesh_p->Points->Owner[i0]  =mpi_info->rank;
            mesh_p->Points->Color[i0] = i0;
	    for (i1 = 0; i1 < numNodes; i1++) {
	      mesh_p->Points->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
	    }
          }
        }

	TMPMEMFREE(tempInts);
      } /* end of Read the nodal element data */

      /* get the name tags */
      if (Finley_noError()) {
        char *remainder, *ptr;
        int tag_key, len, error_code;
        long cur_pos, end_pos;
        if (mpi_info->rank == 0) {	/* Master */
	  /* Read the word 'Tag' */
	  if (! feof(fileHandle_p)) fscanf(fileHandle_p, "%s\n", name);
	  /* Read rest of file in one chunk, after using seek to find length */
          cur_pos = ftell(fileHandle_p);
          fseek(fileHandle_p, 0L, SEEK_END);
          end_pos = ftell(fileHandle_p);
          fseek(fileHandle_p, (long)cur_pos, SEEK_SET);
	  remainder = TMPMEMALLOC(end_pos-cur_pos+1, char);
	  if (! feof(fileHandle_p)) fread(remainder, (size_t) end_pos-cur_pos, sizeof(char), fileHandle_p);
	  remainder[end_pos-cur_pos] = 0;
	  len = strlen(remainder);    
	  while ((--len)>0 && isspace(remainder[len])) remainder[len]=0;
	  len = strlen(remainder);
        }
#ifdef PASO_MPI
        error_code = MPI_Bcast (&len, 1, MPI_INT,  0, mpi_info->comm);
        if (error_code != MPI_SUCCESS) {
          Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of tag len failed");
          return NULL;
        }
	if (mpi_info->rank != 0) {
	  remainder = TMPMEMALLOC(len+1, char);
	  remainder[0] = 0;
	}
        error_code = MPI_Bcast (remainder, len+1, MPI_CHAR,  0, mpi_info->comm);
        if (error_code != MPI_SUCCESS) {
          Finley_setError(PASO_MPI_ERROR, "Finley_Mesh_read: broadcast of tags failed");
          return NULL;
        }
#endif
	if (remainder[0]) {
          ptr = remainder;
          do {
            sscanf(ptr, "%s %d\n", name, &tag_key);
            if (*name) Finley_Mesh_addTagMap(mesh_p,name,tag_key);
            ptr++;
          } while(NULL != (ptr = strchr(ptr, '\n')) && *ptr);
          TMPMEMFREE(remainder);
	}
      }

     }

     /* close file */
     if (mpi_info->rank == 0) fclose(fileHandle_p);

     /*   resolve id's : */
     /* rearrange elements: */

     if (Finley_noError()) Finley_Mesh_resolveNodeIds(mesh_p);
     if (Finley_noError()) Finley_Mesh_prepare(mesh_p, optimize);

     /* that's it */
     #ifdef Finley_TRACE
     printf("timing: reading mesh: %.4e sec\n",Finley_timer()-time0);
     #endif

  /* that's it */
  if (! Finley_noError()) {
       Finley_Mesh_free(mesh_p);
  }
  Paso_MPIInfo_free( mpi_info );
  return mesh_p;
}
