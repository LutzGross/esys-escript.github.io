
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*   Finley: read mesh from file */

/************************************************************************************/

#include <ctype.h>
#include "Mesh.h"

#define FSCANF_CHECK(scan_ret, reason) { if (scan_ret == EOF) { perror(reason); Finley_setError(IO_ERROR,"scan error while reading finley file"); return NULL;} }


Finley_Mesh* Finley_Mesh_read(char* fname,index_t order, index_t reduced_order,  bool_t optimize)
{
    Esys_MPIInfo *mpi_info = NULL;
    dim_t numNodes, numDim, numEle, i0, i1;
    Finley_Mesh *mesh_p=NULL;
    Finley_ReferenceElementSet *refPoints=NULL, *refContactElements=NULL, *refFaceElements=NULL, *refElements=NULL;
    char name[LenString_MAX],element_type[LenString_MAX],frm[20];
    char error_msg[LenErrorMsg_MAX];
    FILE *fileHandle_p = NULL;
    Finley_ElementTypeId typeID=Finley_NoRef;
    int scan_ret;

    Finley_resetError();
    mpi_info = Esys_MPIInfo_alloc( MPI_COMM_WORLD );

    if (mpi_info->rank == 0) {
        /* get file handle */
        fileHandle_p = fopen(fname, "r");
        if (fileHandle_p==NULL) {
            sprintf(error_msg,"Finley_Mesh_read: Opening file %s for reading failed.",fname);
            Finley_setError(IO_ERROR,error_msg);
            Esys_MPIInfo_free( mpi_info );
            return NULL;
        }

        /* read header */
        sprintf(frm,"%%%d[^\n]",LenString_MAX-1);
        scan_ret = fscanf(fileHandle_p, frm, name);
        FSCANF_CHECK(scan_ret, "Finley_Mesh_read")

        /* get the number of nodes */
        scan_ret = fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
        FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
    }

    #ifdef ESYS_MPI
        /* MPI Broadcast numDim, numNodes, name if there are multiple MPI procs*/
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
            MPI_Bcast (temp1, 3, MPI_INT,  0, mpi_info->comm);
            numDim = temp1[0];
            numNodes = temp1[1];
            MPI_Bcast (name, temp1[2], MPI_CHAR, 0, mpi_info->comm);
        }
    #endif

    /* allocate mesh */
    mesh_p = Finley_Mesh_alloc(name,numDim,mpi_info);
	
    if (Finley_noError()) {
		/* Each CPU will get at most chunkSize nodes so the message has to be sufficiently large */
		int chunkSize = numNodes / mpi_info->size + 1, totalNodes=0, chunkNodes=0,  nextCPU=1;
		int *tempInts = new index_t[chunkSize*3+1];        /* Stores the integer message data */
		double *tempCoords = new double[chunkSize*numDim]; /* Stores the double message data */

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
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					if (2 == numDim) {
						scan_ret = fscanf(fileHandle_p, "%d %d %d %le %le\n",
											&tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
											&tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					if (3 == numDim) {
						scan_ret = fscanf(fileHandle_p, "%d %d %d %le %le %le\n",
											&tempInts[0+i1], &tempInts[chunkSize+i1], &tempInts[chunkSize*2+i1],
											&tempCoords[i1*numDim+0], &tempCoords[i1*numDim+1], &tempCoords[i1*numDim+2]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					totalNodes++; /* When do we quit the infinite loop? */
					chunkNodes++; /* How many nodes do we actually have in this chunk? It may be smaller than chunkSize. */
				}
				if (chunkNodes > chunkSize) {
					Finley_setError(ESYS_MPI_ERROR, "Finley_Mesh_read: error reading chunks of mesh, data too large for message size");
					return NULL;
				}
				#ifdef ESYS_MPI
					/* Eventually we'll send chunkSize nodes to each CPU numbered 1 ... mpi_info->size-1, here goes one of them */
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
		else {  /* Worker */
			#ifdef ESYS_MPI
				/* Each worker receives two messages */
				MPI_Status status;
				MPI_Recv(tempInts, chunkSize*3+1, MPI_INT, 0, 81720, mpi_info->comm, &status);
				MPI_Recv(tempCoords, chunkSize*numDim, MPI_DOUBLE, 0, 81721, mpi_info->comm, &status);
				chunkNodes = tempInts[chunkSize*3];   /* How many nodes are in this workers chunk? */
			#endif
		}   /* Worker */

		/* Copy node data from tempMem to mesh_p */
        mesh_p->Nodes->allocTable(chunkNodes);
        
		if (Finley_noError()) {
			#pragma omp parallel for private (i0, i1) schedule(static)
			for (i0=0; i0<chunkNodes; i0++) {
				mesh_p->Nodes->Id[i0]               = tempInts[0+i0];
				mesh_p->Nodes->globalDegreesOfFreedom[i0]       = tempInts[chunkSize+i0];
				mesh_p->Nodes->Tag[i0]              = tempInts[chunkSize*2+i0];
				for (i1=0; i1<numDim; i1++) {
					mesh_p->Nodes->Coordinates[INDEX2(i1,i0,numDim)]  = tempCoords[i0*numDim+i1];
				}
			}
        }
		delete[] tempInts;
		delete[] tempCoords;
	}

	/* ***********************************  read elements ****************************************************************************************/
	if (Finley_noError()) {

		/* Read the element typeID */
		if (mpi_info->rank == 0) {
			scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
			FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
            typeID=Finley_ReferenceElement_getTypeId(element_type);
		}
		#ifdef ESYS_MPI
			if (mpi_info->size > 1) {
				int temp1[2], mpi_error;
				temp1[0] = (int) typeID;
				temp1[1] = numEle;
				mpi_error = MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
				if (mpi_error != MPI_SUCCESS) {
					Finley_setError(ESYS_MPI_ERROR, "Finley_Mesh_read: broadcast of Element typeID failed");
					return NULL;
				}
				typeID = (Finley_ElementTypeId) temp1[0];
				numEle = temp1[1];
			}
		#endif
        if (typeID==Finley_NoRef) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
          }
    }

    /* Allocate the ElementFile */
	if (Finley_noError()) {
		refElements= Finley_ReferenceElementSet_alloc(typeID,order, reduced_order);
		mesh_p->Elements=Finley_ElementFile_alloc(refElements, mpi_info);
		numNodes = mesh_p->Elements->referenceElementSet->numNodes; /* New meaning for numNodes: num nodes per element */
	}

    /* *************************** Read the element data *****************************************************************************************/
    if (Finley_noError()) {

		int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1, chunkEle=0;
		int *tempInts = new index_t[chunkSize*(2+numNodes)+1]; /* Store Id + Tag + node list (+ one int at end for chunkEle) */
		/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
		if (mpi_info->rank == 0) {  /* Master */
			for (;;) {            /* Infinite loop */
				#pragma omp parallel for private (i0) schedule(static)
				for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
        
				chunkEle = 0;
				for (i0=0; i0<chunkSize; i0++) {
					if (totalEle >= numEle) break; /* End inner loop */
					scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					for (i1 = 0; i1 < numNodes; i1++) {
						scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					scan_ret = fscanf(fileHandle_p, "\n");
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					totalEle++;
					chunkEle++;
				}
				#ifdef ESYS_MPI
					/* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
					if (nextCPU < mpi_info->size) {
						tempInts[chunkSize*(2+numNodes)] = chunkEle;
						MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81722, mpi_info->comm);
					}
				#endif
				nextCPU++;
				/* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
				if (nextCPU > mpi_info->size) break; /* End infinite loop */
			} /* Infinite loop */
		}   /* End master */
		else {  /* Worker */
			#ifdef ESYS_MPI
				/* Each worker receives one message */
				MPI_Status status;
				MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81722, mpi_info->comm, &status);
				chunkEle = tempInts[chunkSize*(2+numNodes)];
			#endif
		}   /* Worker */

		
		Finley_ElementFile_allocTable(mesh_p->Elements, chunkEle);

	
		/* Copy Element data from tempInts to mesh_p */
		if (Finley_noError()) {

			mesh_p->Elements->minColor=0;
			mesh_p->Elements->maxColor=chunkEle-1;
			#pragma omp parallel for private (i0, i1) schedule(static)
			for (i0=0; i0<chunkEle; i0++) {
				mesh_p->Elements->Id[i0]    = tempInts[i0*(2+numNodes)+0];
				mesh_p->Elements->Tag[i0]   = tempInts[i0*(2+numNodes)+1];
				mesh_p->Elements->Owner[i0]  =mpi_info->rank;
				mesh_p->Elements->Color[i0] = i0;
				for (i1 = 0; i1 < numNodes; i1++) {
					mesh_p->Elements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
				}
			}
		}

		delete[] tempInts;
	} 
	/* ******************** end of Read the element data ******************************************************/

    /* ********************* read face elements *************************************************************************************/
    if (Finley_noError()) {
		/* Read the element typeID */

		if (mpi_info->rank == 0) {
       	     scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
			 FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
			 typeID=Finley_ReferenceElement_getTypeId(element_type);
		}
		#ifdef ESYS_MPI
			if (mpi_info->size > 1) {
				int temp1[2];
				temp1[0] = (int) typeID;
				temp1[1] = numEle;
				MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
				typeID = (Finley_ElementTypeId) temp1[0];
				numEle = temp1[1];
			}
		#endif
        if (typeID==Finley_NoRef) {
            sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
        }
		if (Finley_noError()) {
			/* Allocate the ElementFile */
			refFaceElements= Finley_ReferenceElementSet_alloc(typeID,order, reduced_order);
			mesh_p->FaceElements=Finley_ElementFile_alloc(refFaceElements, mpi_info);
			numNodes = mesh_p->FaceElements->referenceElementSet->numNodes; /* New meaning for numNodes: num nodes per element */
		}

	}

    /* ********************** Read the face element data ******************************************************************************* */

	if (Finley_noError()) {
		int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1, chunkEle=0;
    	int *tempInts = new index_t[chunkSize*(2+numNodes)+1]; /* Store Id + Tag + node list (+ one int at end for chunkEle) */
		/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
		if (mpi_info->rank == 0) {  /* Master */
			for (;;) {            /* Infinite loop */
				#pragma omp parallel for private (i0) schedule(static)
				for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
        
				chunkEle = 0;
				for (i0=0; i0<chunkSize; i0++) {
					if (totalEle >= numEle) break; /* End inner loop */
					scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					for (i1 = 0; i1 < numNodes; i1++) {
						scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					scan_ret = fscanf(fileHandle_p, "\n");
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
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
				/* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
				if (nextCPU > mpi_info->size) break; /* End infinite loop */
			} /* Infinite loop */
		}   /* End master */
		else {  /* Worker */
			#ifdef ESYS_MPI
				/* Each worker receives one message */
				MPI_Status status;
				MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81723, mpi_info->comm, &status);
				chunkEle = tempInts[chunkSize*(2+numNodes)];
				#endif
		}   /* Worker */
		
		Finley_ElementFile_allocTable(mesh_p->FaceElements, chunkEle);
		
		if (Finley_noError()) {
			/* Copy Element data from tempInts to mesh_p */
			
        	mesh_p->FaceElements->minColor=0;
			mesh_p->FaceElements->maxColor=chunkEle-1;
			#pragma omp parallel for private (i0, i1)
			for (i0=0; i0<chunkEle; i0++) {
				mesh_p->FaceElements->Id[i0]    = tempInts[i0*(2+numNodes)+0];
				mesh_p->FaceElements->Tag[i0]   = tempInts[i0*(2+numNodes)+1];
				mesh_p->FaceElements->Owner[i0]  =mpi_info->rank;
				mesh_p->FaceElements->Color[i0] = i0;
				for (i1 = 0; i1 < numNodes; i1++) {
					mesh_p->FaceElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
				}
			}
        }

		delete[] tempInts;
	} 
	/* ************************************* end of Read the face element data *************************************** */


	/* ************************************* read contact elements ************************************************** */

    /* Read the element typeID */
	if (Finley_noError()) {
		if (mpi_info->rank == 0) {
			scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
			FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
            typeID=Finley_ReferenceElement_getTypeId(element_type);
		}
		#ifdef ESYS_MPI
	  		if (mpi_info->size > 1) {
				int temp1[2];
				temp1[0] = (int) typeID;
				temp1[1] = numEle;
				MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
				typeID = (Finley_ElementTypeId) temp1[0];
				numEle = temp1[1];
			}
		#endif
        if (typeID==Finley_NoRef) {
			sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
         }
    }

	if (Finley_noError()) {
		/* Allocate the ElementFile */
		refContactElements= Finley_ReferenceElementSet_alloc(typeID,order, reduced_order);
		mesh_p->ContactElements=Finley_ElementFile_alloc(refContactElements, mpi_info);
        numNodes = mesh_p->ContactElements->referenceElementSet->numNodes; /* New meaning for numNodes: num nodes per element */
	}
    /* *************************** Read the contact element data ************************************************* */
    if (Finley_noError()) {
		int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1, chunkEle=0;
		int *tempInts = new index_t[chunkSize*(2+numNodes)+1]; /* Store Id + Tag + node list (+ one int at end for chunkEle) */
		/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
		if (mpi_info->rank == 0) {  /* Master */
			for (;;) {            /* Infinite loop */
				#pragma omp parallel for private (i0) schedule(static)
				for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
       
				chunkEle = 0;
				for (i0=0; i0<chunkSize; i0++) {
					if (totalEle >= numEle) break; /* End inner loop */
					scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					for (i1 = 0; i1 < numNodes; i1++) {
						scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					scan_ret = fscanf(fileHandle_p, "\n");
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					totalEle++;
					chunkEle++;
				}
				#ifdef ESYS_MPI
					/* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
					if (nextCPU < mpi_info->size) {
						tempInts[chunkSize*(2+numNodes)] = chunkEle;
						MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81724, mpi_info->comm);
					}
				#endif
				nextCPU++;
				/* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
				if (nextCPU > mpi_info->size) break; /* End infinite loop */
			} /* Infinite loop */
		}   /* End master */
		else {  /* Worker */
			#ifdef ESYS_MPI
				/* Each worker receives one message */
				MPI_Status status;
				MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81724, mpi_info->comm, &status);
				chunkEle = tempInts[chunkSize*(2+numNodes)]	;
			#endif
		}   /* Worker */

		/* Copy Element data from tempInts to mesh_p */
   		 Finley_ElementFile_allocTable(mesh_p->ContactElements, chunkEle);
		 
		 if (Finley_noError()) {		 
			 mesh_p->ContactElements->minColor=0;
			 mesh_p->ContactElements->maxColor=chunkEle-1;
			 #pragma omp parallel for private (i0, i1)
			 for (i0=0; i0<chunkEle; i0++) {
				 mesh_p->ContactElements->Id[i0] = tempInts[i0*(2+numNodes)+0];
				 mesh_p->ContactElements->Tag[i0]    = tempInts[i0*(2+numNodes)+1];
				 mesh_p->ContactElements->Owner[i0]  =mpi_info->rank;
				 mesh_p->ContactElements->Color[i0] = i0;
				 for (i1 = 0; i1 < numNodes; i1++) {
					 mesh_p->ContactElements->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
				 }
			 }
        }
		delete[] tempInts;
	} /* end of Read the contact element data */

	/* ********************************* read nodal elements ****************************************************** */

    /* *******************************  Read the element typeID */

    if (Finley_noError()) {
		if (mpi_info->rank == 0) {
			scan_ret = fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
			FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
            typeID=Finley_ReferenceElement_getTypeId(element_type);
		}
		#ifdef ESYS_MPI
			if (mpi_info->size > 1) {
				int temp1[2];
				temp1[0] = (int) typeID;
				temp1[1] = numEle;
				MPI_Bcast (temp1, 2, MPI_INT,  0, mpi_info->comm);
				typeID = (Finley_ElementTypeId) temp1[0];
				numEle = temp1[1];
			}
		#endif
        if (typeID==Finley_NoRef) {
			sprintf(error_msg, "Finley_Mesh_read: Unidentified element type %s", element_type);
            Finley_setError(VALUE_ERROR, error_msg);
         }
    }

	if (Finley_noError()) {
		/* Allocate the ElementFile */
		refPoints= Finley_ReferenceElementSet_alloc(typeID,order, reduced_order);
		mesh_p->Points=Finley_ElementFile_alloc(refPoints, mpi_info);
		numNodes = mesh_p->Points->referenceElementSet->numNodes; /* New meaning for numNodes: num nodes per element */
	}

	/**********************************  Read the nodal element data **************************************************/
    if (Finley_noError()) {
		int chunkSize = numEle / mpi_info->size + 1, totalEle=0, nextCPU=1, chunkEle=0;
		int *tempInts = new index_t[chunkSize*(2+numNodes)+1]; /* Store Id + Tag + node list (+ one int at end for chunkEle) */
		/* Elements are specified as a list of integers...only need one message instead of two as with the nodes */
		if (mpi_info->rank == 0) {  /* Master */
			for (;;) {            /* Infinite loop */
				#pragma omp parallel for private (i0) schedule(static)
				for (i0=0; i0<chunkSize*(2+numNodes)+1; i0++) tempInts[i0] = -1;
        
				chunkEle = 0;
				for (i0=0; i0<chunkSize; i0++) {
					if (totalEle >= numEle) break; /* End inner loop */
					scan_ret = fscanf(fileHandle_p, "%d %d", &tempInts[i0*(2+numNodes)+0], &tempInts[i0*(2+numNodes)+1]);
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					for (i1 = 0; i1 < numNodes; i1++) {
						scan_ret = fscanf(fileHandle_p, " %d", &tempInts[i0*(2+numNodes)+2+i1]);
						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					}
					scan_ret = fscanf(fileHandle_p, "\n");
					FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
					totalEle++;
					chunkEle++;
				}
				#ifdef ESYS_MPI
					/* Eventually we'll send chunk of elements to each CPU except 0 itself, here goes one of them */
					if (nextCPU < mpi_info->size) {
						tempInts[chunkSize*(2+numNodes)] = chunkEle;
						MPI_Send(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, nextCPU, 81725, mpi_info->comm);
					}
				#endif
				nextCPU++;
				/* Infinite loop ends when I've read a chunk for each of the worker nodes plus one more chunk for the master */
				if (nextCPU > mpi_info->size) break; /* End infinite loop */
			} /* Infinite loop */
		}   /* End master */
		else {  /* Worker */
			#ifdef ESYS_MPI
				/* Each worker receives one message */
				MPI_Status status;
				MPI_Recv(tempInts, chunkSize*(2+numNodes)+1, MPI_INT, 0, 81725, mpi_info->comm, &status);
				chunkEle = tempInts[chunkSize*(2+numNodes)];
			#endif
		}   /* Worker */

		/* Copy Element data from tempInts to mesh_p */
		Finley_ElementFile_allocTable(mesh_p->Points, chunkEle);
        
		if (Finley_noError()) {
			mesh_p->Points->minColor=0;
			mesh_p->Points->maxColor=chunkEle-1;
			#pragma omp parallel for private (i0, i1) schedule(static)
			for (i0=0; i0<chunkEle; i0++) {
				mesh_p->Points->Id[i0]  = tempInts[i0*(2+numNodes)+0];
				mesh_p->Points->Tag[i0] = tempInts[i0*(2+numNodes)+1];
				mesh_p->Points->Owner[i0]  =mpi_info->rank;
				mesh_p->Points->Color[i0] = i0;
				for (i1 = 0; i1 < numNodes; i1++) {
					mesh_p->Points->Nodes[INDEX2(i1, i0, numNodes)] = tempInts[i0*(2+numNodes)+2+i1];
				}
			}
		}

		delete[] tempInts;
	} /* ******************************** end of Read the nodal element data *********************************************************************************** */

	
	/******************  get the name tags *****************************************/
	if (Finley_noError()) {
        char *remainder=0, *ptr;
        size_t len=0;
		#ifdef ESYS_MPI
        	int len_i;
		#endif
        int tag_key;
        if (mpi_info->rank == 0) {  /* Master */
			/* Read the word 'Tag' */
			if (! feof(fileHandle_p)) {
				scan_ret = fscanf(fileHandle_p, "%s\n", name);
				FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
			}

			#if defined(_WIN32)  /* windows ftell lies on unix formatted text files */
				remainder = NULL;
				len=0;
				while (1)
				{
					size_t malloc_chunk = 1024;
					size_t buff_size = 0;
					int ch;

					ch = fgetc(fileHandle_p);
					if( ch == '\r' )
					{
						continue;
					}
					if( len+1 > buff_size )
					{
						TMPMEMREALLOC(remainder,remainder,buff_size+malloc_chunk,char);
					} 
					if( ch == EOF )
					{
						/* hit EOF */
						remainder[len] = (char)0;
						break;
					}
					remainder[len] = (char)ch;
					len++;
				}
			#else
				/* Read rest of file in one chunk, after using seek to find length */
				{
					long cur_pos, end_pos;

					cur_pos = ftell(fileHandle_p);
					fseek(fileHandle_p, 0L, SEEK_END);
					end_pos = ftell(fileHandle_p);
					fseek(fileHandle_p, (long)cur_pos, SEEK_SET);
					remainder = new char[end_pos-cur_pos+1];
					if (! feof(fileHandle_p))
					{
						scan_ret = fread(remainder, (size_t) end_pos-cur_pos,
                        			         sizeof(char), fileHandle_p);

						FSCANF_CHECK(scan_ret, "Finley_Mesh_read")
						remainder[end_pos-cur_pos] = 0;
					}
				}
			#endif
			len = strlen(remainder);
			// trim the string
			while ((len>1) && isspace(remainder[--len])) {remainder[len]=0;}
			len = strlen(remainder);
			// shrink the allocation unit
			//TMPMEMREALLOC(remainder,remainder,len+1,char);
        } /* Master */
		#ifdef ESYS_MPI

        	len_i=(int) len;
			MPI_Bcast (&len_i, 1, MPI_INT,  0, mpi_info->comm);
			len=(size_t) len_i;
			if (mpi_info->rank != 0) {
				remainder = new char[len+1];
				remainder[0] = 0;
			}
			if (MPI_Bcast (remainder, len+1, MPI_CHAR,  0, mpi_info->comm) !=
				   	MPI_SUCCESS)
				Finley_setError(ESYS_MPI_ERROR, "Finley_Mesh_read: broadcast of remainder failed");
		#endif

		if (remainder[0]) {
			ptr = remainder;
			do {
				sscanf(ptr, "%s %d\n", name, &tag_key);
				if (*name) Finley_Mesh_addTagMap(mesh_p,name,tag_key);
				ptr++;
			} while(NULL != (ptr = strchr(ptr, '\n')) && *ptr);
		}
		if (remainder)
			delete[] remainder;
	}

	/* close file */
	if (mpi_info->rank == 0) fclose(fileHandle_p);

	/*   resolve id's : */
	/* rearrange elements: */
	if (Finley_noError()) Finley_Mesh_resolveNodeIds(mesh_p);
	if (Finley_noError()) Finley_Mesh_prepare(mesh_p, optimize);

	/* that's it */
	if (! Finley_noError()) {
		Finley_Mesh_free(mesh_p);
	}
	/* free up memory */
	Finley_ReferenceElementSet_dealloc(refPoints);
	Finley_ReferenceElementSet_dealloc(refContactElements);
	Finley_ReferenceElementSet_dealloc(refFaceElements);
	Finley_ReferenceElementSet_dealloc(refElements);
	Esys_MPIInfo_free( mpi_info );
	return mesh_p;
}

