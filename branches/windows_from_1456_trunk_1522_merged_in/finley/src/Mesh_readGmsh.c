
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

#include "Mesh.h"
#include <stdio.h>


/**************************************************************/

/*  reads a mesh from a Finley file of name fname */

#define MAX_numNodes_gmsh 10

Finley_Mesh* Finley_Mesh_readGmsh(char* fname ,index_t numDim, index_t order, index_t reduced_order, bool_t optimize) {

  double version = 1.0;
  int format = 0, size = sizeof(double);
  dim_t numNodes, totalNumElements=0, numTags=0, numNodesPerElement, numNodesPerElement2, element_dim;
  index_t e, i0, j, gmsh_type, partition_id, itmp, elementary_id;
  index_t numElements=0, numFaceElements=0, *id=NULL, *tag=NULL, *vertices=NULL;
  Finley_Mesh *mesh_p=NULL;
  char line[LenString_MAX+1];
  char error_msg[LenErrorMsg_MAX];
  double rtmp0, rtmp1;
  double time0=Finley_timer();
  FILE * fileHandle_p = NULL;
  ElementTypeId* element_type=NULL;

  Paso_MPIInfo *mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  Finley_resetError();
  if (mpi_info->size>1) {
    Finley_setError(IO_ERROR,"reading GMSH with MPI is not supported yet.");
    Paso_MPIInfo_free( mpi_info );
    return NULL;
  } else {

     /* allocate mesh */
   
     mesh_p = Finley_Mesh_alloc(fname,numDim,order, reduced_order, mpi_info);
     if (! Finley_noError()) return NULL;
   
     /* get file handle */
     fileHandle_p = fopen(fname, "r");
     if (fileHandle_p==NULL) {
       sprintf(error_msg,"Opening Gmsh file %s for reading failed.",fname);
       Finley_setError(IO_ERROR,error_msg);
       Paso_MPIInfo_free( mpi_info );
       return NULL;
     }
   
     /* start reading */
     while(1) {
       if (! Finley_noError()) break;
       /* find line staring with $ */
       do {
         if( ! fgets(line, sizeof(line), fileHandle_p) ) break;
         if(feof(fileHandle_p)) break;
       } while(line[0] != '$');
   
       if (feof(fileHandle_p)) break;
   
   
       /* format */
       if (!strncmp(&line[1], "MeshFormat", 10)) {
         fscanf(fileHandle_p, "%lf %d %d\n", &version, &format, &size);
       }
       /* nodes are read */
       if ( !strncmp(&line[1], "NOD", 3)   ||
            !strncmp(&line[1], "NOE", 3)   ||
            !strncmp(&line[1], "Nodes", 5)    ) {
        
         fscanf(fileHandle_p, "%d", &numNodes);
         if (! Finley_noError()) break;
         Finley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
         if (! Finley_noError()) break;
         for (i0 = 0; i0 < numNodes; i0++) {
            if (1 == numDim) {
   	       fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
   	              &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)], 
                         &rtmp0, 
                         &rtmp1);
            } else if (2 == numDim) {
   	       fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
   	              &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
   	              &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
                         &rtmp0);
            } else if (3 == numDim) {
   	       fscanf(fileHandle_p, "%d %le %le %le\n", &mesh_p->Nodes->Id[i0],
   	              &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
   	              &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
   	              &mesh_p->Nodes->Coordinates[INDEX2(2,i0,numDim)]);
      
          
            }
            mesh_p->Nodes->globalDegreesOfFreedom[i0]=mesh_p->Nodes->Id[i0];
            mesh_p->Nodes->Tag[i0]=0;
         }
       }
       /* elements */
       else if(!strncmp(&line[1], "ELM", 3) ||
   	    !strncmp(&line[1], "Elements", 8)) {
   
         ElementTypeId final_element_type = NoType;
         ElementTypeId final_face_element_type = NoType;
         numElements=0;
         numFaceElements=0;
         fscanf(fileHandle_p, "%d", &totalNumElements);
   
         id=TMPMEMALLOC(totalNumElements,index_t);
         tag=TMPMEMALLOC(totalNumElements,index_t);
   
   
         element_type=TMPMEMALLOC(totalNumElements,ElementTypeId);
         vertices=TMPMEMALLOC(totalNumElements*MAX_numNodes_gmsh,index_t);
         if (! (Finley_checkPtr(id) || Finley_checkPtr(tag) || Finley_checkPtr(element_type) || Finley_checkPtr(vertices) ) ) {
            /* read all in */
            for(e = 0; e < totalNumElements; e++) {
              fscanf(fileHandle_p, "%d %d", &id[e], &gmsh_type);
              switch (gmsh_type) {
                  case 1:  /* line order 1 */
                      element_type[e]=Line2;
                      element_dim=1;
                      numNodesPerElement=2;  
                      break;
                  case 2:  /* traingle order 1 */
                      element_type[e]=Tri3;
                      numNodesPerElement= 3;
                      element_dim=2;
                      break;
                  case 3:  /* quadrilateral order 1 */
                      element_type[e]=Rec4;
                      numNodesPerElement= 4;
                      element_dim=2;
                      break;
                  case 4:  /* tetrahedron order 1 */
                      element_type[e]=Tet4;
                      numNodesPerElement= 4;
                      element_dim=3;
                      break;
                  case 5:  /* hexahedron order 1 */
                      element_type[e]=Hex8;
                      numNodesPerElement= 8;
                      element_dim=3;
                      break;
                  case 8:  /* line order 2 */
                      element_type[e]=Line3;
                      numNodesPerElement= 3;
                      element_dim=1;
                      break;
                  case 9:  /* traingle order 2 */
                      element_type[e]=Tri6;
                      numNodesPerElement= 6;
                      element_dim=2;
                      break;
                  case 10:  /* quadrilateral order 2 */
                      element_type[e]=Rec9;
                      numNodesPerElement= 9;
                      element_dim=2;
                      break;
                  case 11:  /* tetrahedron order 2 */
                      element_type[e]=Tet10;
                      numNodesPerElement= 10;
                      element_dim=3;
                      break;
                  case 15 :  /* point */
                      element_type[e]=Point1;
                      numNodesPerElement= 1;
                      element_dim=0;
                      break;
                  default:
                     element_type[e]=NoType;
                     sprintf(error_msg,"Unexected gmsh element type %d in mesh file %s.",gmsh_type,fname);
                     Finley_setError(IO_ERROR,error_msg);
              }
              if (element_dim == numDim) {
                 if (final_element_type == NoType) {
                    final_element_type = element_type[e];
                 } else if (final_element_type != element_type[e]) {
                     sprintf(error_msg,"Finley can handle a single type of internal elements only.");
                     Finley_setError(IO_ERROR,error_msg);
                     break;
                 }
                 numElements++;
              } else if (element_dim == numDim-1) {
                 if (final_face_element_type == NoType) {
                    final_face_element_type = element_type[e];
                 } else if (final_face_element_type != element_type[e]) {
                     sprintf(error_msg,"Finley can handle a single type of face elements only.");
                     Finley_setError(IO_ERROR,error_msg);
                     break;
                 }
                 numFaceElements++;
              }
              
   	   if(version <= 1.0){
   	     fscanf(fileHandle_p, "%d %d %d", &tag[e], &elementary_id, &numNodesPerElement2);
   	     partition_id = 1;
                if (numNodesPerElement2 != numNodesPerElement) {
                     sprintf(error_msg,"Illegal number of nodes for element %d in mesh file %s.",id[e],fname);
                     Finley_setError(IO_ERROR,error_msg);
                }
   	   } else {
   	     fscanf(fileHandle_p, "%d", &numTags);
   	     elementary_id = tag[e] = partition_id = 1;
                numNodesPerElement2=-1;
   	     for(j = 0; j < numTags; j++){
   	       fscanf(fileHandle_p, "%d", &itmp);	    
   	       if (j == 0) {
   	         tag[e] = itmp;
   	       } else if (j == 1) {
   	         elementary_id = itmp;
   	       } else if (j == 2) {
   	         partition_id = itmp;
                  }
   	       /* ignore any other tags */
   	     }
   	   }
              if (! Finley_noError()) break;
              for(j = 0; j < numNodesPerElement; j++) fscanf(fileHandle_p, "%d", &vertices[INDEX2(j,e,MAX_numNodes_gmsh)]);
              /* for tet10 the last two nodes need to be swapped */
              if (element_type[e]==Tet10) {
                   itmp=vertices[INDEX2(9,e,MAX_numNodes_gmsh)];
                   vertices[INDEX2(9,e,MAX_numNodes_gmsh)]=vertices[INDEX2(8,e,MAX_numNodes_gmsh)];
                   vertices[INDEX2(8,e,MAX_numNodes_gmsh)]=itmp;
              }
            }
            /* all elements have been read, now we have to identify the elements for finley */
        
            if (Finley_noError()) {
              /* first we have to identify the elements to define Elementis and FaceElements */
              if (final_element_type == NoType) {
                 if (numDim==1) {
                    final_element_type=Line2;
                 } else if (numDim==2) {
                    final_element_type=Tri3;
                 } else if (numDim==3) {
                    final_element_type=Tet4;
                 }
              }
              if (final_face_element_type == NoType) {
                 if (numDim==1) {
                    final_face_element_type=Point1;
                 } else if (numDim==2) {
                    final_face_element_type=Line2;
                 } else if (numDim==3) {
                    final_face_element_type=Tri3;
                 }
              }
              mesh_p->Elements=Finley_ElementFile_alloc(final_element_type,mesh_p->order, mesh_p->reduced_order, mpi_info);
              mesh_p->FaceElements=Finley_ElementFile_alloc(final_face_element_type,mesh_p->order, mesh_p->reduced_order, mpi_info);
              mesh_p->ContactElements=Finley_ElementFile_alloc(Point1_Contact,mesh_p->order, mesh_p->reduced_order, mpi_info);
              mesh_p->Points=Finley_ElementFile_alloc(Point1,mesh_p->order, mesh_p->reduced_order, mpi_info);
              if (Finley_noError()) {
                  Finley_ElementFile_allocTable(mesh_p->Elements, numElements);
                  Finley_ElementFile_allocTable(mesh_p->FaceElements, numFaceElements);
                  Finley_ElementFile_allocTable(mesh_p->ContactElements, 0);
                  Finley_ElementFile_allocTable(mesh_p->Points, 0);
                  if (Finley_noError()) {
                      mesh_p->Elements->minColor=0;
                      mesh_p->Elements->maxColor=numElements-1;
                      mesh_p->FaceElements->minColor=0;
                      mesh_p->FaceElements->maxColor=numFaceElements-1;
                      mesh_p->ContactElements->minColor=0;
                      mesh_p->ContactElements->maxColor=0;
                      mesh_p->Points->minColor=0;
                      mesh_p->Points->maxColor=0;
                      numElements=0;
                      numFaceElements=0;
                      for(e = 0; e < totalNumElements; e++) {
                         if (element_type[e] == final_element_type) {
                            mesh_p->Elements->Id[numElements]=id[e];
                            mesh_p->Elements->Tag[numElements]=tag[e];
                            mesh_p->Elements->Color[numElements]=numElements;
                            for (j = 0; j<  mesh_p->Elements->ReferenceElement->Type->numNodes; ++j)  {
                                  mesh_p->Elements->Nodes[INDEX2(j, numElements, mesh_p->Elements->ReferenceElement->Type->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                            }
                            numElements++;
                         } else if (element_type[e] == final_face_element_type) {
                            mesh_p->FaceElements->Id[numFaceElements]=id[e];
                            mesh_p->FaceElements->Tag[numFaceElements]=tag[e];
                            mesh_p->FaceElements->Color[numFaceElements]=numFaceElements;
                            for (j = 0; j<  mesh_p->FaceElements->ReferenceElement->Type->numNodes; ++j) {
                                     mesh_p->FaceElements->Nodes[INDEX2(j, numFaceElements, mesh_p->FaceElements->ReferenceElement->Type->numNodes)]=vertices[INDEX2(j,e,MAX_numNodes_gmsh)];
                            }
                            numFaceElements++;
                         }
                      }
                 }
              }
            }
         }
         /* and clean up */
         TMPMEMFREE(id);
         TMPMEMFREE(tag);
         TMPMEMFREE(element_type);
         TMPMEMFREE(vertices);
      }
      /* serach for end of data block */
      do {
         if (!fgets(line, sizeof(line), fileHandle_p)) {
            sprintf(error_msg,"Unexected end of file in %s",fname);
            Finley_setError(IO_ERROR,error_msg);
         }
         if (feof(fileHandle_p)) {
            sprintf(error_msg,"Unexected end of file in %s",fname);
            Finley_setError(IO_ERROR,error_msg);
         }
         if (! Finley_noError()) break;
       } while(line[0] != '$');
     }
   
     /* close file */
     fclose(fileHandle_p);
     /* clean up */
     if (! Finley_noError()) {
        Finley_Mesh_free(mesh_p);
        return NULL;
     }
     /*   resolve id's : */
     Finley_Mesh_resolveNodeIds(mesh_p);
     /* rearrange elements: */
     Finley_Mesh_prepare(mesh_p, optimize);
     /* that's it */
     #ifdef Finley_TRACE
     printf("timing: reading mesh: %.4e sec\n",Finley_timer()-time0);
     #endif
     Paso_MPIInfo_free( mpi_info );
     return mesh_p;
  }
}
