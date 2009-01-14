
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   writes data and mesh in a vtk file */
/*   nodal data needs to be given on FINLEY_NODES or FINLEY_REDUCED_NODES */

/**************************************************************/


#include "Mesh.h"
#include "Assemble.h"
#include "vtkCellType.h"  /* copied from vtk source directory !!! */
#include "paso/PasoUtil.h"

#define LEN_PRINTED_INT_FORMAT (9+1)
#define INT_FORMAT "%d "
#define INT_NEWLINE_FORMAT "%d\n"
#define FLOAT_SCALAR_FORMAT "%12.6e\n"
#define FLOAT_VECTOR_FORMAT "%12.6e %12.6e %12.6e\n"
#define FLOAT_TENSOR_FORMAT "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n"
#define LEN_PRINTED_FLOAT_SCALAR_FORMAT (12+1)
#define LEN_PRINTED_FLOAT_VECTOR_FORMAT (3*(12+1)+1)
#define LEN_PRINTED_FLOAT_TENSOR_FORMAT (9*(12+1)+1)
#define NEWLINE "\n"
#define LEN_TMP_BUFFER LEN_PRINTED_FLOAT_TENSOR_FORMAT+(MAX_numNodes*LEN_PRINTED_INT_FORMAT+1)+2
#define NCOMP_MAX 9
#define __STRCAT(dest,chunk,dest_in_use)  \
{                  \
  strcpy(&dest[dest_in_use], chunk); \
  dest_in_use+=strlen(chunk); \
}
#define INSIDE_1D(_X_,_C_,_R_) ( ABS((_X_)-(_C_)) <= (_R_) ) 
#define INSIDE_2D(_X_,_Y_,_CX_,_CY_,_R_) ( INSIDE_1D(_X_,_CX_,_R_) &&  INSIDE_1D(_Y_,_CY_,_R_))
#define INSIDE_3D(_X_,_Y_,_Z_,_CX_,_CY_,_CZ_,_R_) ( INSIDE_1D(_X_,_CX_,_R_) &&  INSIDE_1D(_Y_,_CY_,_R_) && INSIDE_1D(_Z_,_CZ_,_R_) )

void Finley_Mesh_saveVTK(const char * filename_p,
                         Finley_Mesh *mesh_p,
                         const dim_t num_data,
                         char* *names_p, 
                         escriptDataC* *data_pp)
{
#ifdef USE_VTK
  char error_msg[LenErrorMsg_MAX], *txt_buffer=NULL, tmp_buffer[LEN_TMP_BUFFER];
  __const double *values;
  double sampleAvg[NCOMP_MAX], *QuadNodes=NULL;
  size_t txt_buffer_in_use;
  dim_t len_txt_buffer,  max_len_names;
  FILE * fileHandle_p = NULL;
  int mpi_size, i, j, l, cellType=0;
  dim_t i_data, hits, hits_old;
  dim_t nDim, globalNumPoints=0, numCells=0, globalNumCells=0, numVTKNodesPerElement=0;
  dim_t myNumPoints=0, numPointsPerSample, rank, nComp, nCompReqd;
  dim_t shape, NN=0, numCellFactor=0, myNumCells=0, max_name_len;
  bool_t *isCellCentered=NULL, write_celldata=FALSE, write_pointdata=FALSE, reduced_elements=FALSE;
  bool_t set_scalar=FALSE, set_vector=FALSE, set_tensor=FALSE;
  index_t myFirstNode=0, myLastNode=0, *globalNodeIndex=NULL;
  index_t k, *node_index, myFirstCell=0;
  #ifdef PASO_MPI
  int ierr;
  /* int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY |  MPI_MODE_SEQUENTIAL;  */
  const int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN; 
  MPI_File mpi_fileHandle_p;
  MPI_Status mpi_status;
  MPI_Request mpi_req;
  MPI_Info mpi_info=MPI_INFO_NULL;
  #endif
  Paso_MPI_rank my_mpi_rank;
  int nodetype=FINLEY_NODES;
  int elementtype=FINLEY_UNKNOWN;
  Finley_NodeMapping *nodeMapping=NULL;
  Finley_ElementFile* elements=NULL;
  ElementTypeId TypeId=NoType;
  
 
  /****************************************/
  /*                                      */
  /*       tags in the vtk file           */

  const char* tags_header="<?xml version=\"1.0\"?>\n" \
                    "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n" \
                    "<UnstructuredGrid>\n" \
                    "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" \
                    "<Points>\n" \
                    "<DataArray NumberOfComponents=\"%d\" type=\"Float64\" format=\"ascii\">\n";
  char* tag_End_DataArray = "</DataArray>\n";
  char* tag_End_PointData = "</PointData>\n";
  char* tag_End_CellData =  "</CellData>\n";
  char* footer = "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  char* tags_End_Points_and_Start_Conn = "</DataArray>\n</Points>\n<Cells>\n<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">\n" ;
  char* tags_End_Conn_and_Start_Offset = "</DataArray>\n<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n";
  char* tags_End_Offset_and_Start_Type = "</DataArray>\n<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n";
  const char* tag_Float_DataArray="<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n";
  char* tags_End_Type_And_Cells = "</DataArray>\n</Cells>\n";

  int VTK_QUADRATIC_HEXAHEDRON_INDEX[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15 };

  /* if there is no mesh we just return */
  if (mesh_p==NULL) return;

  my_mpi_rank = mesh_p->Nodes->MPIInfo->rank;
  mpi_size  = mesh_p->Nodes->MPIInfo->size;
  nDim = mesh_p->Nodes->numDim;

  if (! ( (nDim ==2) || (nDim == 3) ) ) {
        Finley_setError(IO_ERROR, "saveVTK: spatial dimension 2 or 3 is supported only.");
        return;  
  }
  /*************************************************************************************/

  /* open the file and check handle */

  if (mpi_size > 1) {
        #ifdef PASO_MPI
          /* Collective Call */
          #ifdef MPIO_HINTS
            MPI_Info_create(&mpi_info);
            /*  MPI_Info_set(mpi_info, "striping_unit",        "424288"); */
            /*  MPI_Info_set(mpi_info, "striping_factor",      "16"); */
            /*  MPI_Info_set(mpi_info, "collective_buffering", "true"); */
            /*  MPI_Info_set(mpi_info, "cb_block_size",        "131072"); */
            /*  MPI_Info_set(mpi_info, "cb_buffer_size",       "1048567"); */
            /*  MPI_Info_set(mpi_info, "cb_nodes",             "8"); */
            /*    MPI_Info_set(mpi_info, "access_style", "write_once, sequential"); */
          
            /*XFS only */
            /*   MPI_Info_set(mpi_info, "direct_write",          "true"); */
          #endif
          if ( my_mpi_rank == 0) {
              if  (Paso_fileExists(filename_p)) remove(filename_p);
          }
          ierr=MPI_File_open(mesh_p->Nodes->MPIInfo->comm, (char*)filename_p, amode,mpi_info, &mpi_fileHandle_p);
          if (ierr != MPI_SUCCESS) {
	      perror(filename_p);
              sprintf(error_msg, "saveVTK: File %s could not be opened for writing in parallel.", filename_p);
              Finley_setError(IO_ERROR,error_msg);
          } else {
             MPI_File_set_view(mpi_fileHandle_p,MPI_DISPLACEMENT_CURRENT,MPI_CHAR, MPI_CHAR, "native" , mpi_info);
          }
        #endif
  } else {
        fileHandle_p = fopen(filename_p, "w");
        if (fileHandle_p==NULL) {
           sprintf(error_msg, "saveVTK: File %s could not be opened for writing.", filename_p);
           Finley_setError(IO_ERROR,error_msg);
         }
  }
  if (! Paso_MPIInfo_noError(mesh_p->Nodes->MPIInfo) ) return;
  /*************************************************************************************/

  /* find the mesh type to be written */

  isCellCentered=TMPMEMALLOC(num_data,bool_t);
  max_len_names=0;
  if (!Finley_checkPtr(isCellCentered)) {
     reduced_elements=FALSE;
     nodetype=FINLEY_UNKNOWN;
     elementtype=FINLEY_UNKNOWN;
     for (i_data=0;i_data<num_data;++i_data) {
       if (! isEmpty(data_pp[i_data])) {
         switch(getFunctionSpaceType(data_pp[i_data]) ) {
         case FINLEY_NODES:
           nodetype = (nodetype == FINLEY_REDUCED_NODES) ? FINLEY_REDUCED_NODES : FINLEY_NODES;
           isCellCentered[i_data]=FALSE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
             elementtype=FINLEY_ELEMENTS;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_REDUCED_NODES:
           nodetype = FINLEY_REDUCED_NODES;
           isCellCentered[i_data]=FALSE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
             elementtype=FINLEY_ELEMENTS;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_REDUCED_ELEMENTS:
            reduced_elements=TRUE;
         case FINLEY_ELEMENTS:
           isCellCentered[i_data]=TRUE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
             elementtype=FINLEY_ELEMENTS;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_REDUCED_FACE_ELEMENTS:
            reduced_elements=TRUE;
         case FINLEY_FACE_ELEMENTS:
           isCellCentered[i_data]=TRUE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_FACE_ELEMENTS) {
             elementtype=FINLEY_FACE_ELEMENTS;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_POINTS:
           isCellCentered[i_data]=TRUE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_POINTS) {
             elementtype=FINLEY_POINTS;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_REDUCED_CONTACT_ELEMENTS_1:
            reduced_elements=TRUE;
         case FINLEY_CONTACT_ELEMENTS_1:
           isCellCentered[i_data]=TRUE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
             elementtype=FINLEY_CONTACT_ELEMENTS_1;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         case FINLEY_REDUCED_CONTACT_ELEMENTS_2:
            reduced_elements=TRUE;
         case FINLEY_CONTACT_ELEMENTS_2:
           isCellCentered[i_data]=TRUE;
           if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
             elementtype=FINLEY_CONTACT_ELEMENTS_1;
           } else {
             Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
           }
           break;
         default:
           sprintf(error_msg,"saveVTK: unknown function space type %d",getFunctionSpaceType(data_pp[i_data]));
           Finley_setError(TYPE_ERROR,error_msg);
         }
         if (isCellCentered[i_data]) {
           write_celldata=TRUE;
         } else {
           write_pointdata=TRUE;
         }
         max_len_names =MAX(max_len_names,(dim_t)strlen(names_p[i_data]));
       }
     }
     nodetype = (nodetype == FINLEY_UNKNOWN) ? FINLEY_NODES : nodetype;
  }
  if (Finley_noError()) {

     /***************************************/

     /* select number of points and the mesh component */

     if (nodetype == FINLEY_REDUCED_NODES) {
        myFirstNode = Finley_NodeFile_getFirstReducedNode(mesh_p->Nodes);
        myLastNode = Finley_NodeFile_getLastReducedNode(mesh_p->Nodes);
        globalNumPoints = Finley_NodeFile_getGlobalNumReducedNodes(mesh_p->Nodes);
        globalNodeIndex= Finley_NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
     } else {
        myFirstNode = Finley_NodeFile_getFirstNode(mesh_p->Nodes);
        myLastNode = Finley_NodeFile_getLastNode(mesh_p->Nodes);
        globalNumPoints = Finley_NodeFile_getGlobalNumNodes(mesh_p->Nodes);
        globalNodeIndex= Finley_NodeFile_borrowGlobalNodesIndex(mesh_p->Nodes);
     }
     myNumPoints = myLastNode - myFirstNode;
     if (elementtype==FINLEY_UNKNOWN) elementtype=FINLEY_ELEMENTS;
     switch(elementtype) {
       case FINLEY_ELEMENTS:
          elements=mesh_p->Elements;
          break;
        case FINLEY_FACE_ELEMENTS:
          elements=mesh_p->FaceElements;
          break;
        case FINLEY_POINTS:
          elements=mesh_p->Points;
          break;
        case FINLEY_CONTACT_ELEMENTS_1:
          elements=mesh_p->ContactElements;
          break;
     }
     if (elements==NULL) {
       Finley_setError(SYSTEM_ERROR,"saveVTK: undefined element file");
     } else {
       /* map finley element type to VTK element type */
       numCells = elements->numElements;
       globalNumCells = Finley_ElementFile_getGlobalNumElements(elements);
       myNumCells= Finley_ElementFile_getMyNumElements(elements);
       myFirstCell= Finley_ElementFile_getFirstElement(elements);
       NN = elements->numNodes;
       if (nodetype==FINLEY_REDUCED_NODES) {
          TypeId = elements->LinearReferenceElement->Type->TypeId;
          if (reduced_elements) {
              QuadNodes=elements->LinearReferenceElementReducedOrder->QuadNodes;
          } else {
              QuadNodes=elements->LinearReferenceElement->QuadNodes;
          }
       } else {
          TypeId = elements->ReferenceElement->Type->TypeId;
          if (reduced_elements) {
              QuadNodes=elements->ReferenceElementReducedOrder->QuadNodes;
          } else {
              QuadNodes=elements->ReferenceElement->QuadNodes;
          }
       }
       switch(TypeId) {
        case Point1:
        case Line2Face:
        case Line3Face:
        case Point1_Contact:
        case Line2Face_Contact:
        case Line3Face_Contact:
          numCellFactor=1;
          cellType = VTK_VERTEX;
          numVTKNodesPerElement = 1;
          break;
      
        case Line2:
        case Tri3Face:
        case Rec4Face:
        case Line2_Contact:
        case Tri3_Contact:
        case Tri3Face_Contact:
        case Rec4Face_Contact:
          numCellFactor=1;
          cellType = VTK_LINE;
          numVTKNodesPerElement = 2;
          break;
      
        case Tri3:
        case Tet4Face:
        case Tet4Face_Contact:
          numCellFactor=1;
          cellType = VTK_TRIANGLE;
          numVTKNodesPerElement = 3;
          break;
      
        case Rec4:
        case Hex8Face:
        case Rec4_Contact:
        case Hex8Face_Contact:
          numCellFactor=1;
          cellType = VTK_QUAD;
          numVTKNodesPerElement = 4;
          break;

        case Rec9:
          numCellFactor=4;
          cellType = VTK_QUAD;
          numVTKNodesPerElement = 4;
          break;
      
        case Tet4:
          numCellFactor=1;
          cellType = VTK_TETRA;
          numVTKNodesPerElement = 4;
          break;
      
        case Hex8:
          numCellFactor=1;
          cellType = VTK_HEXAHEDRON;
          numVTKNodesPerElement = 8;
          break;
      
        case Line3:
        case Tri6Face:
        case Rec8Face:
        case Line3_Contact:
        case Tri6Face_Contact:
        case Rec8Face_Contact:
          numCellFactor=1;
          cellType = VTK_QUADRATIC_EDGE;
          numVTKNodesPerElement = 3;
          break;
      
        case Tri6:
        case Tet10Face:
        case Tri6_Contact:
        case Tet10Face_Contact:
          numCellFactor=1;
          cellType = VTK_QUADRATIC_TRIANGLE;
          numVTKNodesPerElement = 6;
          break;
      
        case Rec8:
        case Hex20Face:
        case Rec8_Contact:
        case Hex20Face_Contact:
          numCellFactor=1;
          cellType = VTK_QUADRATIC_QUAD;
          numVTKNodesPerElement = 8;
          break;
      
        case Tet10:
          numCellFactor=1;
          cellType = VTK_QUADRATIC_TETRA;
          numVTKNodesPerElement = 10;
          break;
      
        case Hex20:
          numCellFactor=1;
          cellType = VTK_QUADRATIC_HEXAHEDRON;
          numVTKNodesPerElement = 20;
          break;

        case Hex27:
          numCellFactor=8;
          cellType = VTK_HEXAHEDRON;
          numVTKNodesPerElement = 8;
          break;
      
        default:
          sprintf(error_msg, "saveVTK: Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
          Finley_setError(VALUE_ERROR,error_msg);
        }
     }
  }
  /***************************************/

  /***************************************/
  /*                                     */
  /*   allocate text buffer              */
  /*                                     */
  max_name_len=0;
  for (i_data =0 ;i_data<num_data;++i_data) max_name_len=MAX(max_name_len,(dim_t)strlen(names_p[i_data]));
  len_txt_buffer= strlen(tags_header) + 3 * LEN_PRINTED_INT_FORMAT + (30+3*max_name_len); /* header */
  if (mpi_size > 1) len_txt_buffer=MAX(len_txt_buffer, myNumPoints * LEN_TMP_BUFFER);
  if (mpi_size > 1) len_txt_buffer=MAX(len_txt_buffer, numCellFactor*myNumCells*(LEN_PRINTED_INT_FORMAT*numVTKNodesPerElement+1));
  len_txt_buffer=MAX(len_txt_buffer,200+3*max_len_names);
  len_txt_buffer=MAX(len_txt_buffer, (dim_t)strlen(tag_Float_DataArray) + LEN_PRINTED_INT_FORMAT + max_len_names);
  if (mpi_size > 1) len_txt_buffer=MAX(len_txt_buffer, numCellFactor*myNumCells*LEN_PRINTED_FLOAT_TENSOR_FORMAT);
  if (mpi_size > 1) len_txt_buffer=MAX(len_txt_buffer, myNumPoints*LEN_PRINTED_FLOAT_TENSOR_FORMAT);
  txt_buffer=TMPMEMALLOC(len_txt_buffer+1,char);
  Finley_checkPtr(txt_buffer);
  
  if (Finley_noError()) {

     /* select number of points and the mesh component */

     sprintf(txt_buffer,tags_header,globalNumPoints,numCellFactor*globalNumCells,3);

      if (mpi_size > 1) {
          if ( my_mpi_rank == 0) {
            #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,txt_buffer,strlen(txt_buffer),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
            #endif
          }
      } else {
         fprintf(fileHandle_p, "%s", txt_buffer);
      }

      /* write the nodes */
      
      if (mpi_size > 1) {

         txt_buffer[0] = '\0';
         txt_buffer_in_use=0;
         if (nDim==2) {
            for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
               if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                 sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,
                                    mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                                    mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                                    0.);
                 __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
               }
            }      
         } else {
            for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
               if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                 sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,
                                                  mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                                                  mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                                                  mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
                 __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
               }
            }    
   
         }
         #ifdef PASO_MPI
            if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
            MPI_File_write_ordered(mpi_fileHandle_p, txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
         #endif     
      } else {
         if (nDim==2) {
            for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
               if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                 fprintf(fileHandle_p,FLOAT_VECTOR_FORMAT,
                                      mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                                      mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                                      0.);
               }
            }      
         } else {
            for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
               if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                 fprintf(fileHandle_p,FLOAT_VECTOR_FORMAT,
                                              mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                                              mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                                              mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
               }
            }    
   
         }
      }

      /* close the Points and open connectivity */

      if (mpi_size > 1) {
          if ( my_mpi_rank == 0) {
             #ifdef PASO_MPI
                MPI_File_iwrite_shared(mpi_fileHandle_p, tags_End_Points_and_Start_Conn, strlen(tags_End_Points_and_Start_Conn), MPI_CHAR, &mpi_req);
                MPI_Wait(&mpi_req,&mpi_status);
             #endif
          }
      } else {
         fprintf(fileHandle_p, "%s", tags_End_Points_and_Start_Conn);
      }

     /* write the cells */
     if (nodetype == FINLEY_REDUCED_NODES) {
        node_index=elements->ReferenceElement->Type->linearNodes;
     } else if (VTK_QUADRATIC_HEXAHEDRON==cellType) {
        node_index=VTK_QUADRATIC_HEXAHEDRON_INDEX;
     } else if ( (numVTKNodesPerElement!=NN) && (TypeId!=Rec9) && (TypeId!=Hex27) ) {
        node_index=elements->ReferenceElement->Type->geoNodes;
     } else {
        node_index=NULL;
     }

     if ( mpi_size > 1) {
        txt_buffer[0] = '\0';
        txt_buffer_in_use=0;
        if (node_index == NULL) {
           if (TypeId==Rec9) {
              for (i = 0; i < numCells; i++) {
                 if (elements->Owner[i] == my_mpi_rank) {
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(0, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(4, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(7, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(4, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(1, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(5, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(7, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(6, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(3, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(5, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(2, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(6, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)
                }
              }
           } else if (TypeId==Hex27) {
              for (i = 0; i < numCells; i++) {
                 if (elements->Owner[i] == my_mpi_rank) {
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 0, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 8, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(11, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(12, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 8, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 1, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 9, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(13, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(11, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(10, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 3, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(15, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 9, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 2, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(10, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(14, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(12, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 4, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(16, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(19, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(13, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(16, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 5, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(17, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(15, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(19, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(18, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 7, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)

                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(14, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(17, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 6, i, NN)]]);
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(18, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                        __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)
                 }
              }
           } else {
              for (i = 0; i < numCells; i++) {
                 if (elements->Owner[i] == my_mpi_rank) {
                    for (j = 0; j < numVTKNodesPerElement; j++) {
                        sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
                        __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                    }
                    __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)
                 } 
              }
           }
        } else {
           for (i = 0; i < numCells; i++) {
              if (elements->Owner[i] == my_mpi_rank) {
                 for (j = 0; j < numVTKNodesPerElement; j++) {
                     sprintf(tmp_buffer,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(node_index[j], i, NN)]]);
                     __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use)
                 }
                 __STRCAT(txt_buffer,NEWLINE,txt_buffer_in_use)
              }
           }
        }
        #ifdef PASO_MPI
           if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
           MPI_File_write_ordered(mpi_fileHandle_p,txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
        #endif     
     } else {
        if (node_index == NULL) {
           if (TypeId==Rec9) {
              for (i = 0; i < numCells; i++) {
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(0, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(4, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(7, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(4, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(1, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(5, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(7, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(6, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(3, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(8, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(5, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(2, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(6, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);
              }

           } else if (TypeId==Hex27) {
                 for (i = 0; i < numCells; i++) {
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 0, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 8, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(11, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(12, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 8, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 1, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 9, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(13, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(11, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(10, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 3, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(15, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(20, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 9, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 2, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(10, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(14, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(12, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 4, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(16, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(19, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(21, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(13, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(16, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 5, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(17, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(24, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(15, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(19, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(18, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 7, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);

                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(26, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(22, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(14, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(23, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(25, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(17, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2( 6, i, NN)]]);
                        fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(18, i, NN)]]);
                        fprintf(fileHandle_p,NEWLINE);
                 }
              } else {
           for (i = 0; i < numCells; i++) {
                 for (j = 0; j < numVTKNodesPerElement; j++) {
                    fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
                  }
                 fprintf(fileHandle_p,NEWLINE);
              }
          }
        } else {
           for (i = 0; i < numCells; i++) {
              for (j = 0; j < numVTKNodesPerElement; j++) {
                 fprintf(fileHandle_p,INT_FORMAT,globalNodeIndex[elements->Nodes[INDEX2(node_index[j], i, NN)]]);
               }
              fprintf(fileHandle_p,NEWLINE);
           }
        }

     }
     /* finalize the connection and start the offset section */
     if (mpi_size > 1) {
        if( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,tags_End_Conn_and_Start_Offset,strlen(tags_End_Conn_and_Start_Offset),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
     } else {
        fprintf(fileHandle_p, "%s", tags_End_Conn_and_Start_Offset);
     }

    /* write the offsets */
      
     if ( mpi_size > 1) {
        txt_buffer[0] = '\0';
        txt_buffer_in_use=0;
        for (i=numVTKNodesPerElement*(myFirstCell*numCellFactor+1); i<=(myFirstCell+myNumCells)*numVTKNodesPerElement*numCellFactor; i+=numVTKNodesPerElement) {
            sprintf(tmp_buffer, INT_NEWLINE_FORMAT, i);
            __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
         }
         #ifdef PASO_MPI
            if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
            MPI_File_write_ordered(mpi_fileHandle_p,txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
         #endif     
     } else {
        for (i=numVTKNodesPerElement; i<=numCells*numVTKNodesPerElement*numCellFactor; i+=numVTKNodesPerElement) {
           fprintf(fileHandle_p, INT_NEWLINE_FORMAT, i);
        }
     
     }
     /* finalize the offset section and start the type section */
     if ( mpi_size > 1) {
        if ( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,tags_End_Offset_and_Start_Type,strlen(tags_End_Offset_and_Start_Type),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
    } else {
       fprintf(fileHandle_p, "%s", tags_End_Offset_and_Start_Type);
    }
     /* write element type */
     sprintf(tmp_buffer, INT_NEWLINE_FORMAT, cellType);
     if ( mpi_size > 1) {
        txt_buffer[0] = '\0';
        txt_buffer_in_use=0;
        for (i=numVTKNodesPerElement*(myFirstCell*numCellFactor+1); i<=(myFirstCell+myNumCells)*numVTKNodesPerElement*numCellFactor; i+=numVTKNodesPerElement) {
          __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
	}
         #ifdef PASO_MPI
            if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
            MPI_File_write_ordered(mpi_fileHandle_p,txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
         #endif     
     } else {
        for (i=0; i<numCells*numCellFactor; i++) fprintf(fileHandle_p, "%s",  tmp_buffer);
     }
     /* finalize cell information */
     if ( mpi_size > 1) {
        if ( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,tags_End_Type_And_Cells,strlen(tags_End_Type_And_Cells),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
    } else {
       fprintf(fileHandle_p, "%s", tags_End_Type_And_Cells);
    }
 }

 /* Write cell data */
 if (write_celldata && Finley_noError()) {
      /* mark the active data arrays */
      txt_buffer[0] = '\0';
      set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
      strcat(txt_buffer, "<CellData");
      for (i_data =0 ;i_data<num_data;++i_data) {
        if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data]) {
          /* if the rank == 0:   --> scalar data */
          /* if the rank == 1:   --> vector data */
          /* if the rank == 2:   --> tensor data */

          switch(getDataPointRank(data_pp[i_data])) {
          case 0:
            if (! set_scalar) {
              strcat(txt_buffer," Scalars=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_scalar=TRUE;
            }
            break;
          case 1:
            if (! set_vector) {
              strcat(txt_buffer," Vectors=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_vector=TRUE;
            }
            break;
          case 2:
            if (! set_tensor) {
              strcat(txt_buffer," Tensors=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_tensor=TRUE;
            }
            break;
          default:
            sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
            Finley_setError(VALUE_ERROR,error_msg);
            return;
          }
        }
      }
      strcat(txt_buffer, ">\n");
      if ( mpi_size > 1) {
        if ( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,txt_buffer,strlen(txt_buffer),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
      } else {
          fprintf(fileHandle_p, "%s", txt_buffer);
      }
      /* write the arrays */
      for (i_data =0 ;i_data<num_data;++i_data) {
         if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data]) {
	    void* buffer=allocSampleBuffer(data_pp[i_data]);
            txt_buffer[0] = '\0';
            txt_buffer_in_use=0;
            numPointsPerSample=getNumDataPointsPerSample(data_pp[i_data]);
            rank = getDataPointRank(data_pp[i_data]);
            nComp = getDataPointSize(data_pp[i_data]);
            nCompReqd=1;   /* the number of components mpi_required by vtk */
            shape=0;
            if (rank == 0) {
              nCompReqd = 1;
            } else if (rank == 1) {
              shape=getDataPointShape(data_pp[i_data], 0);
              if  (shape>3) {
                Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
              }
              nCompReqd = 3;
            } else {
              shape=getDataPointShape(data_pp[i_data], 0);
              if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1)) {
                Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
              }
              nCompReqd = 9;
            }
            if (Finley_noError()) {
               sprintf(txt_buffer,tag_Float_DataArray,names_p[i_data], nCompReqd);
               if ( mpi_size > 1) {
                 if ( my_mpi_rank == 0) {
                    #ifdef PASO_MPI
                       MPI_File_iwrite_shared(mpi_fileHandle_p,txt_buffer,strlen(txt_buffer),MPI_CHAR,&mpi_req);
                       MPI_Wait(&mpi_req,&mpi_status);
                    #endif
                 }
               } else {
                   fprintf(fileHandle_p, "%s", txt_buffer);
               }

               for (i=0; i<numCells; i++) {
                   if (elements->Owner[i] == my_mpi_rank) {
                      values = getSampleDataRO(data_pp[i_data], i,buffer);
                      for (l=0; l< numCellFactor;++l) {
                         /* averaging over the number of points in the sample */
                         if (isExpanded(data_pp[i_data])) {
                              for (k=0; k<MIN(nComp,NCOMP_MAX); k++) sampleAvg[k]=0;
                              hits=0;
                              for (j=0; j<numPointsPerSample; j++) {
                                 hits_old=hits;
                                 if (TypeId==Rec9) {
                                    switch(l) {
                                      case 0:
                                        if (INSIDE_2D(QuadNodes[2*j],QuadNodes[2*j+1],0.25,0.25,0.25)) hits++;  
                                        break;
                                      case 1:
                                        if (INSIDE_2D(QuadNodes[2*j],QuadNodes[2*j+1],0.75,0.25,0.25)) hits++;  
                                        break;
                                      case 2:
                                        if (INSIDE_2D(QuadNodes[2*j],QuadNodes[2*j+1],0.25,0.75,0.25)) hits++;  
                                        break;
                                      case 3:
                                        if (INSIDE_2D(QuadNodes[2*j],QuadNodes[2*j+1],0.75,0.75,0.25)) hits++;  
                                        break;
                                      }
                                 } else if (TypeId==Hex27) {
                                    switch(l) {
                                      case 0:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.25,0.25,0.25,0.25)) hits++;  
                                        break;
                                      case 1:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.75,0.25,0.25,0.25)) hits++;  
                                        break;
                                      case 2:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.25,0.75,0.25,0.25)) hits++;  
                                        break;
                                      case 3:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.75,0.75,0.25,0.25)) hits++;  
                                        break;
                                      case 4:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.25,0.25,0.75,0.25)) hits++;  
                                        break;
                                      case 5:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.75,0.25,0.75,0.25)) hits++;  
                                        break;
                                      case 6:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.25,0.75,0.75,0.25)) hits++;  
                                        break;
                                      case 7:
                                        if (INSIDE_3D(QuadNodes[3*j],QuadNodes[3*j+1],QuadNodes[3*j+2],0.75,0.75,0.75,0.25)) hits++;  
                                        break;
                                    }
                                 } else {
                                    hits++;
                                 }
                                 if (hits_old<hits) for (k=0; k<MIN(nComp,NCOMP_MAX); k++) {
                                     sampleAvg[k] += values[INDEX2(k,j,nComp)];
                                 }
                              }
                              for (k=0; k<MIN(nComp,NCOMP_MAX); k++) sampleAvg[k] /=MAX(hits,1);
                         } else {
                              for (k=0; k<MIN(nComp,NCOMP_MAX); k++) sampleAvg[k] = values[k];
                         }
                         /* if the number of required components is more than the number
                         * of actual components, pad with zeros
                         */
                         /* probably only need to get shape of first element */
                         /* write the data different ways for scalar, vector and tensor */
                         if (nCompReqd == 1) {
                           sprintf(tmp_buffer,FLOAT_SCALAR_FORMAT,sampleAvg[0]);
                         } else if (nCompReqd == 3) { 
                           if (shape==1) {
                            sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,sampleAvg[0],0.,0.);
                           } else if (shape==2) {
                            sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,sampleAvg[0],sampleAvg[1],0.);
                           } else if (shape==3) {
                            sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,sampleAvg[0],sampleAvg[1],sampleAvg[2]);
                           }
                         } else if (nCompReqd == 9) {
                           if (shape==1) {
                            sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,sampleAvg[0],0.,0.,
                                                                   0.,0.,0.,
                                                                0.,0.,0.);
                           } else if (shape==2) {
                            sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,sampleAvg[0],sampleAvg[1],0.,
                                                                   sampleAvg[2],sampleAvg[3],0.,
                                                                   0.,0.,0.);
                           } else if (shape==3) {
                            sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,sampleAvg[0],sampleAvg[1],sampleAvg[2],
                                                                   sampleAvg[3],sampleAvg[4],sampleAvg[5],
                                                                   sampleAvg[6],sampleAvg[7],sampleAvg[8]);
                           }
                         }
                         /* this needs a bit mor work!!! */
                            if ( mpi_size > 1) {
                              __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
                            } else {
                              fprintf(fileHandle_p, "%s", tmp_buffer);
                            }
                        }
                   }
               }
               if ( mpi_size > 1) {
                     #ifdef PASO_MPI
                        if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
                        MPI_File_write_ordered(mpi_fileHandle_p,txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
                     #endif     
                     if ( my_mpi_rank == 0) {
                        #ifdef PASO_MPI
                           MPI_File_iwrite_shared(mpi_fileHandle_p,tag_End_DataArray,strlen(tag_End_DataArray),MPI_CHAR,&mpi_req);
                           MPI_Wait(&mpi_req,&mpi_status);
                        #endif
                     }
               } else {
                   fprintf(fileHandle_p, "%s", tag_End_DataArray);
               }
            }
	    freeSampleBuffer(buffer);
         }
      }
      if ( mpi_size > 1) {
        if ( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,tag_End_CellData,strlen(tag_End_CellData),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
      } else {
          fprintf(fileHandle_p, "%s", tag_End_CellData);
      }
  }
  /* point data */
  if (write_pointdata && Finley_noError()) {
      /* mark the active data arrays */
      set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
      txt_buffer[0] = '\0';
      strcat(txt_buffer, "<PointData");
      for (i_data =0 ;i_data<num_data;++i_data) {
        if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data]) {
          /* if the rank == 0:   --> scalar data */
          /* if the rank == 1:   --> vector data */
          /* if the rank == 2:   --> tensor data */

          switch(getDataPointRank(data_pp[i_data])) {
          case 0:
            if (! set_scalar) {
              strcat(txt_buffer," Scalars=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_scalar=TRUE;
            }
            break;
          case 1:
            if (! set_vector) {
              strcat(txt_buffer," Vectors=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_vector=TRUE;
            }
            break;
          case 2:
            if (! set_tensor) {
              strcat(txt_buffer," Tensors=\"");
              strcat(txt_buffer,names_p[i_data]);
              strcat(txt_buffer,"\"");
              set_tensor=TRUE;
            }
            break;
          default:
            sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
            Finley_setError(VALUE_ERROR,error_msg);
            return;
          }
        }
      }
      strcat(txt_buffer, ">\n");
      if ( mpi_size > 1) {
        if ( my_mpi_rank == 0) {
           #ifdef PASO_MPI
              MPI_File_iwrite_shared(mpi_fileHandle_p,txt_buffer,strlen(txt_buffer),MPI_CHAR,&mpi_req);
              MPI_Wait(&mpi_req,&mpi_status);
           #endif
        }
      } else {
          fprintf(fileHandle_p, "%s", txt_buffer);
      }
      /* write the arrays */
      for (i_data =0 ;i_data<num_data;++i_data) {
         if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data]) {
	    void* buffer=allocSampleBuffer(data_pp[i_data]);	
            txt_buffer[0] = '\0';
            txt_buffer_in_use=0;
            numPointsPerSample=getNumDataPointsPerSample(data_pp[i_data]);
            rank = getDataPointRank(data_pp[i_data]);
            nComp = getDataPointSize(data_pp[i_data]);
            if (getFunctionSpaceType(data_pp[i_data]) == FINLEY_REDUCED_NODES) {
               nodeMapping=mesh_p->Nodes->reducedNodesMapping;
            } else {
               nodeMapping=mesh_p->Nodes->nodesMapping;
            }
            nCompReqd=1;   /* the number of components mpi_required by vtk */
            shape=0;
            if (rank == 0) {
              nCompReqd = 1;
            } else if (rank == 1) {
              shape=getDataPointShape(data_pp[i_data], 0);
              if  (shape>3) {
                Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
              }
              nCompReqd = 3;
            } else {
              shape=getDataPointShape(data_pp[i_data], 0);
              if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1)) {
                Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
              }
              nCompReqd = 9;
            }
            if (Finley_noError()) {
               sprintf(txt_buffer,tag_Float_DataArray,names_p[i_data], nCompReqd);
               if ( mpi_size > 1) {
                 if ( my_mpi_rank == 0) {
                    #ifdef PASO_MPI
                       MPI_File_iwrite_shared(mpi_fileHandle_p,txt_buffer,strlen(txt_buffer),MPI_CHAR,&mpi_req);
                       MPI_Wait(&mpi_req,&mpi_status);
                    #endif
                 }
               } else {
                   fprintf(fileHandle_p, "%s", txt_buffer);
               }
               for (i=0; i<mesh_p->Nodes->numNodes; i++) {
                  k=globalNodeIndex[i];
                  if ( (myFirstNode <= k) && (k < myLastNode) ) {
                     values = getSampleDataRO(data_pp[i_data], nodeMapping->target[i],buffer);
                     /* if the number of mpi_required components is more than the number
                     * of actual components, pad with zeros
                     */
                     /* probably only need to get shape of first element */
                     /* write the data different ways for scalar, vector and tensor */
                     if (nCompReqd == 1) {
                       sprintf(tmp_buffer,FLOAT_SCALAR_FORMAT,values[0]);
                     } else if (nCompReqd == 3) { 
                       if (shape==1) {
                        sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,values[0],0.,0.);
                       } else if (shape==2) {
                        sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,values[0],values[1],0.);
                       } else if (shape==3) {
                        sprintf(tmp_buffer,FLOAT_VECTOR_FORMAT,values[0],values[1],values[2]);
                       }
                     } else if (nCompReqd == 9) {
                       if (shape==1) {
                        sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,values[0],0.,0.,
                                                               0.,0.,0.,
                                                               0.,0.,0.);
                       } else if (shape==2) {
                        sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,values[0],values[1],0.,
                                                               values[2],values[3],0.,
                                                               0.,0.,0.);
                       } else if (shape==3) {
                        sprintf(tmp_buffer,FLOAT_TENSOR_FORMAT,values[0],values[1],values[2],
                                                               values[3],values[4],values[5],
                                                               values[6],values[7],values[8]);
                       }
                     }
                     if ( mpi_size > 1) {
                       __STRCAT(txt_buffer,tmp_buffer,txt_buffer_in_use);
                     } else {
                       fprintf(fileHandle_p, "%s", tmp_buffer);
                     }
                  }
               }
               if ( mpi_size > 1) {
                   #ifdef PASO_MPI
                     if (txt_buffer_in_use==0) { strcpy(txt_buffer, " "); txt_buffer_in_use = 1; } /* avoid zero-length writes */
                     MPI_File_write_ordered(mpi_fileHandle_p,txt_buffer,txt_buffer_in_use, MPI_CHAR, &mpi_status);
                   #endif     
                   if ( my_mpi_rank == 0) {
                      #ifdef PASO_MPI
                         MPI_File_iwrite_shared(mpi_fileHandle_p,tag_End_DataArray,strlen(tag_End_DataArray),MPI_CHAR,&mpi_req);
                         MPI_Wait(&mpi_req,&mpi_status);
                      #endif
                   }
               } else {
                  fprintf(fileHandle_p, "%s", tag_End_DataArray);
               }
            }
	    freeSampleBuffer(buffer);
          }
        }
        if ( mpi_size > 1) {
          if ( my_mpi_rank == 0) {
             #ifdef PASO_MPI
                MPI_File_iwrite_shared(mpi_fileHandle_p,tag_End_PointData,strlen(tag_End_PointData),MPI_CHAR,&mpi_req);
                MPI_Wait(&mpi_req,&mpi_status);
             #endif
          }
        } else {
            fprintf(fileHandle_p, "%s", tag_End_PointData);
        }
  }
  if (Finley_noError()) {
     if ( mpi_size > 1) {
       if ( my_mpi_rank == 0) {
          #ifdef PASO_MPI
             MPI_File_iwrite_shared(mpi_fileHandle_p,footer,strlen(footer),MPI_CHAR,&mpi_req);
             MPI_Wait(&mpi_req,&mpi_status);
             #ifdef MPIO_HINTS
               MPI_Info_free(&mpi_info);
               #undef MPIO_HINTS
             #endif
          #endif
        }
     } else {
         fprintf(fileHandle_p, "%s", footer);
     }
  }
  if ( mpi_size > 1) {
#ifdef PASO_MPI
     MPI_File_close(&mpi_fileHandle_p);
#endif
  } else {
     fclose(fileHandle_p);
  }
  TMPMEMFREE(isCellCentered);
  TMPMEMFREE(txt_buffer);
#else
  /* Don't kill the job if saveVTK() doesn't work */
  fprintf(stderr, "\n\nsaveVTK warning: VTK is not available\n\n\n");
#endif
}
