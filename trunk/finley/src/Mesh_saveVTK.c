
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

/***************************************************************************/
/*   Writes data and mesh in VTK XML format to a VTU file.                 */
/*   Nodal data needs to be given on FINLEY_NODES or FINLEY_REDUCED_NODES  */
/***************************************************************************/

#include "Mesh.h"
#include "Assemble.h"
#include "vtkCellType.h"  /* copied from vtk source directory */
#include "paso/PasoUtil.h"

#define INT_FORMAT "%d "
#define LEN_INT_FORMAT (9+1)
#define INT_NEWLINE_FORMAT "%d\n"
#define SCALAR_FORMAT "%12.6e\n"
#define VECTOR_FORMAT "%12.6e %12.6e %12.6e\n"
#define TENSOR_FORMAT "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n"
#define LEN_TENSOR_FORMAT (9*(12+1)+1)
#define NEWLINE "\n"
#define LEN_TMP_BUFFER LEN_TENSOR_FORMAT+(MAX_numNodes*LEN_INT_FORMAT+1)+2
#define NCOMP_MAX 9

#define __STRCAT(dest, chunk, dest_in_use) \
do {\
    strcpy(&dest[dest_in_use], chunk);\
    dest_in_use += strlen(chunk);\
} while(0)

#ifdef PASO_MPI

/* writes buffer to file catching the empty buffer case which causes problems
 * with some MPI versions */
#define MPI_WRITE_ORDERED(BUF, LEN) \
do {\
    if (LEN==0) { strcpy(BUF, " "); LEN=1; }\
    MPI_File_write_ordered(mpi_fileHandle_p, BUF, LEN, MPI_CHAR, &mpi_status);\
} while(0)

/* writes buffer to file on master only */
#define MPI_RANK0_WRITE_SHARED(BUF) \
do {\
    if (my_mpi_rank == 0) {\
        MPI_File_iwrite_shared(mpi_fileHandle_p, BUF, strlen(BUF), MPI_CHAR, &mpi_req);\
        MPI_Wait(&mpi_req, &mpi_status);\
    }\
} while(0)

/* For reference only. Investigation needed as to which values may improve
 * performance */
#if 0
void create_MPIInfo(MPI_Info& info)
{
    MPI_Info_create(&info);
    MPI_Info_set(info, "access_style", "write_once, sequential");
    MPI_Info_set(info, "collective_buffering", "true");
    MPI_Info_set(info, "cb_block_size",        "131072");
    MPI_Info_set(info, "cb_buffer_size",       "1048567");
    MPI_Info_set(info, "cb_nodes",             "8");
    MPI_Info_set(info, "striping_factor",      "16");
    MPI_Info_set(info, "striping_unit",        "424288");
}
#endif

#else

#define MPI_WRITE_ORDERED(A, B)
#define MPI_RANK0_WRITE_SHARED(A)

#endif /* PASO_MPI */


/* Returns one if the node given by coords and idx is within the quadrant
 * indexed by q and if the element type is Rec9 or Hex27, zero otherwise */
int nodeInQuadrant(const double *coords, ElementTypeId type, int idx, int q)
{
#define INSIDE_1D(_X_,_C_,_R_) ( ABS((_X_)-(_C_)) <= (_R_) )
#define INSIDE_2D(_X_,_Y_,_CX_,_CY_,_R_) ( INSIDE_1D(_X_,_CX_,_R_) && INSIDE_1D(_Y_,_CY_,_R_))
#define INSIDE_3D(_X_,_Y_,_Z_,_CX_,_CY_,_CZ_,_R_) ( INSIDE_1D(_X_,_CX_,_R_) && INSIDE_1D(_Y_,_CY_,_R_) && INSIDE_1D(_Z_,_CZ_,_R_) )

    int ret;
    if (type == Rec9) {
        if (q==0)
            ret = INSIDE_2D(coords[2*idx], coords[2*idx+1], 0.25, 0.25, 0.25);
        else if (q==1)
            ret = INSIDE_2D(coords[2*idx], coords[2*idx+1], 0.75, 0.25, 0.25);
        else if (q==2)
            ret = INSIDE_2D(coords[2*idx], coords[2*idx+1], 0.25, 0.75, 0.25);
        else if (q==3)
            ret = INSIDE_2D(coords[2*idx], coords[2*idx+1], 0.75, 0.75, 0.25);
        else
            ret = 0;
    } else if (type == Hex27) {
        if (q==0)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.25, 0.25, 0.25, 0.25);
        else if (q==1)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.75, 0.25, 0.25, 0.25);
        else if (q==2)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.25, 0.75, 0.25, 0.25);
        else if (q==3)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.75, 0.75, 0.25, 0.25);
        else if (q==4)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.25, 0.25, 0.75, 0.25);
        else if (q==5)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.75, 0.25, 0.75, 0.25);
        else if (q==6)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.25, 0.75, 0.75, 0.25);
        else if (q==7)
            ret = INSIDE_3D(coords[3*idx], coords[3*idx+1], coords[3*idx+2],
                    0.75, 0.75, 0.75, 0.25);
        else
            ret = 0;
    } else {
        ret = 1;
    }
    return ret;
}

void Finley_Mesh_saveVTK(const char *filename_p,
                         Finley_Mesh *mesh_p,
                         const dim_t num_data,
                         char **names_p,
                         escriptDataC **data_pp)
{
#ifdef PASO_MPI
    MPI_File mpi_fileHandle_p;
    MPI_Status mpi_status;
    MPI_Request mpi_req;
    MPI_Info mpi_info = MPI_INFO_NULL;
#endif
    Paso_MPI_rank my_mpi_rank;
    FILE *fileHandle_p = NULL;
    char errorMsg[LenErrorMsg_MAX], *txtBuffer;
    char tmpBuffer[LEN_TMP_BUFFER];
    size_t txtBufferSize, txtBufferInUse, maxNameLen;
    double *quadNodes_p = NULL;
    dim_t dataIdx, nDim;
    dim_t numCells=0, globalNumCells=0, numVTKNodesPerElement=0;
    dim_t myNumPoints=0, globalNumPoints=0;
    dim_t shape, NN=0, numCellFactor=1, myNumCells=0;
    bool_t *isCellCentered;
    bool_t writeCellData=FALSE, writePointData=FALSE, hasReducedElements=FALSE;
    index_t myFirstNode=0, myLastNode=0, *globalNodeIndex=NULL;
    index_t myFirstCell=0, k;
    int mpi_size, i, j, l;
    int cellType=0, nodeType=FINLEY_NODES, elementType=FINLEY_UNKNOWN;
    Finley_ElementFile *elements = NULL;
    ElementTypeId typeId = NoType;

    const char *vtkHeader = \
      "<?xml version=\"1.0\"?>\n" \
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n" \
      "<UnstructuredGrid>\n" \
      "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" \
      "<Points>\n" \
      "<DataArray NumberOfComponents=\"%d\" type=\"Float64\" format=\"ascii\">\n";
    char *vtkFooter = "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    const char *tag_Float_DataArray="<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n";
    char *tags_End_Points_and_Start_Conn = "</DataArray>\n</Points>\n<Cells>\n<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">\n" ;
    char *tags_End_Conn_and_Start_Offset = "</DataArray>\n<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n";
    char *tags_End_Offset_and_Start_Type = "</DataArray>\n<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n";
    char *tag_End_DataArray = "</DataArray>\n";

    const int VTK_HEX20_INDEX[] =
      { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15 };
    const int VTK_REC9_INDEX[] =
      { 0, 4, 8, 7,  4, 1, 5, 8,  7, 8, 6, 3,  8, 5, 2, 6 };
    const int VTK_HEX27_INDEX[] =
      {  0,  8, 20, 11, 12, 21, 26, 24,
         8,  1,  9, 20, 21, 13, 22, 26,
        11, 20, 10,  3, 24, 26, 23, 15,
        20,  9,  2, 10, 26, 22, 14, 23,
        12, 21, 26, 24,  4, 16, 25, 19,
        21, 13, 22, 26, 16,  5, 17, 25,
        24, 26, 23, 15, 19, 25, 18,  7,
        26, 22, 14, 23, 25, 17,  6, 18 };

    /* if there is no mesh we just return */
    if (mesh_p==NULL) return;

    nDim = mesh_p->Nodes->numDim;

    if (nDim != 2 && nDim != 3) {
        Finley_setError(TYPE_ERROR, "saveVTK: spatial dimension 2 or 3 is supported only.");
        return;
    }
    my_mpi_rank = mesh_p->Nodes->MPIInfo->rank;
    mpi_size = mesh_p->Nodes->MPIInfo->size;

    /************************************************************************/
    /* open the file and check handle */

    if (mpi_size > 1) {
#ifdef PASO_MPI
        const int amode = MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_UNIQUE_OPEN;
        int ierr;
        if (my_mpi_rank == 0 && Paso_fileExists(filename_p)) {
            remove(filename_p);
        }
        ierr = MPI_File_open(mesh_p->Nodes->MPIInfo->comm, (char*)filename_p,
                             amode, mpi_info, &mpi_fileHandle_p);
        if (ierr != MPI_SUCCESS) {
            sprintf(errorMsg, "saveVTK: File %s could not be opened for writing in parallel.", filename_p);
            Finley_setError(IO_ERROR, errorMsg);
        } else {
            MPI_File_set_view(mpi_fileHandle_p, MPI_DISPLACEMENT_CURRENT,
                    MPI_CHAR, MPI_CHAR, "native", mpi_info);
        }
#endif /* PASO_MPI */
    } else {
        fileHandle_p = fopen(filename_p, "w");
        if (fileHandle_p==NULL) {
            sprintf(errorMsg, "saveVTK: File %s could not be opened for writing.", filename_p);
            Finley_setError(IO_ERROR, errorMsg);
        }
    }
    if (!Paso_MPIInfo_noError(mesh_p->Nodes->MPIInfo)) return;

    /* General note: From this point if an error occurs Finley_setError is
     * called and subsequent steps are skipped until the end of this function
     * where allocated memory is freed and the file is closed. */

    /************************************************************************/
    /* find the mesh type to be written */

    isCellCentered = TMPMEMALLOC(num_data, bool_t);
    maxNameLen = 0;
    if (!Finley_checkPtr(isCellCentered)) {
        for (dataIdx=0; dataIdx<num_data; ++dataIdx) {
            if (! isEmpty(data_pp[dataIdx])) {
                switch(getFunctionSpaceType(data_pp[dataIdx]) ) {
                    case FINLEY_NODES:
                        nodeType = (nodeType == FINLEY_REDUCED_NODES) ? FINLEY_REDUCED_NODES : FINLEY_NODES;
                        isCellCentered[dataIdx] = FALSE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_ELEMENTS) {
                            elementType = FINLEY_ELEMENTS;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    case FINLEY_REDUCED_NODES:
                        nodeType = FINLEY_REDUCED_NODES;
                        isCellCentered[dataIdx] = FALSE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_ELEMENTS) {
                            elementType = FINLEY_ELEMENTS;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                    }
                    break;
                    case FINLEY_REDUCED_ELEMENTS:
                        hasReducedElements = TRUE;
                    case FINLEY_ELEMENTS:
                        isCellCentered[dataIdx] = TRUE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_ELEMENTS) {
                            elementType = FINLEY_ELEMENTS;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    case FINLEY_REDUCED_FACE_ELEMENTS:
                        hasReducedElements = TRUE;
                    case FINLEY_FACE_ELEMENTS:
                        isCellCentered[dataIdx] = TRUE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_FACE_ELEMENTS) {
                            elementType = FINLEY_FACE_ELEMENTS;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    case FINLEY_POINTS:
                        isCellCentered[dataIdx]=TRUE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_POINTS) {
                            elementType = FINLEY_POINTS;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    case FINLEY_REDUCED_CONTACT_ELEMENTS_1:
                        hasReducedElements = TRUE;
                    case FINLEY_CONTACT_ELEMENTS_1:
                        isCellCentered[dataIdx] = TRUE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_CONTACT_ELEMENTS_1) {
                            elementType = FINLEY_CONTACT_ELEMENTS_1;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    case FINLEY_REDUCED_CONTACT_ELEMENTS_2:
                        hasReducedElements = TRUE;
                    case FINLEY_CONTACT_ELEMENTS_2:
                        isCellCentered[dataIdx] = TRUE;
                        if (elementType==FINLEY_UNKNOWN || elementType==FINLEY_CONTACT_ELEMENTS_1) {
                            elementType = FINLEY_CONTACT_ELEMENTS_1;
                        } else {
                            Finley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
                        }
                    break;
                    default:
                        sprintf(errorMsg, "saveVTK: unknown function space type %d",getFunctionSpaceType(data_pp[dataIdx]));
                        Finley_setError(TYPE_ERROR, errorMsg);
                }
                if (isCellCentered[dataIdx]) {
                    writeCellData = TRUE;
                } else {
                    writePointData = TRUE;
                }
                maxNameLen = MAX(maxNameLen, strlen(names_p[dataIdx]));
            }
        }
    }

    /************************************************************************/
    /* select number of points and the mesh component */

    if (Finley_noError()) {
        if (nodeType == FINLEY_REDUCED_NODES) {
            myFirstNode = Finley_NodeFile_getFirstReducedNode(mesh_p->Nodes);
            myLastNode = Finley_NodeFile_getLastReducedNode(mesh_p->Nodes);
            globalNumPoints = Finley_NodeFile_getGlobalNumReducedNodes(mesh_p->Nodes);
            globalNodeIndex = Finley_NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
        } else {
            myFirstNode = Finley_NodeFile_getFirstNode(mesh_p->Nodes);
            myLastNode = Finley_NodeFile_getLastNode(mesh_p->Nodes);
            globalNumPoints = Finley_NodeFile_getGlobalNumNodes(mesh_p->Nodes);
            globalNodeIndex = Finley_NodeFile_borrowGlobalNodesIndex(mesh_p->Nodes);
        }
        myNumPoints = myLastNode - myFirstNode;
        if (elementType==FINLEY_UNKNOWN) elementType=FINLEY_ELEMENTS;
        switch(elementType) {
            case FINLEY_ELEMENTS:
                elements = mesh_p->Elements;
            break;
            case FINLEY_FACE_ELEMENTS:
                elements = mesh_p->FaceElements;
            break;
            case FINLEY_POINTS:
                elements = mesh_p->Points;
            break;
            case FINLEY_CONTACT_ELEMENTS_1:
                elements = mesh_p->ContactElements;
            break;
        }
        if (elements==NULL) {
            Finley_setError(SYSTEM_ERROR, "saveVTK: undefined element file");
        } else {
            /* map finley element type to VTK element type */
            numCells = elements->numElements;
            globalNumCells = Finley_ElementFile_getGlobalNumElements(elements);
            myNumCells = Finley_ElementFile_getMyNumElements(elements);
            myFirstCell = Finley_ElementFile_getFirstElement(elements);
            NN = elements->numNodes;
            if (nodeType==FINLEY_REDUCED_NODES) {
                typeId = elements->LinearReferenceElement->Type->TypeId;
                if (hasReducedElements) {
                    quadNodes_p=elements->LinearReferenceElementReducedOrder->QuadNodes;
                } else {
                    quadNodes_p=elements->LinearReferenceElement->QuadNodes;
                }
            } else {
                typeId = elements->ReferenceElement->Type->TypeId;
                if (hasReducedElements)
                    quadNodes_p=elements->ReferenceElementReducedOrder->QuadNodes;
                else
                    quadNodes_p=elements->ReferenceElement->QuadNodes;
            }
            switch (typeId) {
                case Point1:
                case Line2Face:
                case Line3Face:
                case Point1_Contact:
                case Line2Face_Contact:
                case Line3Face_Contact:
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
                    cellType = VTK_LINE;
                    numVTKNodesPerElement = 2;
                break;

                case Tri3:
                case Tet4Face:
                case Tet4Face_Contact:
                    cellType = VTK_TRIANGLE;
                    numVTKNodesPerElement = 3;
                break;

                case Rec4:
                case Hex8Face:
                case Rec4_Contact:
                case Hex8Face_Contact:
                    cellType = VTK_QUAD;
                    numVTKNodesPerElement = 4;
                break;

                case Rec9:
                    numCellFactor = 4;
                    cellType = VTK_QUAD;
                    numVTKNodesPerElement = 4;
                break;

                case Tet4:
                    cellType = VTK_TETRA;
                    numVTKNodesPerElement = 4;
                break;

                case Hex8:
                    cellType = VTK_HEXAHEDRON;
                    numVTKNodesPerElement = 8;
                break;

                case Line3:
                case Tri6Face:
                case Rec8Face:
                case Line3_Contact:
                case Tri6Face_Contact:
                case Rec8Face_Contact:
                    cellType = VTK_QUADRATIC_EDGE;
                    numVTKNodesPerElement = 3;
                break;

                case Tri6:
                case Tet10Face:
                case Tri6_Contact:
                case Tet10Face_Contact:
                    cellType = VTK_QUADRATIC_TRIANGLE;
                    numVTKNodesPerElement = 6;
                break;

                case Rec8:
                case Hex20Face:
                case Rec8_Contact:
                case Hex20Face_Contact:
                    cellType = VTK_QUADRATIC_QUAD;
                    numVTKNodesPerElement = 8;
                break;

                case Tet10:
                    cellType = VTK_QUADRATIC_TETRA;
                    numVTKNodesPerElement = 10;
                break;

                case Hex20:
                    cellType = VTK_QUADRATIC_HEXAHEDRON;
                    numVTKNodesPerElement = 20;
                break;

                case Hex27:
                    numCellFactor = 8;
                    cellType = VTK_HEXAHEDRON;
                    numVTKNodesPerElement = 8;
                break;

                default:
                    sprintf(errorMsg, "saveVTK: Element type %s is not supported by VTK.", elements->ReferenceElement->Type->Name);
                    Finley_setError(VALUE_ERROR, errorMsg);
            }
        }
    }

    /* allocate enough memory for text buffer */

    txtBufferSize = strlen(vtkHeader) + 3*LEN_INT_FORMAT + (30+3*maxNameLen);

    if (mpi_size > 1) {
        txtBufferSize = MAX(txtBufferSize, myNumPoints * LEN_TMP_BUFFER);
        txtBufferSize = MAX(txtBufferSize, numCellFactor * myNumCells *
                (LEN_INT_FORMAT * numVTKNodesPerElement + 1));
        txtBufferSize = MAX(txtBufferSize,
                numCellFactor * myNumCells * LEN_TENSOR_FORMAT);
        txtBufferSize = MAX(txtBufferSize, myNumPoints * LEN_TENSOR_FORMAT);
    }
    txtBuffer = TMPMEMALLOC(txtBufferSize+1, char);

    /* sets error if memory allocation failed */
    Finley_checkPtr(txtBuffer);

    /************************************************************************/
    /* write number of points and the mesh component */

    if (Finley_noError()) {
        const index_t *nodeIndex;
        if (FINLEY_REDUCED_NODES == nodeType) {
            nodeIndex = elements->ReferenceElement->Type->linearNodes;
        } else if (Rec9 == typeId) {
            nodeIndex = VTK_REC9_INDEX;
        } else if (Hex20 == typeId) {
            nodeIndex = VTK_HEX20_INDEX;
        } else if (Hex27 == typeId) {
            nodeIndex = VTK_HEX27_INDEX;
        } else if (numVTKNodesPerElement != NN) {
            nodeIndex = elements->ReferenceElement->Type->geoNodes;
        } else {
            nodeIndex = NULL;
        }

        sprintf(txtBuffer, vtkHeader, globalNumPoints,
                numCellFactor*globalNumCells, 3);

        if (mpi_size > 1) {
            /* write the nodes */
            MPI_RANK0_WRITE_SHARED(txtBuffer);
            txtBuffer[0] = '\0';
            txtBufferInUse = 0;
            if (nDim==2) {
                for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
                    if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                        sprintf(tmpBuffer, VECTOR_FORMAT,
                            mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                            0.);
                        __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                    }
                }
            } else {
                for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
                    if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                        sprintf(tmpBuffer, VECTOR_FORMAT,
                            mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
                        __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                    }
                }
            } /* nDim */
            MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);

            /* write the cells */
            MPI_RANK0_WRITE_SHARED(tags_End_Points_and_Start_Conn);
            txtBuffer[0] = '\0';
            txtBufferInUse = 0;
            if (nodeIndex == NULL) {
                for (i = 0; i < numCells; i++) {
                    if (elements->Owner[i] == my_mpi_rank) {
                        for (j = 0; j < numVTKNodesPerElement; j++) {
                            sprintf(tmpBuffer, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
                            __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                        }
                        __STRCAT(txtBuffer, NEWLINE, txtBufferInUse);
                    }
                }
            } else {
                for (i = 0; i < numCells; i++) {
                    if (elements->Owner[i] == my_mpi_rank) {
                        for (l = 0; l < numCellFactor; l++) {
                            const int* idx=&nodeIndex[l*numVTKNodesPerElement];
                            for (j = 0; j < numVTKNodesPerElement; j++) {
                                sprintf(tmpBuffer, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(idx[j], i, NN)]]);
                                __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                            }
                            __STRCAT(txtBuffer, NEWLINE, txtBufferInUse);
                        }
                    }
                }
            } /* nodeIndex */
            MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);

            /* write the offsets */
            MPI_RANK0_WRITE_SHARED(tags_End_Conn_and_Start_Offset);
            txtBuffer[0] = '\0';
            txtBufferInUse = 0;
            for (i = numVTKNodesPerElement*(myFirstCell*numCellFactor+1);
                    i <= (myFirstCell+myNumCells)*numVTKNodesPerElement*numCellFactor; i += numVTKNodesPerElement)
            {
                sprintf(tmpBuffer, INT_NEWLINE_FORMAT, i);
                __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
            }
            MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);

            /* write element type */
            sprintf(tmpBuffer, INT_NEWLINE_FORMAT, cellType);
            MPI_RANK0_WRITE_SHARED(tags_End_Offset_and_Start_Type);
            txtBuffer[0] = '\0';
            txtBufferInUse = 0;
            for (i = numVTKNodesPerElement*(myFirstCell*numCellFactor+1);
                    i <= (myFirstCell+myNumCells)*numVTKNodesPerElement*numCellFactor; i += numVTKNodesPerElement)
            {
                __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
            }
            MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);
            /* finalize cell information */
            strcpy(txtBuffer, "</DataArray>\n</Cells>\n");
            MPI_RANK0_WRITE_SHARED(txtBuffer);

        } else { /***** mpi_size == 1 *****/

            /* write the nodes */
            fputs(txtBuffer, fileHandle_p);
            if (nDim==2) {
                for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
                    if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                        fprintf(fileHandle_p, VECTOR_FORMAT,
                            mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                            0.);
                    }
                }
            } else {
                for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
                    if ( (myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode) ) {
                        fprintf(fileHandle_p, VECTOR_FORMAT,
                            mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
                            mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
                    }
                }
            } /* nDim */

            /* write the cells */
            fputs(tags_End_Points_and_Start_Conn, fileHandle_p);
            if (nodeIndex == NULL) {
                for (i = 0; i < numCells; i++) {
                    for (j = 0; j < numVTKNodesPerElement; j++) {
                        fprintf(fileHandle_p, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
                    }
                    fprintf(fileHandle_p, NEWLINE);
                }
            } else {
                for (i = 0; i < numCells; i++) {
                    for (l = 0; l < numCellFactor; l++) {
                        const int* idx=&nodeIndex[l*numVTKNodesPerElement];
                        for (j = 0; j < numVTKNodesPerElement; j++) {
                            fprintf(fileHandle_p, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(idx[j], i, NN)]]);
                        }
                        fprintf(fileHandle_p, NEWLINE);
                    }
                }
            } /* nodeIndex */

            /* write the offsets */
            fputs(tags_End_Conn_and_Start_Offset, fileHandle_p);
            for (i = numVTKNodesPerElement; i <= numCells*numVTKNodesPerElement*numCellFactor; i += numVTKNodesPerElement) {
                fprintf(fileHandle_p, INT_NEWLINE_FORMAT, i);
            }

            /* write element type */
            sprintf(tmpBuffer, INT_NEWLINE_FORMAT, cellType);
            fputs(tags_End_Offset_and_Start_Type, fileHandle_p);
            for (i = 0; i < numCells*numCellFactor; i++)
                fputs(tmpBuffer, fileHandle_p);
            /* finalize cell information */
            fputs("</DataArray>\n</Cells>\n", fileHandle_p);
        } /* mpi_size */

    } /* Finley_noError */

    /************************************************************************/
    /* write cell data */

    if (writeCellData && Finley_noError()) {
        bool_t set_scalar=FALSE, set_vector=FALSE, set_tensor=FALSE;
        /* mark the active data arrays */
        strcpy(txtBuffer, "<CellData");
        for (dataIdx = 0; dataIdx < num_data; dataIdx++) {
            if (!isEmpty(data_pp[dataIdx]) && isCellCentered[dataIdx]) {
                /* rank == 0 <--> scalar data */
                /* rank == 1 <--> vector data */
                /* rank == 2 <--> tensor data */
                switch (getDataPointRank(data_pp[dataIdx])) {
                    case 0:
                        if (!set_scalar) {
                            strcat(txtBuffer, " Scalars=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_scalar = TRUE;
                        }
                    break;
                    case 1:
                        if (!set_vector) {
                            strcat(txtBuffer, " Vectors=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_vector = TRUE;
                        }
                    break;
                    case 2:
                        if (!set_tensor) {
                            strcat(txtBuffer, " Tensors=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_tensor = TRUE;
                        }
                    break;
                    default:
                        sprintf(errorMsg, "saveVTK: data %s: VTK supports data with rank <= 2 only.", names_p[dataIdx]);
                        Finley_setError(VALUE_ERROR, errorMsg);
                }
            }
            if (!Finley_noError())
                break;
        }
    }
    /* only continue if no error occurred */
    if (writeCellData && Finley_noError()) {
        strcat(txtBuffer, ">\n");
        if ( mpi_size > 1) {
            MPI_RANK0_WRITE_SHARED(txtBuffer);
        } else {
            fputs(txtBuffer, fileHandle_p);
        }

        /* write the arrays */
        for (dataIdx = 0; dataIdx < num_data; dataIdx++) {
            if (!isEmpty(data_pp[dataIdx]) && isCellCentered[dataIdx]) {
                dim_t numPointsPerSample=getNumDataPointsPerSample(data_pp[dataIdx]);
                dim_t rank = getDataPointRank(data_pp[dataIdx]);
                dim_t nComp = getDataPointSize(data_pp[dataIdx]);
                dim_t nCompReqd = 1; /* number of components mpi_required */
                if (rank == 0) {
                    nCompReqd = 1;
                    shape = 0;
                } else if (rank == 1) {
                    nCompReqd = 3;
                    shape = getDataPointShape(data_pp[dataIdx], 0);
                    if (shape > 3) {
                        Finley_setError(VALUE_ERROR, "saveVTK: rank 1 objects must have 3 components at most.");
                    }
                } else {
                    nCompReqd = 9;
                    shape = getDataPointShape(data_pp[dataIdx], 0);
                    if (shape > 3 || shape != getDataPointShape(data_pp[dataIdx], 1)) {
                        Finley_setError(VALUE_ERROR, "saveVTK: rank 2 objects of shape 2x2 or 3x3 supported only.");
                    }
                }
                /* bail out if an error occurred */
                if (!Finley_noError())
                    break;

                sprintf(txtBuffer, tag_Float_DataArray, names_p[dataIdx], nCompReqd);
                if ( mpi_size > 1) {
                    MPI_RANK0_WRITE_SHARED(txtBuffer);
                } else {
                    fputs(txtBuffer, fileHandle_p);
                }

                txtBuffer[0] = '\0';
                txtBufferInUse = 0;
                for (i=0; i<numCells; i++) {
                    if (elements->Owner[i] == my_mpi_rank) {
                        double *values = getSampleData(data_pp[dataIdx], i);
                        for (l = 0; l < numCellFactor; l++) {
                            double sampleAvg[NCOMP_MAX];
                            dim_t nCompUsed = MIN(nComp, NCOMP_MAX);

                            /* average over number of points in the sample */
                            if (isExpanded(data_pp[dataIdx])) {
                                dim_t hits=0, hits_old;
                                for (k=0; k<nCompUsed; k++) sampleAvg[k]=0;
                                for (j=0; j<numPointsPerSample; j++) {
                                    hits_old=hits;
                                    if (nodeInQuadrant(quadNodes_p, typeId, j, l)) {
                                        hits++;
                                        for (k=0; k<nCompUsed; k++) {
                                            sampleAvg[k] += values[INDEX2(k,j,nComp)];
                                        }
                                    }
                                }
                                for (k=0; k<nCompUsed; k++)
                                    sampleAvg[k] /= MAX(hits, 1);
                            } else {
                                for (k=0; k<nCompUsed; k++)
                                    sampleAvg[k] = values[k];
                            } /* isExpanded */

                            /* if the number of required components is more than
                             * the number of actual components, pad with zeros
                             */
                            /* probably only need to get shape of first element */
                            if (nCompReqd == 1) {
                                sprintf(tmpBuffer, SCALAR_FORMAT, sampleAvg[0]);
                            } else if (nCompReqd == 3) {
                                if (shape==1) {
                                    sprintf(tmpBuffer, VECTOR_FORMAT,
                                        sampleAvg[0], 0.f, 0.f);
                                } else if (shape==2) {
                                    sprintf(tmpBuffer, VECTOR_FORMAT,
                                        sampleAvg[0], sampleAvg[1], 0.f);
                                } else if (shape==3) {
                                    sprintf(tmpBuffer, VECTOR_FORMAT,
                                        sampleAvg[0],sampleAvg[1],sampleAvg[2]);
                                }
                            } else if (nCompReqd == 9) {
                                if (shape==1) {
                                    sprintf(tmpBuffer, TENSOR_FORMAT,
                                        sampleAvg[0], 0.f, 0.f,
                                                 0.f, 0.f, 0.f,
                                                 0.f, 0.f, 0.f);
                                } else if (shape==2) {
                                    sprintf(tmpBuffer, TENSOR_FORMAT,
                                        sampleAvg[0], sampleAvg[1], 0.f,
                                        sampleAvg[2], sampleAvg[3], 0.f,
                                                 0.f,          0.f, 0.f);
                                } else if (shape==3) {
                                    sprintf(tmpBuffer, TENSOR_FORMAT,
                                        sampleAvg[0],sampleAvg[1],sampleAvg[2],
                                        sampleAvg[3],sampleAvg[4],sampleAvg[5],
                                        sampleAvg[6],sampleAvg[7],sampleAvg[8]);
                                }
                            }
                            if ( mpi_size > 1) {
                                __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                            } else {
                                fputs(tmpBuffer, fileHandle_p);
                            }
                        } /* for l (numCellFactor) */
                    } /* if I am the owner */
                } /* for i (numCells) */

                if ( mpi_size > 1) {
                    MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);
                    MPI_RANK0_WRITE_SHARED(tag_End_DataArray);
                } else {
                    fputs(tag_End_DataArray, fileHandle_p);
                }
            } /* !isEmpty && cellCentered */
        } /* for dataIdx */

        strcpy(txtBuffer, "</CellData>\n");
        if ( mpi_size > 1) {
            MPI_RANK0_WRITE_SHARED(txtBuffer);
        } else {
            fputs(txtBuffer, fileHandle_p);
        }
    } /* if noError && writeCellData */

    /************************************************************************/
    /* write point data */

    if (writePointData && Finley_noError()) {
        /* mark the active data arrays */
        bool_t set_scalar=FALSE, set_vector=FALSE, set_tensor=FALSE;
        strcpy(txtBuffer, "<PointData");
        for (dataIdx = 0; dataIdx<num_data; dataIdx++) {
            if (!isEmpty(data_pp[dataIdx]) && !isCellCentered[dataIdx]) {
                switch (getDataPointRank(data_pp[dataIdx])) {
                    case 0:
                        if (!set_scalar) {
                            strcat(txtBuffer, " Scalars=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_scalar = TRUE;
                        }
                    break;
                    case 1:
                        if (!set_vector) {
                            strcat(txtBuffer, " Vectors=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_vector = TRUE;
                        }
                    break;
                    case 2:
                        if (!set_tensor) {
                            strcat(txtBuffer, " Tensors=\"");
                            strcat(txtBuffer, names_p[dataIdx]);
                            strcat(txtBuffer, "\"");
                            set_tensor = TRUE;
                        }
                    break;
                    default:
                        sprintf(errorMsg, "saveVTK: data %s: VTK supports data with rank <= 2 only.", names_p[dataIdx]);
                        Finley_setError(VALUE_ERROR, errorMsg);
                }
            }
            if (!Finley_noError())
                break;
        }
    }
    /* only continue if no error occurred */
    if (writePointData && Finley_noError()) {
        strcat(txtBuffer, ">\n");
        if ( mpi_size > 1) {
            MPI_RANK0_WRITE_SHARED(txtBuffer);
        } else {
            fputs(txtBuffer, fileHandle_p);
        }

        /* write the arrays */
        for (dataIdx = 0; dataIdx < num_data; dataIdx++) {
            if (!isEmpty(data_pp[dataIdx]) && !isCellCentered[dataIdx]) {
                Finley_NodeMapping* nodeMapping;
                dim_t rank = getDataPointRank(data_pp[dataIdx]);
                dim_t nCompReqd = 1; /* number of components mpi_required */
                if (getFunctionSpaceType(data_pp[dataIdx]) == FINLEY_REDUCED_NODES) {
                    nodeMapping = mesh_p->Nodes->reducedNodesMapping;
                } else {
                    nodeMapping = mesh_p->Nodes->nodesMapping;
                }
                if (rank == 0) {
                    nCompReqd = 1;
                    shape = 0;
                } else if (rank == 1) {
                    nCompReqd = 3;
                    shape = getDataPointShape(data_pp[dataIdx], 0);
                    if (shape > 3) {
                        Finley_setError(VALUE_ERROR, "saveVTK: rank 1 objects must have 3 components at most.");
                    }
                } else {
                    nCompReqd = 9;
                    shape=getDataPointShape(data_pp[dataIdx], 0);
                    if (shape > 3 || shape != getDataPointShape(data_pp[dataIdx], 1)) {
                        Finley_setError(VALUE_ERROR, "saveVTK: rank 2 objects of shape 2x2 or 3x3 supported only.");
                    }
                }
                /* bail out if an error occurred */
                if (!Finley_noError())
                    break;

                sprintf(txtBuffer, tag_Float_DataArray, names_p[dataIdx], nCompReqd);
                if ( mpi_size > 1) {
                    MPI_RANK0_WRITE_SHARED(txtBuffer);
                } else {
                    fputs(txtBuffer, fileHandle_p);
                }

                txtBuffer[0] = '\0';
                txtBufferInUse = 0;
                for (i=0; i<mesh_p->Nodes->numNodes; i++) {
                    k = globalNodeIndex[i];
                    if ( (myFirstNode <= k) && (k < myLastNode) ) {
                        double *values = getSampleData(data_pp[dataIdx], nodeMapping->target[i]);
                        /* if the number of mpi_required components is more than
                         * the number of actual components, pad with zeros.
                         * Probably only need to get shape of first element */
                        if (nCompReqd == 1) {
                            sprintf(tmpBuffer, SCALAR_FORMAT, values[0]);
                        } else if (nCompReqd == 3) {
                            if (shape==1) {
                                sprintf(tmpBuffer, VECTOR_FORMAT,
                                        values[0], 0.f, 0.f);
                            } else if (shape==2) {
                                sprintf(tmpBuffer, VECTOR_FORMAT,
                                        values[0], values[1], 0.f);
                            } else if (shape==3) {
                                sprintf(tmpBuffer, VECTOR_FORMAT,
                                        values[0], values[1], values[2]);
                            }
                        } else if (nCompReqd == 9) {
                            if (shape==1) {
                                sprintf(tmpBuffer, TENSOR_FORMAT,
                                    values[0], 0.f, 0.f,
                                          0.f, 0.f, 0.f,
                                          0.f, 0.f, 0.f);
                            } else if (shape==2) {
                                sprintf(tmpBuffer, TENSOR_FORMAT,
                                    values[0], values[1], 0.f,
                                    values[2], values[3], 0.f,
                                          0.f,       0.f, 0.f);
                            } else if (shape==3) {
                                sprintf(tmpBuffer, TENSOR_FORMAT,
                                    values[0], values[1], values[2],
                                    values[3], values[4], values[5],
                                    values[6], values[7], values[8]);
                            }
                        }
                        if ( mpi_size > 1) {
                            __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
                        } else {
                            fputs(tmpBuffer, fileHandle_p);
                        }
                    } /* if this is my node */
                } /* for i (numNodes) */

                if ( mpi_size > 1) {
                    MPI_WRITE_ORDERED(txtBuffer, txtBufferInUse);
                    MPI_RANK0_WRITE_SHARED(tag_End_DataArray);
                } else {
                    fputs(tag_End_DataArray, fileHandle_p);
                }
            } /* !isEmpty && !isCellCentered */
        } /* for dataIdx */

        strcpy(txtBuffer, "</PointData>\n");
        if ( mpi_size > 1) {
            MPI_RANK0_WRITE_SHARED(txtBuffer);
        } else {
            fputs(txtBuffer, fileHandle_p);
        }
    } /* if noError && writePointData */

    /* Final write to VTK file */
    if (Finley_noError()) {
        if (mpi_size > 1) {
            MPI_RANK0_WRITE_SHARED(vtkFooter);
        } else {
            fputs(vtkFooter, fileHandle_p);
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
    TMPMEMFREE(txtBuffer);
}

