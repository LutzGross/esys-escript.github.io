
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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
/*   Nodal data needs to be given on DUDLEY_NODES or DUDLEY_REDUCED_NODES  */
/***************************************************************************/

#include "Mesh.h"
#include "Assemble.h"
#include "vtkCellType.h"	/* copied from vtk source directory */
#include "paso/PasoUtil.h"

#define INT_FORMAT "%d "
#define LEN_INT_FORMAT (unsigned int)(9+2)
#define INT_NEWLINE_FORMAT "%d\n"
#define SCALAR_FORMAT "%12.6e\n"
#define VECTOR_FORMAT "%12.6e %12.6e %12.6e\n"
#define TENSOR_FORMAT "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n"
/* strlen("-1.234567e+789 ") == 15 */
#define LEN_TENSOR_FORMAT (unsigned int)(9*15+2)
#define NEWLINE "\n"
/* This value is pulled from finley's ReferenceElements.h */
#define MAX_numNodes 64
#define LEN_TMP_BUFFER LEN_TENSOR_FORMAT+(MAX_numNodes*LEN_INT_FORMAT+1)+2
#define NCOMP_MAX (unsigned int)9

#define __STRCAT(dest, chunk, dest_in_use) \
do {\
    strcpy(&dest[dest_in_use], chunk);\
    dest_in_use += strlen(chunk);\
} while(0)

#ifdef ESYS_MPI
/* writes buffer to file catching the empty buffer case which causes problems
 * with some MPI versions */
#define MPI_WRITE_ORDERED(BUF) \
do {\
    int LLEN=0; \
    LLEN=(int) strlen(BUF); \
    if (LLEN==0) { strcpy(BUF, ""); LLEN=0; }\
    MPI_File_write_ordered(mpi_fileHandle_p, BUF, LLEN, MPI_CHAR, &mpi_status);\
} while(0)

/* writes buffer to file on master only */
#define MPI_RANK0_WRITE_SHARED(BUF) \
do {\
    int LLEN=0; \
    if (my_mpi_rank == 0) {\
        LLEN=(int) strlen(BUF); \
        if (LLEN==0) { strcpy(BUF,""); LLEN=0; }\
        MPI_File_iwrite_shared(mpi_fileHandle_p, BUF, LLEN, MPI_CHAR, &mpi_req);\
        MPI_Wait(&mpi_req, &mpi_status);\
    }\
} while(0)

/* For reference only. Investigation needed as to which values may improve
 * performance */
#if 0
void create_MPIInfo(MPI_Info & info)
{
    MPI_Info_create(&info);
    MPI_Info_set(info, "access_style", "write_once, sequential");
    MPI_Info_set(info, "collective_buffering", "true");
    MPI_Info_set(info, "cb_block_size", "131072");
    MPI_Info_set(info, "cb_buffer_size", "1048567");
    MPI_Info_set(info, "cb_nodes", "8");
    MPI_Info_set(info, "striping_factor", "16");
    MPI_Info_set(info, "striping_unit", "424288");
}
#endif

#else

#define MPI_WRITE_ORDERED(A)
#define MPI_RANK0_WRITE_SHARED(A)

#endif				/* ESYS_MPI */

#include "ShapeTable.h"

void Dudley_Mesh_saveVTK(const char *filename_p,
			 Dudley_Mesh * mesh_p,
			 const dim_t num_data,
			 char **names_p, escriptDataC ** data_pp, const char *metadata, const char *metadata_schema)
{
#ifdef ESYS_MPI
    MPI_File mpi_fileHandle_p;
    MPI_Status mpi_status;
    MPI_Request mpi_req;
    MPI_Info mpi_info = MPI_INFO_NULL;
#endif
    Esys_MPI_rank my_mpi_rank;
    FILE *fileHandle_p = NULL;
    char errorMsg[LenErrorMsg_MAX], *txtBuffer;
    char tmpBuffer[LEN_TMP_BUFFER];
    size_t txtBufferSize, txtBufferInUse, maxNameLen;
    const double *quadNodes_p = NULL;
    dim_t dataIdx, nDim;
    dim_t numCells = 0, globalNumCells = 0, numVTKNodesPerElement = 0;
    dim_t myNumPoints = 0, globalNumPoints = 0;
    dim_t shape, NN = 0, numCellFactor = 1, myNumCells = 0;
    bool_t *isCellCentered;
    bool_t writeCellData = FALSE, writePointData = FALSE, hasReducedElements = FALSE;
    index_t myFirstNode = 0, myLastNode = 0, *globalNodeIndex = NULL;
    index_t myFirstCell = 0, k;
    int mpi_size, i, j, l;
    int cellType = 0, nodeType = DUDLEY_NODES, elementType = DUDLEY_UNKNOWN;
    Dudley_ElementFile *elements = NULL;
    Dudley_ElementTypeId typeId = Dudley_NoRef;

    const char *vtkHeader =
	"<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"%s%s>\n%s%s"
	"<UnstructuredGrid>\n"
	"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n"
	"<Points>\n" "<DataArray NumberOfComponents=\"%d\" type=\"Float64\" format=\"ascii\">\n";
    char *vtkFooter = "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    const char *tag_Float_DataArray =
	"<DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n";
    char *tags_End_Points_and_Start_Conn =
	"</DataArray>\n</Points>\n<Cells>\n<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">\n";
    char *tags_End_Conn_and_Start_Offset =
	"</DataArray>\n<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n";
    char *tags_End_Offset_and_Start_Type = "</DataArray>\n<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n";
    char *tag_End_DataArray = "</DataArray>\n";

    /* if there is no mesh we just return */
    if (mesh_p == NULL)
	return;

    nDim = mesh_p->Nodes->numDim;

    if (nDim != 2 && nDim != 3)
    {
	Dudley_setError(TYPE_ERROR, "saveVTK: spatial dimension 2 or 3 is supported only.");
	return;
    }
    my_mpi_rank = mesh_p->Nodes->MPIInfo->rank;
    mpi_size = mesh_p->Nodes->MPIInfo->size;

    /************************************************************************
     * open the file and check handle *
     */
    if (mpi_size > 1)
    {
#ifdef ESYS_MPI
	const int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN;
	int ierr;
	if (my_mpi_rank == 0 && Paso_fileExists(filename_p))
	{
	    remove(filename_p);
	}
	ierr = MPI_File_open(mesh_p->Nodes->MPIInfo->comm, (char *)filename_p, amode, mpi_info, &mpi_fileHandle_p);
	if (ierr != MPI_SUCCESS)
	{
	    sprintf(errorMsg, "saveVTK: File %s could not be opened for writing in parallel.", filename_p);
	    Dudley_setError(IO_ERROR, errorMsg);
	}
	else
	{
	    ierr = MPI_File_set_view(mpi_fileHandle_p, MPI_DISPLACEMENT_CURRENT,
				     MPI_CHAR, MPI_CHAR, "native", mpi_info);
	}
#endif				/* ESYS_MPI */
    }
    else
    {
	fileHandle_p = fopen(filename_p, "w");
	if (fileHandle_p == NULL)
	{
	    sprintf(errorMsg, "saveVTK: File %s could not be opened for writing.", filename_p);
	    Dudley_setError(IO_ERROR, errorMsg);
	}
    }
    if (!Esys_MPIInfo_noError(mesh_p->Nodes->MPIInfo))
	return;

    /* General note: From this point if an error occurs Dudley_setError is
     * called and subsequent steps are skipped until the end of this function
     * where allocated memory is freed and the file is closed. */

    /************************************************************************/
    /* find the mesh type to be written */

    isCellCentered = TMPMEMALLOC(num_data, bool_t);
    maxNameLen = 0;
    if (!Dudley_checkPtr(isCellCentered))
    {
	for (dataIdx = 0; dataIdx < num_data; ++dataIdx)
	{
	    if (!isEmpty(data_pp[dataIdx]))
	    {
		switch (getFunctionSpaceType(data_pp[dataIdx]))
		{
		case DUDLEY_NODES:
		    nodeType = (nodeType == DUDLEY_REDUCED_NODES) ? DUDLEY_REDUCED_NODES : DUDLEY_NODES;
		    isCellCentered[dataIdx] = FALSE;
		    if (elementType == DUDLEY_UNKNOWN || elementType == DUDLEY_ELEMENTS)
		    {
			elementType = DUDLEY_ELEMENTS;
		    }
		    else
		    {
			Dudley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
		    }
		    break;
		case DUDLEY_REDUCED_NODES:
		    nodeType = DUDLEY_REDUCED_NODES;
		    isCellCentered[dataIdx] = FALSE;
		    if (elementType == DUDLEY_UNKNOWN || elementType == DUDLEY_ELEMENTS)
		    {
			elementType = DUDLEY_ELEMENTS;
		    }
		    else
		    {
			Dudley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
		    }
		    break;
		case DUDLEY_REDUCED_ELEMENTS:
		    hasReducedElements = TRUE;
		case DUDLEY_ELEMENTS:
		    isCellCentered[dataIdx] = TRUE;
		    if (elementType == DUDLEY_UNKNOWN || elementType == DUDLEY_ELEMENTS)
		    {
			elementType = DUDLEY_ELEMENTS;
		    }
		    else
		    {
			Dudley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
		    }
		    break;
		case DUDLEY_REDUCED_FACE_ELEMENTS:
		    hasReducedElements = TRUE;
		case DUDLEY_FACE_ELEMENTS:
		    isCellCentered[dataIdx] = TRUE;
		    if (elementType == DUDLEY_UNKNOWN || elementType == DUDLEY_FACE_ELEMENTS)
		    {
			elementType = DUDLEY_FACE_ELEMENTS;
		    }
		    else
		    {
			Dudley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
		    }
		    break;
		case DUDLEY_POINTS:
		    isCellCentered[dataIdx] = TRUE;
		    if (elementType == DUDLEY_UNKNOWN || elementType == DUDLEY_POINTS)
		    {
			elementType = DUDLEY_POINTS;
		    }
		    else
		    {
			Dudley_setError(TYPE_ERROR, "saveVTK: cannot write given data in single file.");
		    }
		    break;
		default:
		    sprintf(errorMsg, "saveVTK: unknown function space type %d",
			    getFunctionSpaceType(data_pp[dataIdx]));
		    Dudley_setError(TYPE_ERROR, errorMsg);
		}
		if (isCellCentered[dataIdx])
		{
		    writeCellData = TRUE;
		}
		else
		{
		    writePointData = TRUE;
		}
		maxNameLen = MAX(maxNameLen, strlen(names_p[dataIdx]));
	    }
	}
    }

    /************************************************************************/
    /* select number of points and the mesh component */

    if (Dudley_noError())
    {
	if (nodeType == DUDLEY_REDUCED_NODES)
	{
	    myFirstNode = Dudley_NodeFile_getFirstReducedNode(mesh_p->Nodes);
	    myLastNode = Dudley_NodeFile_getLastReducedNode(mesh_p->Nodes);
	    globalNumPoints = Dudley_NodeFile_getGlobalNumReducedNodes(mesh_p->Nodes);
	    globalNodeIndex = Dudley_NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
	}
	else
	{
	    myFirstNode = Dudley_NodeFile_getFirstNode(mesh_p->Nodes);
	    myLastNode = Dudley_NodeFile_getLastNode(mesh_p->Nodes);
	    globalNumPoints = Dudley_NodeFile_getGlobalNumNodes(mesh_p->Nodes);
	    globalNodeIndex = Dudley_NodeFile_borrowGlobalNodesIndex(mesh_p->Nodes);
	}
	myNumPoints = myLastNode - myFirstNode;
	if (elementType == DUDLEY_UNKNOWN)
	    elementType = DUDLEY_ELEMENTS;
	switch (elementType)
	{
	case DUDLEY_ELEMENTS:
	    elements = mesh_p->Elements;
	    break;
	case DUDLEY_FACE_ELEMENTS:
	    elements = mesh_p->FaceElements;
	    break;
	case DUDLEY_POINTS:
	    elements = mesh_p->Points;
	    break;
	}
	if (elements == NULL)
	{
	    Dudley_setError(SYSTEM_ERROR, "saveVTK: undefined element file");
	}
	else
	{
	    /* map dudley element type to VTK element type */
	    numCells = elements->numElements;
	    globalNumCells = Dudley_ElementFile_getGlobalNumElements(elements);
	    myNumCells = Dudley_ElementFile_getMyNumElements(elements);
	    myFirstCell = Dudley_ElementFile_getFirstElement(elements);
	    NN = elements->numNodes;
	    if (!getQuadShape(elements->numLocalDim, hasReducedElements, &quadNodes_p))
	    {
		Dudley_setError(TYPE_ERROR, "Unable to locate shape function.");
	    }

/*            if (hasReducedElements) {

                quadNodes_p=elements->referenceElementSet->referenceElementReducedQuadrature->BasisFunctions->QuadNodes;
            } else {
                quadNodes_p=elements->referenceElementSet->referenceElement->BasisFunctions->QuadNodes;
            }*/
	    if (nodeType == DUDLEY_REDUCED_NODES)
	    {
		typeId = elements->etype;	/*referenceElementSet->referenceElement->Type->TypeId;*/
	    }
	    else
	    {
		typeId = elements->etype;	/*referenceElementSet->referenceElement->Type->TypeId;*/
	    }
	    switch (typeId)
	    {
	    case Dudley_Point1:
	    case Dudley_Line2Face:
		cellType = VTK_VERTEX;
		numVTKNodesPerElement = 1;
		break;

	    case Dudley_Line2:
	    case Dudley_Tri3Face:
		cellType = VTK_LINE;
		numVTKNodesPerElement = 2;
		break;

	    case Dudley_Tri3:
	    case Dudley_Tet4Face:
		cellType = VTK_TRIANGLE;
		numVTKNodesPerElement = 3;
		break;

	    case Dudley_Tet4:
		cellType = VTK_TETRA;
		numVTKNodesPerElement = 4;
		break;

	    default:
		sprintf(errorMsg, "saveVTK: Element type %s is not supported by VTK.",
			elements->ename /*referenceElementSet->referenceElement->Type->Name */ );
		Dudley_setError(VALUE_ERROR, errorMsg);
	    }
	}
    }

    /* allocate enough memory for text buffer */

    txtBufferSize =
	strlen(vtkHeader) + 3 * LEN_INT_FORMAT + (30 + 3 * maxNameLen) + strlen(metadata) + strlen(metadata_schema);
    if (mpi_size > 1)
    {
	txtBufferSize = MAX(txtBufferSize, myNumPoints * LEN_TMP_BUFFER);
	txtBufferSize = MAX(txtBufferSize, numCellFactor * myNumCells * (LEN_INT_FORMAT * numVTKNodesPerElement + 1));
	txtBufferSize = MAX(txtBufferSize, numCellFactor * myNumCells * LEN_TENSOR_FORMAT);
	txtBufferSize = MAX(txtBufferSize, myNumPoints * LEN_TENSOR_FORMAT);
    }
    txtBuffer = TMPMEMALLOC(txtBufferSize + 1, char);

    /* sets error if memory allocation failed */
    Dudley_checkPtr(txtBuffer);

    /************************************************************************/
    /* write number of points and the mesh component */

    if (Dudley_noError())
    {
	int flag1 = 0;
	if (DUDLEY_REDUCED_NODES == nodeType)
	{
	    flag1 = 1;
	}
	else if (numVTKNodesPerElement !=
		 elements->numNodes /*referenceElementSet->referenceElement->Type->numNodes */ )
	{
	    flag1 = 1;
	}
	if (strlen(metadata) > 0)
	{
	    if (strlen(metadata_schema) > 0)
	    {
		sprintf(txtBuffer, vtkHeader, " ", metadata_schema, metadata, "\n", globalNumPoints,
			numCellFactor * globalNumCells, 3);
	    }
	    else
	    {
		sprintf(txtBuffer, vtkHeader, "", "", metadata, "\n", globalNumPoints, numCellFactor * globalNumCells,
			3);
	    }
	}
	else
	{
	    if (strlen(metadata_schema) > 0)
	    {
		sprintf(txtBuffer, vtkHeader, " ", metadata_schema, "", "", globalNumPoints,
			numCellFactor * globalNumCells, 3);
	    }
	    else
	    {
		sprintf(txtBuffer, vtkHeader, "", "", "", "", globalNumPoints, numCellFactor * globalNumCells, 3);
	    }
	}

	if (mpi_size > 1)
	{
	    /* write the nodes */
	    MPI_RANK0_WRITE_SHARED(txtBuffer);
	    txtBuffer[0] = '\0';
	    txtBufferInUse = 0;
	    if (nDim == 2)
	    {
		for (i = 0; i < mesh_p->Nodes->numNodes; i++)
		{
		    if ((myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode))
		    {
			sprintf(tmpBuffer, VECTOR_FORMAT,
				mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)], 0.);
			__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
		    }
		}
	    }
	    else
	    {
		for (i = 0; i < mesh_p->Nodes->numNodes; i++)
		{
		    if ((myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode))
		    {
			sprintf(tmpBuffer, VECTOR_FORMAT,
				mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
			__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
		    }
		}
	    }			/* nDim */
	    MPI_WRITE_ORDERED(txtBuffer);

	    /* write the cells */
	    MPI_RANK0_WRITE_SHARED(tags_End_Points_and_Start_Conn);
	    txtBuffer[0] = '\0';
	    txtBufferInUse = 0;
	    if (!flag1)
	    {
		for (i = 0; i < numCells; i++)
		{
		    if (elements->Owner[i] == my_mpi_rank)
		    {
			for (j = 0; j < numVTKNodesPerElement; j++)
			{
			    sprintf(tmpBuffer, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
			    __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
			}
			__STRCAT(txtBuffer, NEWLINE, txtBufferInUse);
		    }
		}
	    }
	    else
	    {
		for (i = 0; i < numCells; i++)
		{
		    if (elements->Owner[i] == my_mpi_rank)
		    {
			for (l = 0; l < numCellFactor; l++)
			{
			    const int idx = l * numVTKNodesPerElement;
			    for (j = 0; j < numVTKNodesPerElement; j++)
			    {
				sprintf(tmpBuffer, INT_FORMAT,
					globalNodeIndex[elements->Nodes[INDEX2(idx + j, i, NN)]]);
				__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
			    }
			    __STRCAT(txtBuffer, NEWLINE, txtBufferInUse);
			}
		    }
		}
	    }			/* nodeIndex */
	    MPI_WRITE_ORDERED(txtBuffer);

	    /* write the offsets */
	    MPI_RANK0_WRITE_SHARED(tags_End_Conn_and_Start_Offset);
	    txtBuffer[0] = '\0';
	    txtBufferInUse = 0;
	    for (i = numVTKNodesPerElement * (myFirstCell * numCellFactor + 1);
		 i <= (myFirstCell + myNumCells) * numVTKNodesPerElement * numCellFactor; i += numVTKNodesPerElement)
	    {
		sprintf(tmpBuffer, INT_NEWLINE_FORMAT, i);
		__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
	    }
	    MPI_WRITE_ORDERED(txtBuffer);

	    /* write element type */
	    sprintf(tmpBuffer, INT_NEWLINE_FORMAT, cellType);
	    MPI_RANK0_WRITE_SHARED(tags_End_Offset_and_Start_Type);
	    txtBuffer[0] = '\0';
	    txtBufferInUse = 0;
	    for (i = numVTKNodesPerElement * (myFirstCell * numCellFactor + 1);
		 i <= (myFirstCell + myNumCells) * numVTKNodesPerElement * numCellFactor; i += numVTKNodesPerElement)
	    {
		__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
	    }
	    MPI_WRITE_ORDERED(txtBuffer);
	    /* finalize cell information */
	    strcpy(txtBuffer, "</DataArray>\n</Cells>\n");
	    MPI_RANK0_WRITE_SHARED(txtBuffer);

	}
	else
	{	 /***** mpi_size == 1 *****/

	    /* write the nodes */
	    fputs(txtBuffer, fileHandle_p);
	    if (nDim == 2)
	    {
		for (i = 0; i < mesh_p->Nodes->numNodes; i++)
		{
		    if ((myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode))
		    {
			fprintf(fileHandle_p, VECTOR_FORMAT,
				mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)], 0.);
		    }
		}
	    }
	    else
	    {
		for (i = 0; i < mesh_p->Nodes->numNodes; i++)
		{
		    if ((myFirstNode <= globalNodeIndex[i]) && (globalNodeIndex[i] < myLastNode))
		    {
			fprintf(fileHandle_p, VECTOR_FORMAT,
				mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(1, i, nDim)],
				mesh_p->Nodes->Coordinates[INDEX2(2, i, nDim)]);
		    }
		}
	    }			/* nDim */

	    /* write the cells */
	    fputs(tags_End_Points_and_Start_Conn, fileHandle_p);
	    if (!flag1)
	    {
		for (i = 0; i < numCells; i++)
		{
		    for (j = 0; j < numVTKNodesPerElement; j++)
		    {
			fprintf(fileHandle_p, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(j, i, NN)]]);
		    }
		    fprintf(fileHandle_p, NEWLINE);
		}
	    }
	    else
	    {
		for (i = 0; i < numCells; i++)
		{
		    for (l = 0; l < numCellFactor; l++)
		    {
			const int idx = l * numVTKNodesPerElement;
			for (j = 0; j < numVTKNodesPerElement; j++)
			{
			    fprintf(fileHandle_p, INT_FORMAT, globalNodeIndex[elements->Nodes[INDEX2(idx + j, i, NN)]]);
			}
			fprintf(fileHandle_p, NEWLINE);
		    }
		}
	    }			/* nodeIndex */

	    /* write the offsets */
	    fputs(tags_End_Conn_and_Start_Offset, fileHandle_p);
	    for (i = numVTKNodesPerElement; i <= numCells * numVTKNodesPerElement * numCellFactor;
		 i += numVTKNodesPerElement)
	    {
		fprintf(fileHandle_p, INT_NEWLINE_FORMAT, i);
	    }

	    /* write element type */
	    sprintf(tmpBuffer, INT_NEWLINE_FORMAT, cellType);
	    fputs(tags_End_Offset_and_Start_Type, fileHandle_p);
	    for (i = 0; i < numCells * numCellFactor; i++)
		fputs(tmpBuffer, fileHandle_p);
	    /* finalize cell information */
	    fputs("</DataArray>\n</Cells>\n", fileHandle_p);
	}			/* mpi_size */

    }

    /* Dudley_noError */
 /************************************************************************/
    /* write cell data */
    if (writeCellData && Dudley_noError())
    {
	bool_t set_scalar = FALSE, set_vector = FALSE, set_tensor = FALSE;
	/* mark the active data arrays */
	strcpy(txtBuffer, "<CellData");
	for (dataIdx = 0; dataIdx < num_data; dataIdx++)
	{
	    if (!isEmpty(data_pp[dataIdx]) && isCellCentered[dataIdx])
	    {
		/* rank == 0 <--> scalar data */
		/* rank == 1 <--> vector data */
		/* rank == 2 <--> tensor data */
		switch (getDataPointRank(data_pp[dataIdx]))
		{
		case 0:
		    if (!set_scalar)
		    {
			strcat(txtBuffer, " Scalars=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_scalar = TRUE;
		    }
		    break;
		case 1:
		    if (!set_vector)
		    {
			strcat(txtBuffer, " Vectors=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_vector = TRUE;
		    }
		    break;
		case 2:
		    if (!set_tensor)
		    {
			strcat(txtBuffer, " Tensors=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_tensor = TRUE;
		    }
		    break;
		default:
		    sprintf(errorMsg, "saveVTK: data %s: VTK supports data with rank <= 2 only.", names_p[dataIdx]);
		    Dudley_setError(VALUE_ERROR, errorMsg);
		}
	    }
	    if (!Dudley_noError())
		break;
	}
    }
    /* only continue if no error occurred */
    if (writeCellData && Dudley_noError())
    {
	strcat(txtBuffer, ">\n");
	if (mpi_size > 1)
	{
	    MPI_RANK0_WRITE_SHARED(txtBuffer);
	}
	else
	{
	    fputs(txtBuffer, fileHandle_p);
	}

	/* write the arrays */
	for (dataIdx = 0; dataIdx < num_data; dataIdx++)
	{
	    if (!isEmpty(data_pp[dataIdx]) && isCellCentered[dataIdx])
	    {
		dim_t numPointsPerSample = getNumDataPointsPerSample(data_pp[dataIdx]);
		dim_t rank = getDataPointRank(data_pp[dataIdx]);
		dim_t nComp = getDataPointSize(data_pp[dataIdx]);
		dim_t nCompReqd = 1;	/* number of components mpi_required */
		if (rank == 0)
		{
		    nCompReqd = 1;
		    shape = 0;
		}
		else if (rank == 1)
		{
		    nCompReqd = 3;
		    shape = getDataPointShape(data_pp[dataIdx], 0);
		    if (shape > 3)
		    {
			Dudley_setError(VALUE_ERROR, "saveVTK: rank 1 objects must have 3 components at most.");
		    }
		}
		else
		{
		    nCompReqd = 9;
		    shape = getDataPointShape(data_pp[dataIdx], 0);
		    if (shape > 3 || shape != getDataPointShape(data_pp[dataIdx], 1))
		    {
			Dudley_setError(VALUE_ERROR, "saveVTK: rank 2 objects of shape 2x2 or 3x3 supported only.");
		    }
		}
		/* bail out if an error occurred */
		if (!Dudley_noError())
		    break;

		sprintf(txtBuffer, tag_Float_DataArray, names_p[dataIdx], nCompReqd);
		if (mpi_size > 1)
		{
		    MPI_RANK0_WRITE_SHARED(txtBuffer);
		}
		else
		{
		    fputs(txtBuffer, fileHandle_p);
		}

		txtBuffer[0] = '\0';
		txtBufferInUse = 0;
		for (i = 0; i < numCells; i++)
		{
		    if (elements->Owner[i] == my_mpi_rank)
		    {
			__const double *values = getSampleDataRO(data_pp[dataIdx], i);
			for (l = 0; l < numCellFactor; l++)
			{
			    double sampleAvg[NCOMP_MAX];
			    dim_t nCompUsed = MIN(nComp, NCOMP_MAX);

			    /* average over number of points in the sample */
			    if (isExpanded(data_pp[dataIdx]))
			    {
				dim_t hits = 0;
				for (k = 0; k < nCompUsed; k++)
				    sampleAvg[k] = 0;
				for (j = 0; j < numPointsPerSample; j++)
				{
					hits++;
					for (k = 0; k < nCompUsed; k++)
					{
					    sampleAvg[k] += values[INDEX2(k, j, nComp)];
					}
				}
				for (k = 0; k < nCompUsed; k++)
				    sampleAvg[k] /= MAX(hits, 1);
			    }
			    else
			    {
				for (k = 0; k < nCompUsed; k++)
				    sampleAvg[k] = values[k];
			    }	/* isExpanded */

			    /* if the number of required components is more than
			     * the number of actual components, pad with zeros
			     */
			    /* probably only need to get shape of first element */
			    if (nCompReqd == 1)
			    {
				sprintf(tmpBuffer, SCALAR_FORMAT, sampleAvg[0]);
			    }
			    else if (nCompReqd == 3)
			    {
				if (shape == 1)
				{
				    sprintf(tmpBuffer, VECTOR_FORMAT, sampleAvg[0], 0.f, 0.f);
				}
				else if (shape == 2)
				{
				    sprintf(tmpBuffer, VECTOR_FORMAT, sampleAvg[0], sampleAvg[1], 0.f);
				}
				else if (shape == 3)
				{
				    sprintf(tmpBuffer, VECTOR_FORMAT, sampleAvg[0], sampleAvg[1], sampleAvg[2]);
				}
			    }
			    else if (nCompReqd == 9)
			    {
				if (shape == 1)
				{
				    sprintf(tmpBuffer, TENSOR_FORMAT,
					    sampleAvg[0], 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
				}
				else if (shape == 2)
				{
				    sprintf(tmpBuffer, TENSOR_FORMAT,
					    sampleAvg[0], sampleAvg[1], 0.f,
					    sampleAvg[2], sampleAvg[3], 0.f, 0.f, 0.f, 0.f);
				}
				else if (shape == 3)
				{
				    sprintf(tmpBuffer, TENSOR_FORMAT,
					    sampleAvg[0], sampleAvg[1], sampleAvg[2],
					    sampleAvg[3], sampleAvg[4], sampleAvg[5],
					    sampleAvg[6], sampleAvg[7], sampleAvg[8]);
				}
			    }
			    if (mpi_size > 1)
			    {
				__STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
			    }
			    else
			    {
				fputs(tmpBuffer, fileHandle_p);
			    }
			}	/* for l (numCellFactor) */
		    }		/* if I am the owner */
		}		/* for i (numCells) */

		if (mpi_size > 1)
		{
		    MPI_WRITE_ORDERED(txtBuffer);
		    MPI_RANK0_WRITE_SHARED(tag_End_DataArray);
		}
		else
		{
		    fputs(tag_End_DataArray, fileHandle_p);
		}
	    }			/* !isEmpty && cellCentered */
	}			/* for dataIdx */

	strcpy(txtBuffer, "</CellData>\n");
	if (mpi_size > 1)
	{
	    MPI_RANK0_WRITE_SHARED(txtBuffer);
	}
	else
	{
	    fputs(txtBuffer, fileHandle_p);
	}
    }

    /* if noError && writeCellData */
 /************************************************************************/
    /* write point data */
    if (writePointData && Dudley_noError())
    {
	/* mark the active data arrays */
	bool_t set_scalar = FALSE, set_vector = FALSE, set_tensor = FALSE;
	strcpy(txtBuffer, "<PointData");
	for (dataIdx = 0; dataIdx < num_data; dataIdx++)
	{
	    if (!isEmpty(data_pp[dataIdx]) && !isCellCentered[dataIdx])
	    {
		switch (getDataPointRank(data_pp[dataIdx]))
		{
		case 0:
		    if (!set_scalar)
		    {
			strcat(txtBuffer, " Scalars=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_scalar = TRUE;
		    }
		    break;
		case 1:
		    if (!set_vector)
		    {
			strcat(txtBuffer, " Vectors=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_vector = TRUE;
		    }
		    break;
		case 2:
		    if (!set_tensor)
		    {
			strcat(txtBuffer, " Tensors=\"");
			strcat(txtBuffer, names_p[dataIdx]);
			strcat(txtBuffer, "\"");
			set_tensor = TRUE;
		    }
		    break;
		default:
		    sprintf(errorMsg, "saveVTK: data %s: VTK supports data with rank <= 2 only.", names_p[dataIdx]);
		    Dudley_setError(VALUE_ERROR, errorMsg);
		}
	    }
	    if (!Dudley_noError())
		break;
	}
    }
    /* only continue if no error occurred */
    if (writePointData && Dudley_noError())
    {
	strcat(txtBuffer, ">\n");
	if (mpi_size > 1)
	{
	    MPI_RANK0_WRITE_SHARED(txtBuffer);
	}
	else
	{
	    fputs(txtBuffer, fileHandle_p);
	}

	/* write the arrays */
	for (dataIdx = 0; dataIdx < num_data; dataIdx++)
	{
	    if (!isEmpty(data_pp[dataIdx]) && !isCellCentered[dataIdx])
	    {
		Dudley_NodeMapping *nodeMapping;
		dim_t rank = getDataPointRank(data_pp[dataIdx]);
		dim_t nCompReqd = 1;	/* number of components mpi_required */
		if (getFunctionSpaceType(data_pp[dataIdx]) == DUDLEY_REDUCED_NODES)
		{
		    nodeMapping = mesh_p->Nodes->reducedNodesMapping;
		}
		else
		{
		    nodeMapping = mesh_p->Nodes->nodesMapping;
		}
		if (rank == 0)
		{
		    nCompReqd = 1;
		    shape = 0;
		}
		else if (rank == 1)
		{
		    nCompReqd = 3;
		    shape = getDataPointShape(data_pp[dataIdx], 0);
		    if (shape > 3)
		    {
			Dudley_setError(VALUE_ERROR, "saveVTK: rank 1 objects must have 3 components at most.");
		    }
		}
		else
		{
		    nCompReqd = 9;
		    shape = getDataPointShape(data_pp[dataIdx], 0);
		    if (shape > 3 || shape != getDataPointShape(data_pp[dataIdx], 1))
		    {
			Dudley_setError(VALUE_ERROR, "saveVTK: rank 2 objects of shape 2x2 or 3x3 supported only.");
		    }
		}
		/* bail out if an error occurred */
		if (!Dudley_noError())
		    break;

		sprintf(txtBuffer, tag_Float_DataArray, names_p[dataIdx], nCompReqd);
		if (mpi_size > 1)
		{
		    MPI_RANK0_WRITE_SHARED(txtBuffer);
		}
		else
		{
		    fputs(txtBuffer, fileHandle_p);
		}

		txtBuffer[0] = '\0';
		txtBufferInUse = 0;
		for (i = 0; i < mesh_p->Nodes->numNodes; i++)
		{
		    k = globalNodeIndex[i];
		    if ((myFirstNode <= k) && (k < myLastNode))
		    {
			__const double *values = getSampleDataRO(data_pp[dataIdx], nodeMapping->target[i]);
			/* if the number of mpi_required components is more than
			 * the number of actual components, pad with zeros.
			 * Probably only need to get shape of first element */
			if (nCompReqd == 1)
			{
			    sprintf(tmpBuffer, SCALAR_FORMAT, values[0]);
			}
			else if (nCompReqd == 3)
			{
			    if (shape == 1)
			    {
				sprintf(tmpBuffer, VECTOR_FORMAT, values[0], 0.f, 0.f);
			    }
			    else if (shape == 2)
			    {
				sprintf(tmpBuffer, VECTOR_FORMAT, values[0], values[1], 0.f);
			    }
			    else if (shape == 3)
			    {
				sprintf(tmpBuffer, VECTOR_FORMAT, values[0], values[1], values[2]);
			    }
			}
			else if (nCompReqd == 9)
			{
			    if (shape == 1)
			    {
				sprintf(tmpBuffer, TENSOR_FORMAT, values[0], 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
			    }
			    else if (shape == 2)
			    {
				sprintf(tmpBuffer, TENSOR_FORMAT,
					values[0], values[1], 0.f, values[2], values[3], 0.f, 0.f, 0.f, 0.f);
			    }
			    else if (shape == 3)
			    {
				sprintf(tmpBuffer, TENSOR_FORMAT,
					values[0], values[1], values[2],
					values[3], values[4], values[5], values[6], values[7], values[8]);
			    }
			}
			if (mpi_size > 1)
			{
			    __STRCAT(txtBuffer, tmpBuffer, txtBufferInUse);
			}
			else
			{
			    fputs(tmpBuffer, fileHandle_p);
			}
		    }		/* if this is my node */
		}		/* for i (numNodes) */

		if (mpi_size > 1)
		{
		    MPI_WRITE_ORDERED(txtBuffer);
		    MPI_RANK0_WRITE_SHARED(tag_End_DataArray);
		}
		else
		{
		    fputs(tag_End_DataArray, fileHandle_p);
		}
	    }			/* !isEmpty && !isCellCentered */
	}			/* for dataIdx */

	strcpy(txtBuffer, "</PointData>\n");
	if (mpi_size > 1)
	{
	    MPI_RANK0_WRITE_SHARED(txtBuffer);
	}
	else
	{
	    fputs(txtBuffer, fileHandle_p);
	}
    }

    /* if noError && writePointData */
    /* Final write to VTK file */
    if (Dudley_noError())
    {
	if (mpi_size > 1)
	{
	    MPI_RANK0_WRITE_SHARED(vtkFooter);
	}
	else
	{
	    fputs(vtkFooter, fileHandle_p);
	}
    }

    if (mpi_size > 1)
    {
#ifdef ESYS_MPI
	MPI_File_close(&mpi_fileHandle_p);
	MPI_Barrier(mesh_p->Nodes->MPIInfo->comm);
#endif
    }
    else
    {
	fclose(fileHandle_p);
    }
    TMPMEMFREE(isCellCentered);
    TMPMEMFREE(txtBuffer);
}
