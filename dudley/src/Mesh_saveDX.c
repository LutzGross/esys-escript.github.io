
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/**************************************************************/

/*   writes data and mesh in an opendx file                      */
/*   the input data needs to be cell centered or on reducedNodes */

/**************************************************************/

#include "Mesh.h"
#include "Assemble.h"

/**************************************************************/

void Dudley_Mesh_saveDX(const char *filename_p, Dudley_Mesh * mesh_p, const dim_t num_data, char **names_p,
			escriptDataC * *data_pp)
{
    char error_msg[LenErrorMsg_MAX], elemTypeStr[32];
    /* some tables needed for reordering */
    int resort[6][9] = {
	{0, 1},			/* line */
	{0, 1, 2},		/* triagle */
	{0, 1, 2, 3},		/* tetrahedron */
	{0, 3, 1, 2},		/* quadrilateral */
	{3, 0, 7, 4, 2, 1, 6, 5},	/* hexahedron */
    };
    FILE *fileHandle_p = NULL;
    int i, j, k, i_data, elementtype, numPoints = 0, nDim, *resortIndex = NULL, p,
	numDXNodesPerElement = 0, numCells, NN, object_count, rank, nComp, numPointsPerSample;
    __const double *values;
    double rtmp;
    bool_t *isCellCentered = NULL;
    Dudley_ElementFile *elements = NULL;
    Dudley_ElementTypeId TypeId;
    /* open the file  and check handel */

    /* if there is no mesh we just return */
    if (mesh_p == NULL)
	return;
    isCellCentered = MEMALLOC(num_data, bool_t);
    if (Dudley_checkPtr(isCellCentered))
	return;

    fileHandle_p = fopen(filename_p, "w");
    if (fileHandle_p == NULL)
    {
	sprintf(error_msg, "File %s could not be opened for writing.", filename_p);
	MEMFREE(isCellCentered);
	fclose(fileHandle_p);
	Dudley_setError(IO_ERROR, error_msg);
	return;
    }
    /* find the mesh type to be written */
    elementtype = DUDLEY_UNKNOWN;
    for (i_data = 0; i_data < num_data; ++i_data)
    {
	if (!isEmpty(data_pp[i_data]))
	{
	    switch (getFunctionSpaceType(data_pp[i_data]))
	    {
	    case DUDLEY_REDUCED_NODES:
		if (elementtype == DUDLEY_UNKNOWN || elementtype == DUDLEY_ELEMENTS)
		{
		    elementtype = DUDLEY_ELEMENTS;
		}
		else
		{
		    Dudley_setError(TYPE_ERROR, "saveDX: cannot write given data into single file.");
		    MEMFREE(isCellCentered);
		    fclose(fileHandle_p);
		    return;
		}
		isCellCentered[i_data] = FALSE;
		break;
	    case DUDLEY_ELEMENTS:
	    case DUDLEY_REDUCED_ELEMENTS:
		if (elementtype == DUDLEY_UNKNOWN || elementtype == DUDLEY_ELEMENTS)
		{
		    elementtype = DUDLEY_ELEMENTS;
		}
		else
		{
		    Dudley_setError(TYPE_ERROR, "saveDX: cannot write given data into single file.");
		    MEMFREE(isCellCentered);
		    fclose(fileHandle_p);
		    return;
		}
		isCellCentered[i_data] = TRUE;
		break;
	    case DUDLEY_FACE_ELEMENTS:
	    case DUDLEY_REDUCED_FACE_ELEMENTS:
		if (elementtype == DUDLEY_UNKNOWN || elementtype == DUDLEY_FACE_ELEMENTS)
		{
		    elementtype = DUDLEY_FACE_ELEMENTS;
		}
		else
		{
		    Dudley_setError(TYPE_ERROR, "saveDX: cannot write given data into single file.");
		    MEMFREE(isCellCentered);
		    fclose(fileHandle_p);
		    return;
		}
		isCellCentered[i_data] = TRUE;
		break;
	    case DUDLEY_POINTS:
		if (elementtype == DUDLEY_UNKNOWN || elementtype == DUDLEY_POINTS)
		{
		    elementtype = DUDLEY_POINTS;
		}
		else
		{
		    Dudley_setError(TYPE_ERROR, "saveDX: cannot write given data into single file.");
		    MEMFREE(isCellCentered);
		    fclose(fileHandle_p);
		    return;
		}
		isCellCentered[i_data] = TRUE;
		break;
	    default:
		sprintf(error_msg, "saveDX: I don't know how to handle function space type %d",
			getFunctionSpaceType(data_pp[i_data]));
		Dudley_setError(TYPE_ERROR, error_msg);
		MEMFREE(isCellCentered);
		fclose(fileHandle_p);
		return;
	    }
	}
    }
    /* select number of points and the mesh component */
    numPoints = mesh_p->Nodes->reducedNodesMapping->numTargets;
    nDim = mesh_p->Nodes->numDim;
    if (elementtype == DUDLEY_UNKNOWN)
	elementtype = DUDLEY_ELEMENTS;
    switch (elementtype)
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
	Dudley_setError(SYSTEM_ERROR, "saveDX: undefined element file");
	MEMFREE(isCellCentered);
	fclose(fileHandle_p);
	return;
    }

    /* map dudley element type to DX element type */
    TypeId = elements->etype;
    numDXNodesPerElement = 0;
    numCells = elements->numElements;
    if (TypeId == Dudley_Line2)
    {
	numDXNodesPerElement = 2;
	resortIndex = resort[0];
	strcpy(elemTypeStr, "lines");
    }
    else if (TypeId == Dudley_Tri3)
    {
	numDXNodesPerElement = 3;
	resortIndex = resort[1];
	strcpy(elemTypeStr, "triangles");
/*   } else if (TypeId==Rec4 || TypeId==Rec8 || TypeId==Rec9 || TypeId==Rec12 || TypeId==Rec16 || TypeId==Rec9Macro ) {
     numDXNodesPerElement = 4;
     resortIndex=resort[3];
     strcpy(elemTypeStr, "quads");*/
    }
    else if (TypeId == Dudley_Tet4)
    {
	numDXNodesPerElement = 4;
	resortIndex = resort[2];
	strcpy(elemTypeStr, "tetrahedra");
/*   } else if (TypeId==Hex8 || TypeId==Hex20 || TypeId==Hex32 ) {
     numDXNodesPerElement = 8;
     resortIndex=resort[4];
     strcpy(elemTypeStr, "cubes");*/
    }
    else
    {
	sprintf(error_msg, "saveDX: Element type %s is not supported by DX", elements->ename);
	Dudley_setError(VALUE_ERROR, error_msg);
	MEMFREE(isCellCentered);
	fclose(fileHandle_p);
	return;
    }

    /* positions */
    fprintf(fileHandle_p, "object 1 class array type float rank 1 shape %d items %d data follows\n", nDim, numPoints);
    for (i = 0; i < numPoints; i++)
    {
	p = mesh_p->Nodes->reducedNodesMapping->map[i];
	for (j = 0; j < nDim; j++)
	    fprintf(fileHandle_p, " %g", mesh_p->Nodes->Coordinates[INDEX2(j, p, nDim)]);
	fprintf(fileHandle_p, "\n");
    }
    /* connection table */
    NN = elements->numNodes;
    fprintf(fileHandle_p, "object 2 class array type int rank 1 shape %d items %d data follows\n", numDXNodesPerElement,
	    numCells);
    for (i = 0; i < numCells; i++)
    {
	for (j = 0; j < numDXNodesPerElement; j++)
	    fprintf(fileHandle_p, " %d",
		    mesh_p->Nodes->reducedNodesMapping->target[elements->Nodes[INDEX2(resortIndex[j], i, NN)]]);
	fprintf(fileHandle_p, "\n");
    }
    fprintf(fileHandle_p, "attribute \"element type\" string \"%s\"\n", elemTypeStr);
    fprintf(fileHandle_p, "attribute \"ref\" string \"positions\"\n");

    /* data */
    object_count = 2;
    for (i_data = 0; i_data < num_data; ++i_data)
    {
	if (!isEmpty(data_pp[i_data]))
	{
	    object_count++;
	    rank = getDataPointRank(data_pp[i_data]);
	    nComp = getDataPointSize(data_pp[i_data]);
	    fprintf(fileHandle_p, "object %d class array type float rank %d ", object_count, rank);
	    if (0 < rank)
	    {
		fprintf(fileHandle_p, "shape ");
		for (i = 0; i < rank; i++)
		    fprintf(fileHandle_p, "%d ", getDataPointShape(data_pp[i_data], i));
	    }
	    if (isCellCentered[i_data])
	    {
		if (Dudley_Assemble_reducedIntegrationOrder(data_pp[i_data]))
		{
		    numPointsPerSample = (elements->numLocalDim == 0) ? 0 : 1;
		}
		else
		{
		    numPointsPerSample = (elements->numLocalDim == 0) ? 0 : (elements->numLocalDim + 1);
		}
		if (numPointsPerSample > 0)
		{
		    fprintf(fileHandle_p, "items %d data follows\n", numCells);
		    for (i = 0; i < elements->numElements; i++)
		    {
			values = getSampleDataRO(data_pp[i_data], i);
			for (k = 0; k < nComp; k++)
			{
			    if (isExpanded(data_pp[i_data]))
			    {
				rtmp = 0.;
				for (j = 0; j < numPointsPerSample; j++)
				    rtmp += values[INDEX2(k, j, nComp)];
				fprintf(fileHandle_p, " %g", rtmp / numPointsPerSample);
			    }
			    else
			    {
				fprintf(fileHandle_p, " %g", values[k]);
			    }
			}
			fprintf(fileHandle_p, "\n");
		    }
		    fprintf(fileHandle_p, "attribute \"dep\" string \"connections\"\n");
		}
	    }
	    else
	    {
		fprintf(fileHandle_p, "items %d data follows\n", numPoints);
		for (i = 0; i < numPoints; i++)
		{
		    values = getSampleDataRO(data_pp[i_data], i);
		    for (k = 0; k < nComp; k++)
			fprintf(fileHandle_p, " %g", values[k]);
		    fprintf(fileHandle_p, "\n");
		}
		fprintf(fileHandle_p, "attribute \"dep\" string \"positions\"\n");
	    }
	}
    }

    /* and finish it up */
    if (num_data == 0)
    {
	fprintf(fileHandle_p, "object %d class field\n", object_count + 1);
	fprintf(fileHandle_p, "component \"positions\" value 1\n");
	fprintf(fileHandle_p, "component \"connections\" value 2\n");
    }
    else
    {
	object_count = 2;
	for (i_data = 0; i_data < num_data; ++i_data)
	{
	    if (!isEmpty(data_pp[i_data]))
	    {
		object_count++;
		fprintf(fileHandle_p, "object \"%s\" class field\n", names_p[i_data]);
		fprintf(fileHandle_p, "component \"positions\" value 1\n");
		fprintf(fileHandle_p, "component \"connections\" value 2\n");
		fprintf(fileHandle_p, "component \"data\" value %d\n", object_count);
	    }
	}
    }
    fprintf(fileHandle_p, "end\n");
    /* close the file */
    fclose(fileHandle_p);
    MEMFREE(isCellCentered);
    return;
}
