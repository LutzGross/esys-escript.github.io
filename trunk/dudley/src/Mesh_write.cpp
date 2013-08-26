
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

/*   Dudley: write Mesh */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

/*  writes the mesh to the external file fname using the Dudley file format: */

void Dudley_Mesh_write(Dudley_Mesh * in, char *fname)
{
    char error_msg[LenErrorMsg_MAX];
    FILE *f;
    int NN, i, j, numDim;
    Dudley_TagMap *tag_map = in->TagMap;

    if (in->MPIInfo->size > 1)
    {
	Dudley_setError(IO_ERROR, "Mesh_write: only single processor runs are supported.");
	return;

    }
    /* open file */
    f = fopen(fname, "w");
    if (f == NULL)
    {
	sprintf(error_msg, "Mesh_write: Opening file %s for writing failed.", fname);
	Dudley_setError(IO_ERROR, error_msg);
	return;
    }

    /* write header */

    fprintf(f, "%s\n", in->Name);

    /*  write nodes: */

    if (in->Nodes != NULL)
    {
	numDim = Dudley_Mesh_getDim(in);
	fprintf(f, "%1dD-Nodes %d\n", numDim, in->Nodes->numNodes);
	for (i = 0; i < in->Nodes->numNodes; i++)
	{
	    fprintf(f, "%d %d %d", in->Nodes->Id[i], in->Nodes->globalDegreesOfFreedom[i], in->Nodes->Tag[i]);
	    for (j = 0; j < numDim; j++)
		fprintf(f, " %20.15e", in->Nodes->Coordinates[INDEX2(j, i, numDim)]);
	    fprintf(f, "\n");
	}
    }
    else
    {
	fprintf(f, "0D-Nodes 0\n");
    }

    /*  write elements: */

    if (in->Elements != NULL)
    {
	fprintf(f, "%s %d\n", in->Elements->ename /*referenceElementSet->referenceElement->Type->Name */ ,
		in->Elements->numElements);
	NN = in->Elements->numNodes;
	for (i = 0; i < in->Elements->numElements; i++)
	{
	    fprintf(f, "%d %d", in->Elements->Id[i], in->Elements->Tag[i]);
	    for (j = 0; j < NN; j++)
		fprintf(f, " %d", in->Nodes->Id[in->Elements->Nodes[INDEX2(j, i, NN)]]);
	    fprintf(f, "\n");
	}
    }
    else
    {
	fprintf(f, "Tet4 0\n");
    }

    /*  write face elements: */
    if (in->FaceElements != NULL)
    {
	fprintf(f, "%s %d\n", in->FaceElements->ename /*referenceElementSet->referenceElement->Type->Name */ ,
		in->FaceElements->numElements);
	NN = in->FaceElements->numNodes;
	for (i = 0; i < in->FaceElements->numElements; i++)
	{
	    fprintf(f, "%d %d", in->FaceElements->Id[i], in->FaceElements->Tag[i]);
	    for (j = 0; j < NN; j++)
		fprintf(f, " %d", in->Nodes->Id[in->FaceElements->Nodes[INDEX2(j, i, NN)]]);
	    fprintf(f, "\n");
	}
    }
    else
    {
	fprintf(f, "Tri3 0\n");
    }

    /*  write points: */
    if (in->Points != NULL)
    {
	fprintf(f, "%s %d\n", in->Points->ename /*referenceElementSet->referenceElement->Type->Name */ ,
		in->Points->numElements);
	for (i = 0; i < in->Points->numElements; i++)
	{
	    fprintf(f, "%d %d %d\n", in->Points->Id[i], in->Points->Tag[i],
		    in->Nodes->Id[in->Points->Nodes[INDEX2(0, i, 1)]]);
	}
    }
    else
    {
	fprintf(f, "Point1 0\n");
    }

    /*  write tags: */
    if (tag_map)
    {
	fprintf(f, "Tags\n");
	while (tag_map)
	{
	    fprintf(f, "%s %d\n", tag_map->name, tag_map->tag_key);
	    tag_map = tag_map->next;
	}
    }
    fclose(f);
#ifdef Dudley_TRACE
    printf("mesh %s has been written to file %s\n", in->Name, fname);
#endif
}

void Dudley_PrintMesh_Info(Dudley_Mesh * in, bool full)
{
    int NN, i, j, numDim;
    Dudley_TagMap *tag_map = in->TagMap;

    fprintf(stdout, "Dudley_PrintMesh_Info running on CPU %d of %d\n", in->MPIInfo->rank, in->MPIInfo->size);
    fprintf(stdout, "\tMesh name '%s'\n", in->Name);
    fprintf(stdout, "\tApproximation order %d\n", in->approximationOrder);
    fprintf(stdout, "\tReduced Approximation order %d\n", in->reducedApproximationOrder);
    fprintf(stdout, "\tIntegration order %d\n", in->integrationOrder);
    fprintf(stdout, "\tReduced Integration order %d\n", in->reducedIntegrationOrder);

    /* write nodes: */
    if (in->Nodes != NULL)
    {
	numDim = Dudley_Mesh_getDim(in);
	fprintf(stdout, "\tNodes: %1dD-Nodes %d\n", numDim, in->Nodes->numNodes);
	if (full)
	{
	    fprintf(stdout, "\t     Id   Tag  gDOF   gNI grDfI  grNI:  Coordinates\n");
	    for (i = 0; i < in->Nodes->numNodes; i++)
	    {
		fprintf(stdout, "\t  %5d %5d %5d %5d %5d %5d: ", in->Nodes->Id[i], in->Nodes->Tag[i],
			in->Nodes->globalDegreesOfFreedom[i], in->Nodes->globalNodesIndex[i],
			in->Nodes->globalReducedDOFIndex[i], in->Nodes->globalReducedNodesIndex[i]);
		for (j = 0; j < numDim; j++)
		    fprintf(stdout, " %20.15e", in->Nodes->Coordinates[INDEX2(j, i, numDim)]);
		fprintf(stdout, "\n");
	    }
	}
    }
    else
    {
	fprintf(stdout, "\tNodes: 0D-Nodes 0\n");
    }

    /* write elements: */
    if (in->Elements != NULL)
    {
	int mine = 0, overlap = 0;
	for (i = 0; i < in->Elements->numElements; i++)
	{
	    if (in->Elements->Owner[i] == in->MPIInfo->rank)
		mine++;
	    else
		overlap++;
	}
	fprintf(stdout, "\tElements: %s %d (TypeId=%d) owner=%d overlap=%d\n",
		in->Elements->ename /*referenceElementSet->referenceElement->Type->Name */ , in->Elements->numElements,
		in->Elements->etype /*referenceElementSet->referenceElement->Type->TypeId */ , mine, overlap);
	NN = in->Elements->numNodes;
	if (full)
	{
	    fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
	    for (i = 0; i < in->Elements->numElements; i++)
	    {
		fprintf(stdout, "\t  %5d %5d %5d %5d: ", in->Elements->Id[i], in->Elements->Tag[i],
			in->Elements->Owner[i], in->Elements->Color[i]);
		for (j = 0; j < NN; j++)
		    fprintf(stdout, " %5d", in->Nodes->Id[in->Elements->Nodes[INDEX2(j, i, NN)]]);
		fprintf(stdout, "\n");
	    }
	}
    }
    else
    {
	fprintf(stdout, "\tElements: Tet4 0\n");
    }

    /* write face elements: */
    if (in->FaceElements != NULL)
    {
	int mine = 0, overlap = 0;
	for (i = 0; i < in->FaceElements->numElements; i++)
	{
	    if (in->FaceElements->Owner[i] == in->MPIInfo->rank)
		mine++;
	    else
		overlap++;
	}
	fprintf(stdout, "\tFace elements: %s %d (TypeId=%d) owner=%d overlap=%d\n",
		in->FaceElements->ename /*referenceElementSet->referenceElement->Type->Name */ ,
		in->FaceElements->numElements,
		in->FaceElements->etype /*->referenceElementSet->referenceElement->Type->TypeId*/ , mine, overlap);
	NN = in->FaceElements->numNodes;
	if (full)
	{
	    fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
	    for (i = 0; i < in->FaceElements->numElements; i++)
	    {
		fprintf(stdout, "\t  %5d %5d %5d %5d: ", in->FaceElements->Id[i], in->FaceElements->Tag[i],
			in->FaceElements->Owner[i], in->FaceElements->Color[i]);
		for (j = 0; j < NN; j++)
		    fprintf(stdout, " %5d", in->Nodes->Id[in->FaceElements->Nodes[INDEX2(j, i, NN)]]);
		fprintf(stdout, "\n");
	    }
	}
    }
    else
    {
	fprintf(stdout, "\tFace elements: Tri3 0\n");
    }

    /* write points: */
    if (in->Points != NULL)
    {
	int mine = 0, overlap = 0;
	for (i = 0; i < in->Points->numElements; i++)
	{
	    if (in->Points->Owner[i] == in->MPIInfo->rank)
		mine++;
	    else
		overlap++;
	}
	fprintf(stdout, "\tPoints: %s %d (TypeId=%d) owner=%d overlap=%d\n",
		in->Points->ename /*->referenceElementSet->referenceElement->Type->Name*/ , in->Points->numElements,
		in->Points->etype /*referenceElementSet->referenceElement->Type->TypeId */ , mine, overlap);
	if (full)
	{
	    fprintf(stdout, "\t     Id   Tag Owner Color:  Nodes\n");
	    for (i = 0; i < in->Points->numElements; i++)
	    {
		fprintf(stdout, "\t  %5d %5d %5d %5d %5d\n", in->Points->Id[i], in->Points->Tag[i],
			in->Points->Owner[i], in->Points->Color[i], in->Nodes->Id[in->Points->Nodes[INDEX2(0, i, 1)]]);
	    }
	}
    }
    else
    {
	fprintf(stdout, "\tPoints: Point1 0\n");
    }

    /* write tags: */
    if (tag_map)
    {
	fprintf(stdout, "\tTags:\n");
	while (tag_map)
	{
	    fprintf(stdout, "\t  %5d %s\n", tag_map->tag_key, tag_map->name);
	    tag_map = tag_map->next;
	}
    }
}
