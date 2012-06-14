
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

/************************************************************************/

/*   Dudley: Mesh: optimizes the labeling of the DOFs on each processor */

/************************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/**************************************************************/

void Dudley_Mesh_optimizeDOFLabeling(Dudley_Mesh * in, dim_t * distribution)
{

    index_t myFirstVertex, myLastVertex, *newGlobalDOFID = NULL, firstVertex, lastVertex;
    register index_t k;
    dim_t mpiSize, myNumVertices, len, p, i;
    Paso_Pattern *pattern = NULL;
    Esys_MPI_rank myRank, current_rank;
    Dudley_IndexList *index_list = NULL;
#ifdef ESYS_MPI
    Esys_MPI_rank dest, source;
    MPI_Status status;
#endif

    if (in == NULL)
	return;
    if (in->Nodes == NULL)
	return;

    myRank = in->MPIInfo->rank;
    mpiSize = in->MPIInfo->size;
    myFirstVertex = distribution[myRank];
    myLastVertex = distribution[myRank + 1];
    myNumVertices = myLastVertex - myFirstVertex;
    len = 0;
    for (p = 0; p < mpiSize; ++p)
	len = MAX(len, distribution[p + 1] - distribution[p]);

    index_list = TMPMEMALLOC(myNumVertices, Dudley_IndexList);
    newGlobalDOFID = TMPMEMALLOC(len, index_t);
    /* create the adjacency structure xadj and adjncy */
    if (!(Dudley_checkPtr(index_list) || Dudley_checkPtr(newGlobalDOFID)))
    {
#pragma omp parallel private(i)
	{
#pragma omp for schedule(static)
	    for (i = 0; i < myNumVertices; ++i)
	    {
		index_list[i].extension = NULL;
		index_list[i].n = 0;
	    }
	    /*  insert contributions from element matrices into colums index index_list: */
	    Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
								      in->Elements, in->Nodes->globalDegreesOfFreedom,
								      in->Nodes->globalDegreesOfFreedom);
	    Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
								      in->FaceElements,
								      in->Nodes->globalDegreesOfFreedom,
								      in->Nodes->globalDegreesOfFreedom);
	    Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
								      in->Points, in->Nodes->globalDegreesOfFreedom,
								      in->Nodes->globalDegreesOfFreedom);
	}
	/* create the local matrix pattern */
	pattern =
	    Dudley_IndexList_createPattern(0, myNumVertices, index_list, myFirstVertex, myLastVertex, -myFirstVertex);

	/* clean up index list */
	if (index_list != NULL)
	{
#pragma omp parallel for private(i)
	    for (i = 0; i < myNumVertices; ++i)
		Dudley_IndexList_free(index_list[i].extension);
	}

	if (Dudley_noError())
	    Paso_Pattern_reduceBandwidth(pattern, newGlobalDOFID);

	Paso_Pattern_free(pattern);
    }
    Esys_MPIInfo_noError(in->MPIInfo);
    if (Dudley_noError())
    {
	/* shift new labeling to create a global id */
#pragma omp parallel for private(i)
	for (i = 0; i < myNumVertices; ++i)
	    newGlobalDOFID[i] += myFirstVertex;

	/* distribute new labeling to other processors */
#ifdef ESYS_MPI
	dest = Esys_MPIInfo_mod(mpiSize, myRank + 1);
	source = Esys_MPIInfo_mod(mpiSize, myRank - 1);
#endif
	current_rank = myRank;
	for (p = 0; p < mpiSize; ++p)
	{
	    firstVertex = distribution[current_rank];
	    lastVertex = distribution[current_rank + 1];
#pragma omp parallel for private(i,k)
	    for (i = 0; i < in->Nodes->numNodes; ++i)
	    {
		k = in->Nodes->globalDegreesOfFreedom[i];
		if ((firstVertex <= k) && (k < lastVertex))
		{
		    in->Nodes->globalDegreesOfFreedom[i] = newGlobalDOFID[k - firstVertex];
		}
	    }

	    if (p < mpiSize - 1)
	    {			/* the final send can be skipped */
#ifdef ESYS_MPI
		MPI_Sendrecv_replace(newGlobalDOFID, len, MPI_INT,
				     dest, in->MPIInfo->msg_tag_counter,
				     source, in->MPIInfo->msg_tag_counter, in->MPIInfo->comm, &status);
#endif
		in->MPIInfo->msg_tag_counter++;
		current_rank = Esys_MPIInfo_mod(mpiSize, current_rank - 1);
	    }
	}
    }
    TMPMEMFREE(index_list);
    TMPMEMFREE(newGlobalDOFID);
#if 0
    for (i = 0; i < in->Nodes->numNodes; ++i)
	printf("%d ", in->Nodes->globalDegreesOfFreedom[i]);
    printf("\n");
#endif
    return;
}
