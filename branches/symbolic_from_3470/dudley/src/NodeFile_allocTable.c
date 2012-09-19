
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

/**************************************************************/

/*   Dudley: Mesh: NodeFile */

/**************************************************************/

#include "NodeFile.h"
#include "Util.h"

/**************************************************************/

/*  allocates the node table within an node file to hold numNodes of nodes. The LinearTo mapping, if it exists, */
/*  is frees. use Dudley_Mesh_setLinearMesh to create a new one. */

void Dudley_NodeFile_allocTable(Dudley_NodeFile * in, dim_t numNodes)
{
    index_t *Id2 = NULL, *Tag2 = NULL, *globalDegreesOfFreedom2 = NULL, *globalReducedDOFIndex2 = NULL,
	*globalReducedNodesIndex2 = NULL, *globalNodesIndex2 = NULL, *reducedNodesId2 = NULL, *degreesOfFreedomId2 =
	NULL, *reducedDegreesOfFreedomId2 = NULL;
    double *Coordinates2 = NULL;
    dim_t n, i;

    /*  allocate memory: */
    Id2 = MEMALLOC(numNodes, index_t);
    Coordinates2 = MEMALLOC(numNodes * in->numDim, double);
    Tag2 = MEMALLOC(numNodes, index_t);
    globalDegreesOfFreedom2 = MEMALLOC(numNodes, index_t);
    globalReducedDOFIndex2 = MEMALLOC(numNodes, index_t);
    globalReducedNodesIndex2 = MEMALLOC(numNodes, index_t);
    globalNodesIndex2 = MEMALLOC(numNodes, index_t);
    reducedNodesId2 = MEMALLOC(numNodes, index_t);
    degreesOfFreedomId2 = MEMALLOC(numNodes, index_t);
    reducedDegreesOfFreedomId2 = MEMALLOC(numNodes, index_t);

    /*  if fine, freeate the old table and replace by new: */
    if (Dudley_checkPtr(Id2) || Dudley_checkPtr(Coordinates2) || Dudley_checkPtr(Tag2)
	|| Dudley_checkPtr(globalDegreesOfFreedom2)
	|| Dudley_checkPtr(globalReducedDOFIndex2)
	|| Dudley_checkPtr(globalReducedNodesIndex2)
	|| Dudley_checkPtr(globalNodesIndex2)
	|| Dudley_checkPtr(reducedNodesId2) || Dudley_checkPtr(degreesOfFreedomId2))
    {
	MEMFREE(Id2);
	MEMFREE(Coordinates2);
	MEMFREE(Tag2);
	MEMFREE(globalDegreesOfFreedom2);
	MEMFREE(globalReducedDOFIndex2);
	MEMFREE(globalReducedNodesIndex2);
	MEMFREE(globalNodesIndex2);
	MEMFREE(reducedNodesId2);
	MEMFREE(degreesOfFreedomId2);
	MEMFREE(reducedDegreesOfFreedomId2);
    }
    else
    {
	Dudley_NodeFile_freeTable(in);
	in->Id = Id2;
	in->Coordinates = Coordinates2;
	in->globalDegreesOfFreedom = globalDegreesOfFreedom2;
	in->Tag = Tag2;
	in->globalReducedDOFIndex = globalReducedDOFIndex2;
	in->globalReducedNodesIndex = globalReducedNodesIndex2;
	in->globalNodesIndex = globalNodesIndex2;
	in->reducedNodesId = reducedNodesId2;
	in->degreesOfFreedomId = degreesOfFreedomId2;
	in->reducedDegreesOfFreedomId = reducedDegreesOfFreedomId2;
	in->numNodes = numNodes;
	/* this initialization makes sure that data are located on the right processor */
#pragma omp parallel for private(n,i) schedule(static)
	for (n = 0; n < numNodes; n++)
	{
	    in->Id[n] = -1;
	    for (i = 0; i < in->numDim; i++)
		in->Coordinates[INDEX2(i, n, in->numDim)] = 0.;
	    in->Tag[n] = -1;
	    in->globalDegreesOfFreedom[n] = -1;
	    in->globalReducedDOFIndex[n] = -1;
	    in->globalReducedNodesIndex[n] = -1;
	    in->globalNodesIndex[n] = -1;
	    in->reducedNodesId[n] = -1;
	    in->degreesOfFreedomId[n] = -1;
	    in->reducedDegreesOfFreedomId[n] = -1;
	}
    }
    return;
}

/*  frees the node table within an node file: */

void Dudley_NodeFile_freeTable(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	MEMFREE(in->Id);
	MEMFREE(in->Coordinates);
	MEMFREE(in->globalDegreesOfFreedom);
	MEMFREE(in->globalReducedDOFIndex);
	MEMFREE(in->globalReducedNodesIndex);
	MEMFREE(in->globalNodesIndex);
	MEMFREE(in->Tag);
	MEMFREE(in->reducedNodesId);
	MEMFREE(in->degreesOfFreedomId);
	MEMFREE(in->reducedDegreesOfFreedomId);
	MEMFREE(in->tagsInUse);
	in->numTagsInUse = 0;
	Dudley_NodeMapping_free(in->nodesMapping);
	in->nodesMapping = NULL;
	Dudley_NodeMapping_free(in->reducedNodesMapping);
	in->reducedNodesMapping = NULL;
	Dudley_NodeMapping_free(in->degreesOfFreedomMapping);
	in->degreesOfFreedomMapping = NULL;
	Dudley_NodeMapping_free(in->reducedDegreesOfFreedomMapping);
	in->reducedDegreesOfFreedomMapping = NULL;
	Paso_Distribution_free(in->nodesDistribution);
	in->nodesDistribution = NULL;
	Paso_Distribution_free(in->reducedNodesDistribution);
	in->nodesDistribution = NULL;
	Paso_Distribution_free(in->degreesOfFreedomDistribution);
	in->degreesOfFreedomDistribution = NULL;
	Paso_Distribution_free(in->reducedDegreesOfFreedomDistribution);
	in->reducedDegreesOfFreedomDistribution = NULL;
	Paso_Connector_free(in->degreesOfFreedomConnector);
	in->degreesOfFreedomConnector = NULL;
	Paso_Connector_free(in->reducedDegreesOfFreedomConnector);
	in->reducedDegreesOfFreedomConnector = NULL;

	in->numTagsInUse = 0;
	in->numNodes = 0;
    }
}

void Dudley_NodeFile_setTagsInUse(Dudley_NodeFile * in)
{
    index_t *tagsInUse = NULL;
    dim_t numTagsInUse;
    if (in != NULL)
    {
	Dudley_Util_setValuesInUse(in->Tag, in->numNodes, &numTagsInUse, &tagsInUse, in->MPIInfo);
	if (Dudley_noError())
	{
	    MEMFREE(in->tagsInUse);
	    in->tagsInUse = tagsInUse;
	    in->numTagsInUse = numTagsInUse;
	}
    }
}
