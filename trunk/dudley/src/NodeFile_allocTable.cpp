
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

/*   Dudley: Mesh: NodeFile */

/************************************************************************************/

#include "NodeFile.h"
#include "Util.h"

/************************************************************************************/

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
    Id2 = new  index_t[numNodes];
    Coordinates2 = new  double[numNodes * in->numDim];
    Tag2 = new  index_t[numNodes];
    globalDegreesOfFreedom2 = new  index_t[numNodes];
    globalReducedDOFIndex2 = new  index_t[numNodes];
    globalReducedNodesIndex2 = new  index_t[numNodes];
    globalNodesIndex2 = new  index_t[numNodes];
    reducedNodesId2 = new  index_t[numNodes];
    degreesOfFreedomId2 = new  index_t[numNodes];
    reducedDegreesOfFreedomId2 = new  index_t[numNodes];

    /*  if fine, free the old table and replace by new: */
    if (Dudley_checkPtr(Id2) || Dudley_checkPtr(Coordinates2) || Dudley_checkPtr(Tag2)
	|| Dudley_checkPtr(globalDegreesOfFreedom2)
	|| Dudley_checkPtr(globalReducedDOFIndex2)
	|| Dudley_checkPtr(globalReducedNodesIndex2)
	|| Dudley_checkPtr(globalNodesIndex2)
	|| Dudley_checkPtr(reducedNodesId2) || Dudley_checkPtr(degreesOfFreedomId2))
    {
	delete[] Id2;
	delete[] Coordinates2;
	delete[] Tag2;
	delete[] globalDegreesOfFreedom2;
	delete[] globalReducedDOFIndex2;
	delete[] globalReducedNodesIndex2;
	delete[] globalNodesIndex2;
	delete[] reducedNodesId2;
	delete[] degreesOfFreedomId2;
	delete[] reducedDegreesOfFreedomId2;
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
	delete[] in->Id;
	delete[] in->Coordinates;
	delete[] in->globalDegreesOfFreedom;
	delete[] in->globalReducedDOFIndex;
	delete[] in->globalReducedNodesIndex;
	delete[] in->globalNodesIndex;
	delete[] in->Tag;
	delete[] in->reducedNodesId;
	delete[] in->degreesOfFreedomId;
	delete[] in->reducedDegreesOfFreedomId;
	delete[] in->tagsInUse;
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
	    delete[] in->tagsInUse;
	    in->tagsInUse = tagsInUse;
	    in->numTagsInUse = numTagsInUse;
	}
    }
}
