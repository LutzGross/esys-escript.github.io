
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
/*                                                             */
/*   Dudley: Mesh : NodeFile */
/*                                                             */
/*   allocates and freeates node files                      */
/*                                                             */
/************************************************************************************/

#include "NodeFile.h"

/************************************************************************************/

/*   allocates a node file to hold nodes */
/*   use Dudley_NodeFile_allocTable to allocate the node table (Id,Coordinatess). */

Dudley_NodeFile *Dudley_NodeFile_alloc(dim_t numDim, Esys_MPIInfo * MPIInfo)
{
    Dudley_NodeFile *out;

    /*  allocate the return value */

    out = new Dudley_NodeFile;
    if (Dudley_checkPtr(out))
	return NULL;
    out->numNodes = 0;
    out->numDim = numDim;
    out->numTagsInUse = 0;
    out->Id = NULL;
    out->globalDegreesOfFreedom = NULL;
    out->Tag = NULL;
    out->Coordinates = NULL;
    out->status = DUDLEY_INITIAL_STATUS;

    out->nodesMapping = NULL;
    out->reducedNodesMapping = NULL;
    out->degreesOfFreedomMapping = NULL;
    out->reducedDegreesOfFreedomMapping = NULL;

    out->globalReducedDOFIndex = NULL;
    out->globalReducedNodesIndex = NULL;
    out->globalNodesIndex = NULL;
    out->reducedNodesId = NULL;
    out->degreesOfFreedomId = NULL;
    out->reducedDegreesOfFreedomId = NULL;
    out->nodesDistribution = NULL;
    out->reducedNodesDistribution = NULL;
    out->degreesOfFreedomDistribution = NULL;
    out->reducedDegreesOfFreedomDistribution = NULL;
    out->degreesOfFreedomConnector = NULL;
    out->reducedDegreesOfFreedomConnector = NULL;
    out->tagsInUse = NULL;

    out->MPIInfo = Esys_MPIInfo_getReference(MPIInfo);
    return out;
}

/*  frees a node file: */

void Dudley_NodeFile_free(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	Dudley_NodeFile_freeTable(in);
	Esys_MPIInfo_free(in->MPIInfo);
	delete in;
    }
}

index_t Dudley_NodeFile_getFirstReducedNode(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getFirstComponent(in->reducedNodesDistribution);
    }
    else
    {
	return 0;
    }
}

index_t Dudley_NodeFile_getLastReducedNode(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getLastComponent(in->reducedNodesDistribution);
    }
    else
    {
	return 0;
    }

}

dim_t Dudley_NodeFile_getGlobalNumReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getGlobalNumComponents(in->reducedNodesDistribution);
    }
    else
    {
	return 0;
    }

}

index_t *Dudley_NodeFile_borrowGlobalReducedNodesIndex(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->globalReducedNodesIndex;
    }
    else
    {
	return NULL;
    }
}

index_t Dudley_NodeFile_getFirstNode(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getFirstComponent(in->nodesDistribution);
    }
    else
    {
	return 0;
    }
}

index_t Dudley_NodeFile_getLastNode(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getLastComponent(in->nodesDistribution);
    }
    else
    {
	return 0;
    }

}

dim_t Dudley_NodeFile_getGlobalNumNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getGlobalNumComponents(in->nodesDistribution);
    }
    else
    {
	return 0;
    }

}

index_t *Dudley_NodeFile_borrowGlobalNodesIndex(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->globalNodesIndex;
    }
    else
    {
	return NULL;
    }
}

dim_t Dudley_NodeFile_getNumReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->reducedNodesMapping->numTargets;
    }
    else
    {
	return 0;
    }

}

dim_t Dudley_NodeFile_getNumDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getMyNumComponents(in->degreesOfFreedomDistribution);
    }
    else
    {
	return 0;
    }
}

dim_t Dudley_NodeFile_getNumNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->nodesMapping->numNodes;
    }
    else
    {
	return 0;
    }
}

dim_t Dudley_NodeFile_getNumReducedDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return Paso_Distribution_getMyNumComponents(in->reducedDegreesOfFreedomDistribution);
    }
    else
    {
	return 0;
    }
}

index_t *Dudley_NodeFile_borrowTargetReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->reducedNodesMapping->target;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowTargetDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->degreesOfFreedomMapping->target;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowTargetNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->nodesMapping->target;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowTargetReducedDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->reducedDegreesOfFreedomMapping->target;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowReducedNodesTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->reducedNodesMapping->map;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowDegreesOfFreedomTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->degreesOfFreedomMapping->map;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowNodesTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->nodesMapping->map;
    }
    else
    {
	return NULL;
    }
}

index_t *Dudley_NodeFile_borrowReducedDegreesOfFreedomTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
	return in->reducedDegreesOfFreedomMapping->map;
    }
    else
    {
	return NULL;
    }
}
