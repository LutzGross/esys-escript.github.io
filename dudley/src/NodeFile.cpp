
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************/
/*                                                             */
/*   Dudley: Mesh : NodeFile                                   */
/*                                                             */
/*   allocates and frees node files                            */
/*                                                             */
/****************************************************************************/

#include "NodeFile.h"

namespace dudley {

/*   allocates a node file to hold nodes */
/*   use Dudley_NodeFile_allocTable to allocate the node table (Id,Coordinates). */
Dudley_NodeFile *Dudley_NodeFile_alloc(dim_t numDim, escript::JMPI& MPIInfo)
{
    Dudley_NodeFile *out;

    out = new Dudley_NodeFile;
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
    out->tagsInUse = NULL;

    out->MPIInfo = MPIInfo;
    return out;
}

/*  frees a node file: */

void Dudley_NodeFile_free(Dudley_NodeFile * in)
{
    if (in != NULL)
    {
        Dudley_NodeFile_freeTable(in);
        delete in;
    }
}

index_t Dudley_NodeFile_getFirstReducedNode(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesDistribution->getFirstComponent();
    return 0;
}

index_t Dudley_NodeFile_getLastReducedNode(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesDistribution->getLastComponent();
    return 0;
}

dim_t Dudley_NodeFile_getGlobalNumReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesDistribution->getGlobalNumComponents();
    return 0;
}

index_t *Dudley_NodeFile_borrowGlobalReducedNodesIndex(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->globalReducedNodesIndex;
    return NULL;
}

index_t Dudley_NodeFile_getFirstNode(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesDistribution->getFirstComponent();
    return 0;
}

index_t Dudley_NodeFile_getLastNode(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesDistribution->getLastComponent();
    return 0;
}

dim_t Dudley_NodeFile_getGlobalNumNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesDistribution->getGlobalNumComponents();
    return 0;
}

index_t *Dudley_NodeFile_borrowGlobalNodesIndex(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->globalNodesIndex;
    return NULL;
}

dim_t Dudley_NodeFile_getNumReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesMapping->numTargets;
    return 0;
}

dim_t Dudley_NodeFile_getNumDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->degreesOfFreedomDistribution->getMyNumComponents();
    return 0;
}

dim_t Dudley_NodeFile_getNumNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesMapping->numNodes;
    return 0;
}

dim_t Dudley_NodeFile_getNumReducedDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedDegreesOfFreedomDistribution->getMyNumComponents();
    return 0;
}

index_t *Dudley_NodeFile_borrowTargetReducedNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesMapping->target;
    return NULL;
}

index_t *Dudley_NodeFile_borrowTargetDegreesOfFreedom(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->degreesOfFreedomMapping->target;
    return NULL;
}

index_t *Dudley_NodeFile_borrowTargetNodes(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesMapping->target;
    return NULL;
}

index_t *Dudley_NodeFile_borrowTargetReducedDegreesOfFreedom(Dudley_NodeFile* in)
{
    if (in != NULL)
        return in->reducedDegreesOfFreedomMapping->target;
    return NULL;
}

index_t *Dudley_NodeFile_borrowReducedNodesTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedNodesMapping->map;
    return NULL;
}

index_t *Dudley_NodeFile_borrowDegreesOfFreedomTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->degreesOfFreedomMapping->map;
    return NULL;
}

index_t *Dudley_NodeFile_borrowNodesTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->nodesMapping->map;
    return NULL;
}

index_t *Dudley_NodeFile_borrowReducedDegreesOfFreedomTarget(Dudley_NodeFile * in)
{
    if (in != NULL)
        return in->reducedDegreesOfFreedomMapping->map;
    return NULL;
}

} // namespace dudley

