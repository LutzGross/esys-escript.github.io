
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

/*   Dudley: Mesh: prepares the mesh for further calculations  */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

void Dudley_Mesh_prepare(Dudley_Mesh * in, bool_t optimize)
{
    dim_t newGlobalNumDOFs = 0, numReducedNodes = 0, i;
    index_t *distribution = NULL, *maskReducedNodes = NULL, *indexReducedNodes = NULL, *node_distribution = NULL;
    if (in == NULL)
	return;
    if (in->Nodes == NULL)
	return;

    Dudley_Mesh_setOrders(in);

    /* first step is to distribute the elements according to a global distribution of DOF */

    distribution = TMPMEMALLOC(in->MPIInfo->size + 1, index_t);
    node_distribution = TMPMEMALLOC(in->MPIInfo->size + 1, index_t);
    if (!(Dudley_checkPtr(distribution) || Dudley_checkPtr(node_distribution)))
    {
	/* first we create dense labeling for the DOFs */

	newGlobalNumDOFs = Dudley_NodeFile_createDenseDOFLabeling(in->Nodes);

	/* create a distribution of the global DOFs and determine
	   the MPI_rank controlling the DOFs on this processor      */
	Esys_MPIInfo_setDistribution(in->MPIInfo, 0, newGlobalNumDOFs - 1, distribution);

	/* now the mesh is re-distributed according to the mpiRankOfDOF vector */
	/* this will redistribute the Nodes and Elements including overlap and will create an element coloring 
	   but will not create any mappings (see later in this function)                                   */
	if (Dudley_noError())
	    Dudley_Mesh_distributeByRankOfDOF(in, distribution);
    }

    /* at this stage we are able to start an optimization of the DOF distribution using ParMetis */
    /* on return distribution is altered and new DOF ids have been assigned */
    if (Dudley_noError() && optimize && in->MPIInfo->size > 1)
    {
	Dudley_Mesh_optimizeDOFDistribution(in, distribution);
	if (Dudley_noError())
	    Dudley_Mesh_distributeByRankOfDOF(in, distribution);
    }
    /* the local labeling of the degrees of free is optimized */
    if (Dudley_noError() && optimize)
    {
	Dudley_Mesh_optimizeDOFLabeling(in, distribution);
    }
    /* rearrange elements with the attempt to bring elements closer to memory locations of the nodes (distributed shared memory!): */
    if (Dudley_noError())
	Dudley_Mesh_optimizeElementOrdering(in);

    /* create the global indices */
    if (Dudley_noError())
    {

	maskReducedNodes = TMPMEMALLOC(in->Nodes->numNodes, index_t);
	indexReducedNodes = TMPMEMALLOC(in->Nodes->numNodes, index_t);
	if (!(Dudley_checkPtr(maskReducedNodes) || Dudley_checkPtr(indexReducedNodes)))
	{

/* useful DEBUG:
{index_t MIN_id,MAX_id;
printf("Mesh_prepare: global DOF : %d\n",newGlobalNumDOFs);
Dudley_NodeFile_setGlobalIdRange(&MIN_id,&MAX_id,in->Nodes);
printf("Mesh_prepare: global node id range = %d :%d\n", MIN_id,MAX_id);
Dudley_NodeFile_setIdRange(&MIN_id,&MAX_id,in->Nodes);
printf("Mesh_prepare: local node id range = %d :%d\n", MIN_id,MAX_id);
}
*/
#pragma omp parallel for private(i) schedule(static)
	    for (i = 0; i < in->Nodes->numNodes; ++i)
		maskReducedNodes[i] = -1;

	    Dudley_Mesh_markNodes(maskReducedNodes, 0, in, TRUE);

	    numReducedNodes = Dudley_Util_packMask(in->Nodes->numNodes, maskReducedNodes, indexReducedNodes);

	    Dudley_NodeFile_createDenseNodeLabeling(in->Nodes, node_distribution, distribution);
	    Dudley_NodeFile_createDenseReducedDOFLabeling(in->Nodes, maskReducedNodes);
	    Dudley_NodeFile_createDenseReducedNodeLabeling(in->Nodes, maskReducedNodes);
	    /* create the missing mappings */

	    if (Dudley_noError())
		Dudley_Mesh_createNodeFileMappings(in, numReducedNodes, indexReducedNodes, distribution,
						   node_distribution);
	}

	TMPMEMFREE(maskReducedNodes);
	TMPMEMFREE(indexReducedNodes);
    }

    TMPMEMFREE(distribution);
    TMPMEMFREE(node_distribution);

    Dudley_Mesh_setTagsInUse(in);
    return;
}

/*                                                      */
/*  tries to reduce the coloring for all element files: */
/*                                                      */
void Dudley_Mesh_createColoring(Dudley_Mesh * in, index_t * node_localDOF_map)
{
    if (Dudley_noError())
	Dudley_ElementFile_createColoring(in->Elements, in->Nodes->numNodes, node_localDOF_map);
    if (Dudley_noError())
	Dudley_ElementFile_createColoring(in->FaceElements, in->Nodes->numNodes, node_localDOF_map);
    if (Dudley_noError())
	Dudley_ElementFile_createColoring(in->Points, in->Nodes->numNodes, node_localDOF_map);
}

/*                                                                    */
/*  redistribute elements to minimize communication during assemblage */
/*                                                                    */
void Dudley_Mesh_optimizeElementOrdering(Dudley_Mesh * in)
{
    if (Dudley_noError())
	Dudley_ElementFile_optimizeOrdering(&(in->Elements));
    if (Dudley_noError())
	Dudley_ElementFile_optimizeOrdering(&(in->FaceElements));
    if (Dudley_noError())
	Dudley_ElementFile_optimizeOrdering(&(in->Points));
}

/*                                                                    */
/*  redistribute elements to minimize communication during assemblage */
void Dudley_Mesh_setTagsInUse(Dudley_Mesh * in)
{
    if (Dudley_noError())
	Dudley_NodeFile_setTagsInUse(in->Nodes);
    if (Dudley_noError())
	Dudley_ElementFile_setTagsInUse(in->Elements);
    if (Dudley_noError())
	Dudley_ElementFile_setTagsInUse(in->FaceElements);
    if (Dudley_noError())
	Dudley_ElementFile_setTagsInUse(in->Points);
}
