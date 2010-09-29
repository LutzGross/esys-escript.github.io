
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

/*   Dudley: Mesh: this will redistribute the Nodes and Elements including overlap */
/*   according to the dof_distribiution. It will create an element coloring but will not create any mappings. */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_distributeByRankOfDOF(Dudley_Mesh * self, index_t * dof_distribution)
{

    index_t min_dof_id, max_dof_id, *tmp_node_localDOF_map = NULL, *tmp_node_localDOF_mask = NULL;
    Paso_MPI_rank *mpiRankOfDOF = NULL;
    register index_t k;
    dim_t len, n, numDOFs;

    if (self == NULL)
	return;
    mpiRankOfDOF = TMPMEMALLOC(self->Nodes->numNodes, Paso_MPI_rank);
    if (!Dudley_checkPtr(mpiRankOfDOF))
    {

	Dudley_NodeFile_assignMPIRankToDOFs(self->Nodes, mpiRankOfDOF, dof_distribution);

	/* first the elements are redistributed according to mpiRankOfDOF */
	/* at the input the Node tables refering to a the local labeling of the nodes */
	/* while at the output they refer to the global labeling which is rectified in the next step */
	if (Dudley_noError())
	    Dudley_ElementFile_distributeByRankOfDOF(self->Elements, mpiRankOfDOF, self->Nodes->Id);
	if (Dudley_noError())
	    Dudley_ElementFile_distributeByRankOfDOF(self->FaceElements, mpiRankOfDOF, self->Nodes->Id);
	if (Dudley_noError())
	    Dudley_ElementFile_distributeByRankOfDOF(self->Points, mpiRankOfDOF, self->Nodes->Id);

	/* resolve the node ids */
	if (Dudley_noError())
	    Dudley_Mesh_resolveNodeIds(self);

	/* create a local labeling of the DOFs */
	Dudley_NodeFile_setDOFRange(&min_dof_id, &max_dof_id, self->Nodes);
	len = max_dof_id - min_dof_id + 1;
	tmp_node_localDOF_mask = TMPMEMALLOC(len, index_t);	/* local mask for used nodes */
	tmp_node_localDOF_map = TMPMEMALLOC(self->Nodes->numNodes, index_t);
	if (!((Dudley_checkPtr(tmp_node_localDOF_mask) && Dudley_checkPtr(tmp_node_localDOF_map))))
	{

#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < len; n++)
		tmp_node_localDOF_mask[n] = -1;

#pragma omp parallel for private (n) schedule(static)
	    for (n = 0; n < self->Nodes->numNodes; n++)
		tmp_node_localDOF_map[n] = -1;

#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < self->Nodes->numNodes; n++)
	    {
#ifdef BOUNDS_CHECK
		if ((self->Nodes->globalDegreesOfFreedom[n] - min_dof_id) >= len
		    || (self->Nodes->globalDegreesOfFreedom[n] - min_dof_id) < 0)
		{
		    printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
		    exit(1);
		}
#endif
		tmp_node_localDOF_mask[self->Nodes->globalDegreesOfFreedom[n] - min_dof_id] = n;
	    }

	    numDOFs = 0;
	    for (n = 0; n < len; n++)
	    {
		k = tmp_node_localDOF_mask[n];
		if (k >= 0)
		{
		    tmp_node_localDOF_mask[n] = numDOFs;
		    numDOFs++;
		}
	    }
#pragma omp parallel for private (n,k)
	    for (n = 0; n < self->Nodes->numNodes; n++)
	    {
		k = tmp_node_localDOF_mask[self->Nodes->globalDegreesOfFreedom[n] - min_dof_id];
		tmp_node_localDOF_map[n] = k;
	    }
	    /* create element coloring */
	    if (Dudley_noError())
		Dudley_Mesh_createColoring(self, tmp_node_localDOF_map);

	}
	TMPMEMFREE(tmp_node_localDOF_mask);
	TMPMEMFREE(tmp_node_localDOF_map);
    }
    TMPMEMFREE(mpiRankOfDOF);
    return;
}
