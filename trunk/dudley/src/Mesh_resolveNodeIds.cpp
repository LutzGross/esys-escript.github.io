
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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

/*   Dudley: Mesh */

/*   at input the element nodes refers to the numbering defined the global Id assigned to the nodes in the */
/*   NodeFile. It is also not ensured that all nodes referred to by an element are actually available */
/*   on the process.  At the output, a local node labelling is used and all nodes are available */
/*   In particular the numbering of the element nodes is between 0 and in->NodeFile->numNodes */
/*   The function does not create a distribution of the degrees of freedom. */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Mesh.h"
#include "Util.h"

/************************************************************************************/

void Dudley_Mesh_resolveNodeIds(Dudley_Mesh * in)
{

    index_t min_id, max_id, min_id2, max_id2, global_min_id, global_max_id,
	*globalToNewLocalNodeLabels = NULL, *newLocalToGlobalNodeLabels = NULL;
    dim_t len, n, newNumNodes, numDim;
    Dudley_NodeFile *newNodeFile = NULL;
#ifdef ESYS_MPI
    index_t id_range[2], global_id_range[2];
#endif
    numDim = Dudley_Mesh_getDim(in);
    /*  find the minimum and maximum id used by elements: */
    min_id = INDEX_T_MAX;
    max_id = -INDEX_T_MAX;
    Dudley_ElementFile_setNodeRange(&min_id2, &max_id2, in->Elements);
    max_id = MAX(max_id, max_id2);
    min_id = MIN(min_id, min_id2);
    Dudley_ElementFile_setNodeRange(&min_id2, &max_id2, in->FaceElements);
    max_id = MAX(max_id, max_id2);
    min_id = MIN(min_id, min_id2);
    Dudley_ElementFile_setNodeRange(&min_id2, &max_id2, in->Points);
    max_id = MAX(max_id, max_id2);
    min_id = MIN(min_id, min_id2);
#ifdef ESYS_MPI
    id_range[0] = -min_id;
    id_range[1] = max_id;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    global_min_id = -global_id_range[0];
    global_max_id = global_id_range[1];
#else
    global_min_id = min_id;
    global_max_id = max_id;
#endif
#ifdef Dudley_TRACE
    printf("Node id range used by elements is %d:%d\n", global_min_id, global_max_id);
#else
    /* avoid unused var warning if Dudley_TRACE is not defined */
    (void)global_min_id;
    (void)global_max_id;
#endif
    if (min_id > max_id)
    {
	max_id = -1;
	min_id = 0;
    }

    /* allocate mappings for new local node labelling to global node labelling (newLocalToGlobalNodeLabels)
       and global node labelling to the new local node labelling (globalToNewLocalNodeLabels[i-min_id] is the 
       new local id of global node i) */
    len = (max_id >= min_id) ? max_id - min_id + 1 : 0;
    globalToNewLocalNodeLabels = new  index_t[len];	/* local mask for used nodes */
    newLocalToGlobalNodeLabels = new  index_t[len];
    if (!((Dudley_checkPtr(globalToNewLocalNodeLabels) && Dudley_checkPtr(newLocalToGlobalNodeLabels))))
    {

#pragma omp parallel
	{
#pragma omp for private(n) schedule(static)
	    for (n = 0; n < len; n++)
		newLocalToGlobalNodeLabels[n] = -1;
#pragma omp for private(n) schedule(static)
	    for (n = 0; n < len; n++)
		globalToNewLocalNodeLabels[n] = -1;
	}

	/*  mark the nodes referred by elements in globalToNewLocalNodeLabels which is currently used as a mask: */
	Dudley_Mesh_markNodes(globalToNewLocalNodeLabels, min_id, in, FALSE);

	/* create a local labelling newLocalToGlobalNodeLabels of the local nodes by packing the mask globalToNewLocalNodeLabels */

	newNumNodes = Dudley_Util_packMask(len, globalToNewLocalNodeLabels, newLocalToGlobalNodeLabels);

	/* invert the new labelling and shift the index newLocalToGlobalNodeLabels to global node ids */
#pragma omp parallel for private(n) schedule(static)
	for (n = 0; n < newNumNodes; n++)
	{
#ifdef BOUNDS_CHECK
	    if (n >= len || n < 0)
	    {
		printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n);
		exit(1);
	    }
	    if (newLocalToGlobalNodeLabels[n] >= len || newLocalToGlobalNodeLabels[n] < 0)
	    {
		printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n);
		exit(1);
	    }
#endif
	    globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]] = n;
	    newLocalToGlobalNodeLabels[n] += min_id;
	}
	/* create a new table */
	newNodeFile = Dudley_NodeFile_alloc(numDim, in->MPIInfo);
	if (Dudley_noError())
	{
	    Dudley_NodeFile_allocTable(newNodeFile, newNumNodes);
	}
	if (Dudley_noError())
	{
	    Dudley_NodeFile_gather_global(newLocalToGlobalNodeLabels, in->Nodes, newNodeFile);
	}
	if (Dudley_noError())
	{
	    Dudley_NodeFile_free(in->Nodes);
	    in->Nodes = newNodeFile;
	    /*  relabel nodes of the elements: */
	    Dudley_Mesh_relableElementNodes(globalToNewLocalNodeLabels, min_id, in);
	}
    }
    delete[] globalToNewLocalNodeLabels;
    delete[] newLocalToGlobalNodeLabels;
    if (!Dudley_noError())
    {
	Dudley_NodeFile_free(newNodeFile);
    }
}
