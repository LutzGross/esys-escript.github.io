
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

/*   Dudley: NodeFile : creates the mappings using the indexReducedNodes */
/*                 no distribution is happening                          */

/************************************************************************************/

#include "Mesh.h"
#define UNUSED -1

/************************************************************************************/

void Dudley_Mesh_createDOFMappingAndCoupling(Dudley_Mesh * in, bool_t use_reduced_elements)
{
    index_t min_DOF, max_DOF, *shared = NULL, *offsetInShared = NULL, *locDOFMask =
	NULL, i, k, myFirstDOF, myLastDOF, *nodeMask = NULL, firstDOF, lastDOF, *globalDOFIndex, *wanted_DOFs = NULL;
    dim_t mpiSize, len_loc_dof, numNeighbors, n, lastn, numNodes, *rcv_len = NULL, *snd_len = NULL, count;
    Esys_MPI_rank myRank, p, p_min, p_max, *neighbor = NULL;
    Paso_SharedComponents *rcv_shcomp = NULL, *snd_shcomp = NULL;
    Dudley_NodeMapping *this_mapping = NULL;
    Paso_Connector *this_connector = NULL;
    Paso_Distribution *dof_distribution;
    Esys_MPIInfo *mpi_info = in->MPIInfo;
#ifdef ESYS_MPI
    MPI_Request *mpi_requests = NULL;
    MPI_Status *mpi_stati = NULL;
#else
    int *mpi_requests = NULL, *mpi_stati = NULL;
#endif

    numNodes = in->Nodes->numNodes;
    if (use_reduced_elements)
    {
	dof_distribution = in->Nodes->reducedDegreesOfFreedomDistribution;
	globalDOFIndex = in->Nodes->globalReducedDOFIndex;
    }
    else
    {
	dof_distribution = in->Nodes->degreesOfFreedomDistribution;
	globalDOFIndex = in->Nodes->globalDegreesOfFreedom;
    }
    myFirstDOF = Paso_Distribution_getFirstComponent(dof_distribution);
    myLastDOF = Paso_Distribution_getLastComponent(dof_distribution);

    mpiSize = mpi_info->size;
    myRank = mpi_info->rank;

    min_DOF = Dudley_Util_getFlaggedMinInt(1, numNodes, globalDOFIndex, -1);
    max_DOF = Dudley_Util_getFlaggedMaxInt(1, numNodes, globalDOFIndex, -1);

    if (max_DOF < min_DOF)
    {
	min_DOF = myFirstDOF;
	max_DOF = myLastDOF - 1;
    }

    p_min = mpiSize;
    p_max = -1;
    if (max_DOF >= min_DOF)
    {
	for (p = 0; p < mpiSize; ++p)
	{
	    if (dof_distribution->first_component[p] <= min_DOF)
		p_min = p;
	    if (dof_distribution->first_component[p] <= max_DOF)
		p_max = p;
	}
    }

    len_loc_dof = max_DOF - min_DOF + 1;
    if (!((min_DOF <= myFirstDOF) && (myLastDOF - 1 <= max_DOF)))
    {
	Dudley_setError(SYSTEM_ERROR, "Local elements do not span local degrees of freedom.");
	return;
    }
    rcv_len = new  dim_t[mpiSize];
    snd_len = new  dim_t[mpiSize];
#ifdef ESYS_MPI
    mpi_requests = MEMALLOC(mpiSize * 2, MPI_Request);
    mpi_stati = MEMALLOC(mpiSize * 2, MPI_Status);
#else
    mpi_requests = MEMALLOC(mpiSize * 2, int);
    mpi_stati = MEMALLOC(mpiSize * 2, int);
#endif
    wanted_DOFs = new  index_t[numNodes];
    nodeMask = new  index_t[numNodes];
    neighbor = new  Esys_MPI_rank[mpiSize];
    shared = new  index_t[numNodes * (p_max - p_min + 1)];
    offsetInShared = new  index_t[mpiSize + 1];
    locDOFMask = new  index_t[len_loc_dof];
    if (!
	(Dudley_checkPtr(neighbor) || Dudley_checkPtr(shared) || Dudley_checkPtr(offsetInShared)
	 || Dudley_checkPtr(locDOFMask) || Dudley_checkPtr(nodeMask) || Dudley_checkPtr(rcv_len)
	 || Dudley_checkPtr(snd_len) || Dudley_checkPtr(mpi_requests) || Dudley_checkPtr(mpi_stati)
	 || Dudley_checkPtr(mpi_stati)))
    {

	memset(rcv_len, 0, sizeof(dim_t) * mpiSize);
#pragma omp parallel
	{
#pragma omp for private(i) schedule(static)
	    for (i = 0; i < len_loc_dof; ++i)
		locDOFMask[i] = UNUSED;
#pragma omp for private(i) schedule(static)
	    for (i = 0; i < numNodes; ++i)
		nodeMask[i] = UNUSED;
#pragma omp for private(i,k) schedule(static)
	    for (i = 0; i < numNodes; ++i)
	    {
		k = globalDOFIndex[i];
		if (k > -1)
		{
		    locDOFMask[k - min_DOF] = UNUSED - 1;
#ifdef BOUNDS_CHECK
		    if ((k - min_DOF) >= len_loc_dof)
		    {
			printf("BOUNDS_CHECK %s %d i=%d k=%d min_DOF=%d\n", __FILE__, __LINE__, i, k, min_DOF);
			exit(1);
		    }
#endif
		}
	    }

#pragma omp for private(i) schedule(static)
	    for (i = myFirstDOF - min_DOF; i < myLastDOF - min_DOF; ++i)
	    {
		locDOFMask[i] = i - myFirstDOF + min_DOF;
#ifdef BOUNDS_CHECK
		if (i < 0 || i >= len_loc_dof)
		{
		    printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i);
		    exit(1);
		}
#endif
	    }
	}

	numNeighbors = 0;
	n = 0;
	lastn = n;
	for (p = p_min; p <= p_max; ++p)
	{
	    firstDOF = MAX(min_DOF, dof_distribution->first_component[p]);
	    lastDOF = MIN(max_DOF + 1, dof_distribution->first_component[p + 1]);
	    if (p != myRank)
	    {
		for (i = firstDOF - min_DOF; i < lastDOF - min_DOF; ++i)
		{
#ifdef BOUNDS_CHECK
		    if (i < 0 || i >= len_loc_dof)
		    {
			printf("BOUNDS_CHECK %s %d p=%d i=%d\n", __FILE__, __LINE__, p, i);
			exit(1);
		    }
#endif
		    if (locDOFMask[i] == UNUSED - 1)
		    {
			locDOFMask[i] = myLastDOF - myFirstDOF + n;
			wanted_DOFs[n] = i + min_DOF;
			++n;
		    }
		}
		if (n > lastn)
		{
		    rcv_len[p] = n - lastn;
		    neighbor[numNeighbors] = p;
#ifdef BOUNDS_CHECK
		    if (numNeighbors < 0 || numNeighbors >= mpiSize + 1)
		    {
			printf("BOUNDS_CHECK %s %d p=%d numNeighbors=%d n=%d\n", __FILE__, __LINE__, p, numNeighbors,
			       n);
			exit(1);
		    }
#endif
		    offsetInShared[numNeighbors] = lastn;
		    numNeighbors++;
		    lastn = n;
		}
	    }
	}
#ifdef BOUNDS_CHECK
	if (numNeighbors < 0 || numNeighbors >= mpiSize + 1)
	{
	    printf("BOUNDS_CHECK %s %d numNeighbors=%d\n", __FILE__, __LINE__, numNeighbors);
	    exit(1);
	}
#endif
	offsetInShared[numNeighbors] = lastn;

	/* assign new DOF labels to nodes */
#pragma omp parallel for private(i,k) schedule(static)
	for (i = 0; i < numNodes; ++i)
	{
	    k = globalDOFIndex[i];
	    if (k > -1)
		nodeMask[i] = locDOFMask[k - min_DOF];
	}

	/* now we can set the mapping from nodes to local DOFs */
	this_mapping = Dudley_NodeMapping_alloc(numNodes, nodeMask, UNUSED);
	/* define how to get DOF values for controlled bu other processors */
#ifdef BOUNDS_CHECK
	for (i = 0; i < offsetInShared[numNeighbors]; ++i)
	{
	    if (i < 0 || i >= numNodes * (p_max - p_min + 1))
	    {
		printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i);
		exit(1);
	    }
	}
#endif
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < offsetInShared[numNeighbors]; ++i)
	    shared[i] = myLastDOF - myFirstDOF + i;

	rcv_shcomp =
	    Paso_SharedComponents_alloc(myLastDOF - myFirstDOF, numNeighbors, neighbor, shared, offsetInShared, 1, 0,
					mpi_info);

	/*
	 *    now we build the sender
	 */
#ifdef ESYS_MPI
	MPI_Alltoall(rcv_len, 1, MPI_INT, snd_len, 1, MPI_INT, mpi_info->comm);
#else
	for (p = 0; p < mpiSize; ++p)
	    snd_len[p] = rcv_len[p];
#endif
	count = 0;
	for (p = 0; p < rcv_shcomp->numNeighbors; p++)
	{
#ifdef ESYS_MPI
	    MPI_Isend(&(wanted_DOFs[rcv_shcomp->offsetInShared[p]]),
		      rcv_shcomp->offsetInShared[p + 1] - rcv_shcomp->offsetInShared[p], MPI_INT,
		      rcv_shcomp->neighbor[p], mpi_info->msg_tag_counter + myRank, mpi_info->comm,
		      &mpi_requests[count]);
#endif
	    count++;
	}
	n = 0;
	numNeighbors = 0;
	for (p = 0; p < mpiSize; p++)
	{
	    if (snd_len[p] > 0)
	    {
#ifdef ESYS_MPI
		MPI_Irecv(&(shared[n]), snd_len[p],
			  MPI_INT, p, mpi_info->msg_tag_counter + p, mpi_info->comm, &mpi_requests[count]);
#endif
		count++;
		neighbor[numNeighbors] = p;
		offsetInShared[numNeighbors] = n;
		numNeighbors++;
		n += snd_len[p];
	    }
	}
	mpi_info->msg_tag_counter += mpi_info->size;
	offsetInShared[numNeighbors] = n;
#ifdef ESYS_MPI
	MPI_Waitall(count, mpi_requests, mpi_stati);
#endif
	/* map global ids to local id's */
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < offsetInShared[numNeighbors]; ++i)
	{
	    shared[i] = locDOFMask[shared[i] - min_DOF];
	}

	snd_shcomp =
	    Paso_SharedComponents_alloc(myLastDOF - myFirstDOF, numNeighbors, neighbor, shared, offsetInShared, 1, 0,
					dof_distribution->mpi_info);

	if (Dudley_noError())
	    this_connector = Paso_Connector_alloc(snd_shcomp, rcv_shcomp);
	/* assign new DOF labels to nodes */
	Paso_SharedComponents_free(rcv_shcomp);
	Paso_SharedComponents_free(snd_shcomp);
    }
    delete[] rcv_len;
    delete[] snd_len;
    delete[] mpi_requests;
    delete[] mpi_stati;
    delete[] wanted_DOFs;
    delete[] nodeMask;
    delete[] neighbor;
    delete[] shared;
    delete[] offsetInShared;
    delete[] locDOFMask;
    if (Dudley_noError())
    {
	if (use_reduced_elements)
	{
	    in->Nodes->reducedDegreesOfFreedomMapping = this_mapping;
	    in->Nodes->reducedDegreesOfFreedomConnector = this_connector;
	}
	else
	{
	    in->Nodes->degreesOfFreedomMapping = this_mapping;
	    in->Nodes->degreesOfFreedomConnector = this_connector;
	}
    }
    else
    {
	Dudley_NodeMapping_free(this_mapping);
	Paso_Connector_free(this_connector);

    }
}

void Dudley_Mesh_createMappings(Dudley_Mesh * mesh, index_t * dof_distribution, index_t * node_distribution)
{
    int i;
    index_t *maskReducedNodes = NULL, *indexReducedNodes = NULL;
    dim_t numReducedNodes;

    maskReducedNodes = new  index_t[mesh->Nodes->numNodes];
    indexReducedNodes = new  index_t[mesh->Nodes->numNodes];

    if (!(Dudley_checkPtr(maskReducedNodes) || Dudley_checkPtr(indexReducedNodes)))
    {
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < mesh->Nodes->numNodes; ++i)
	    maskReducedNodes[i] = -1;
	Dudley_Mesh_markNodes(maskReducedNodes, 0, mesh, TRUE);

	numReducedNodes = Dudley_Util_packMask(mesh->Nodes->numNodes, maskReducedNodes, indexReducedNodes);
	if (Dudley_noError())
	    Dudley_Mesh_createNodeFileMappings(mesh, numReducedNodes, indexReducedNodes, dof_distribution,
					       node_distribution);
    }

    delete[] maskReducedNodes;
    delete[] indexReducedNodes;
}

void Dudley_Mesh_createNodeFileMappings(Dudley_Mesh * in, dim_t numReducedNodes, index_t * indexReducedNodes,
					index_t * dof_first_component, index_t * nodes_first_component)
{

    index_t myFirstDOF, myLastDOF, myFirstNode, myLastNode, *reduced_dof_first_component = NULL, *nodeMask = NULL,
	*reduced_nodes_first_component = NULL, k, *maskMyReducedDOF = NULL, *indexMyReducedDOF =
	NULL, *maskMyReducedNodes = NULL, *indexMyReducedNodes = NULL;
    dim_t myNumDOF, myNumNodes, myNumReducedNodes, myNumReducedDOF, globalNumReducedNodes, globalNumReducedDOF, i,
	mpiSize;
    Esys_MPI_rank myRank;

    mpiSize = in->Nodes->MPIInfo->size;
    myRank = in->Nodes->MPIInfo->rank;

    /* mark the nodes used by the reduced mesh */

    reduced_dof_first_component = new  index_t[mpiSize + 1];
    reduced_nodes_first_component = new  index_t[mpiSize + 1];

    if (!(Dudley_checkPtr(reduced_dof_first_component) || Dudley_checkPtr(reduced_nodes_first_component)))
    {

	myFirstDOF = dof_first_component[myRank];
	myLastDOF = dof_first_component[myRank + 1];
	myNumDOF = myLastDOF - myFirstDOF;

	myFirstNode = nodes_first_component[myRank];
	myLastNode = nodes_first_component[myRank + 1];
	myNumNodes = myLastNode - myFirstNode;

	maskMyReducedDOF = new  index_t[myNumDOF];
	indexMyReducedDOF = new  index_t[myNumDOF];
	maskMyReducedNodes = new  index_t[myNumNodes];
	indexMyReducedNodes = new  index_t[myNumNodes];

	if (!
	    (Dudley_checkPtr(maskMyReducedDOF) || Dudley_checkPtr(indexMyReducedDOF)
	     || Dudley_checkPtr(maskMyReducedNodes) || Dudley_checkPtr(indexMyReducedNodes)))
	{

#pragma omp parallel private(i)
	    {
#pragma omp for schedule(static)
		for (i = 0; i < myNumNodes; ++i)
		    maskMyReducedNodes[i] = -1;
#pragma omp for schedule(static)
		for (i = 0; i < myNumDOF; ++i)
		    maskMyReducedDOF[i] = -1;
#pragma omp for private(k) schedule(static)
		for (i = 0; i < numReducedNodes; ++i)
		{
		    k = in->Nodes->globalNodesIndex[indexReducedNodes[i]];
		    if ((k >= myFirstNode) && (myLastNode > k))
			maskMyReducedNodes[k - myFirstNode] = i;
		    k = in->Nodes->globalDegreesOfFreedom[indexReducedNodes[i]];
		    if ((k >= myFirstDOF) && (myLastDOF > k))
		    {
			maskMyReducedDOF[k - myFirstDOF] = i;
		    }
		}
	    }
	    myNumReducedNodes = Dudley_Util_packMask(myNumNodes, maskMyReducedNodes, indexMyReducedNodes);
	    myNumReducedDOF = Dudley_Util_packMask(myNumDOF, maskMyReducedDOF, indexMyReducedDOF);

#ifdef ESYS_MPI
	    MPI_Allgather(&myNumReducedNodes, 1, MPI_INT, reduced_nodes_first_component, 1, MPI_INT,
			  in->Nodes->MPIInfo->comm);
	    MPI_Allgather(&myNumReducedDOF, 1, MPI_INT, reduced_dof_first_component, 1, MPI_INT,
			  in->Nodes->MPIInfo->comm);
#else
	    reduced_nodes_first_component[0] = myNumReducedNodes;
	    reduced_dof_first_component[0] = myNumReducedDOF;
#endif
	    globalNumReducedNodes = 0;
	    globalNumReducedDOF = 0;
	    for (i = 0; i < mpiSize; ++i)
	    {
		k = reduced_nodes_first_component[i];
		reduced_nodes_first_component[i] = globalNumReducedNodes;
		globalNumReducedNodes += k;

		k = reduced_dof_first_component[i];
		reduced_dof_first_component[i] = globalNumReducedDOF;
		globalNumReducedDOF += k;
	    }
	    reduced_nodes_first_component[mpiSize] = globalNumReducedNodes;
	    reduced_dof_first_component[mpiSize] = globalNumReducedDOF;
	    /* ==== distribution of Nodes =============================== */
	    in->Nodes->nodesDistribution = Paso_Distribution_alloc(in->Nodes->MPIInfo, nodes_first_component, 1, 0);

	    /* ==== distribution of DOFs =============================== */
	    in->Nodes->degreesOfFreedomDistribution =
		Paso_Distribution_alloc(in->Nodes->MPIInfo, dof_first_component, 1, 0);

	    /* ==== distribution of reduced Nodes =============================== */
	    in->Nodes->reducedNodesDistribution =
		Paso_Distribution_alloc(in->Nodes->MPIInfo, reduced_nodes_first_component, 1, 0);

	    /* ==== distribution of reduced DOF =============================== */
	    in->Nodes->reducedDegreesOfFreedomDistribution =
		Paso_Distribution_alloc(in->Nodes->MPIInfo, reduced_dof_first_component, 1, 0);
	}
	delete[] maskMyReducedDOF;
	delete[] indexMyReducedDOF;
	delete[] maskMyReducedNodes;
	delete[] indexMyReducedNodes;
    }
    delete[] reduced_dof_first_component;
    delete[] reduced_nodes_first_component;

    nodeMask = new  index_t[in->Nodes->numNodes];
    if (!Dudley_checkPtr(nodeMask) && Dudley_noError())
    {

	/* ==== nodes mapping which is a dummy structure ======== */
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < in->Nodes->numNodes; ++i)
	    nodeMask[i] = i;
	in->Nodes->nodesMapping = Dudley_NodeMapping_alloc(in->Nodes->numNodes, nodeMask, UNUSED);

	/* ==== mapping between nodes and reduced nodes ========== */
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < in->Nodes->numNodes; ++i)
	    nodeMask[i] = UNUSED;
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < numReducedNodes; ++i)
	    nodeMask[indexReducedNodes[i]] = i;
	in->Nodes->reducedNodesMapping = Dudley_NodeMapping_alloc(in->Nodes->numNodes, nodeMask, UNUSED);
    }
    delete[] nodeMask;
    /* ==== mapping between nodes and DOFs + DOF connector ========== */
    if (Dudley_noError())
        Dudley_Mesh_createDOFMappingAndCoupling(in, FALSE);
    /* ==== mapping between nodes and reduced DOFs + reduced DOF connector ========== */
    if (Dudley_noError())
        Dudley_Mesh_createDOFMappingAndCoupling(in, TRUE);

    /* get the Ids for DOFs and reduced nodes */
    if (Dudley_noError())
    {
#pragma omp parallel private(i)
	{
#pragma omp for
	    for (i = 0; i < in->Nodes->reducedNodesMapping->numTargets; ++i)
		in->Nodes->reducedNodesId[i] = in->Nodes->Id[in->Nodes->reducedNodesMapping->map[i]];
#pragma omp for
	    for (i = 0; i < in->Nodes->degreesOfFreedomMapping->numTargets; ++i)
		in->Nodes->degreesOfFreedomId[i] = in->Nodes->Id[in->Nodes->degreesOfFreedomMapping->map[i]];
#pragma omp for
	    for (i = 0; i < in->Nodes->reducedDegreesOfFreedomMapping->numTargets; ++i)
		in->Nodes->reducedDegreesOfFreedomId[i] =
		    in->Nodes->Id[in->Nodes->reducedDegreesOfFreedomMapping->map[i]];
	}
    }
    else
    {
	Dudley_NodeMapping_free(in->Nodes->nodesMapping);
	Dudley_NodeMapping_free(in->Nodes->reducedNodesMapping);
	Dudley_NodeMapping_free(in->Nodes->degreesOfFreedomMapping);
	Dudley_NodeMapping_free(in->Nodes->reducedDegreesOfFreedomMapping);
	Paso_Distribution_free(in->Nodes->nodesDistribution);
	Paso_Distribution_free(in->Nodes->reducedNodesDistribution);
	Paso_Distribution_free(in->Nodes->degreesOfFreedomDistribution);
	Paso_Distribution_free(in->Nodes->reducedDegreesOfFreedomDistribution);
	Paso_Connector_free(in->Nodes->degreesOfFreedomConnector);
	Paso_Connector_free(in->Nodes->reducedDegreesOfFreedomConnector);
	in->Nodes->nodesMapping = NULL;
	in->Nodes->reducedNodesMapping = NULL;
	in->Nodes->degreesOfFreedomMapping = NULL;
	in->Nodes->reducedDegreesOfFreedomMapping = NULL;
	in->Nodes->nodesDistribution = NULL;
	in->Nodes->reducedNodesDistribution = NULL;
	in->Nodes->degreesOfFreedomDistribution = NULL;
	in->Nodes->reducedDegreesOfFreedomDistribution = NULL;
	in->Nodes->degreesOfFreedomConnector = NULL;
	in->Nodes->reducedDegreesOfFreedomConnector = NULL;
    }
}
