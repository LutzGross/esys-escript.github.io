
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

/*   Dudley: Mesh: NodeFile                                   */

/*   creates a dense labeling of the global degrees of freedom  */
/*   and returns the new number of  global degrees of freedom  */

/************************************************************************************/

#include "NodeFile.h"

/************************************************************************************/

dim_t Dudley_NodeFile_createDenseDOFLabeling(Dudley_NodeFile * in)
{
    index_t min_dof, max_dof, unset_dof = -1, set_dof = 1, dof_0, dof_1, *DOF_buffer = NULL, k;
    Esys_MPI_rank buffer_rank, *distribution = NULL;
    dim_t p, buffer_len, n, myDOFs, *offsets = NULL, *loc_offsets = NULL, new_numGlobalDOFs = 0, myNewDOFs;
    bool *set_new_DOF = NULL;
#ifdef ESYS_MPI
    Esys_MPI_rank dest, source;
    MPI_Status status;
#endif

    /* get the global range of node ids */
    Dudley_NodeFile_setGlobalDOFRange(&min_dof, &max_dof, in);

    distribution = new  index_t[in->MPIInfo->size + 1];
    offsets = new  dim_t[in->MPIInfo->size];
    loc_offsets = new  dim_t[in->MPIInfo->size];
    set_new_DOF = new  bool[in->numNodes];

    if (!
	(Dudley_checkPtr(distribution) || Dudley_checkPtr(offsets) || Dudley_checkPtr(loc_offsets)
	 || Dudley_checkPtr(set_new_DOF)))
    {
	/* distribute the range of node ids */
	buffer_len = Esys_MPIInfo_setDistribution(in->MPIInfo, min_dof, max_dof, distribution);
	myDOFs = distribution[in->MPIInfo->rank + 1] - distribution[in->MPIInfo->rank];
	/* allocate buffers */
	DOF_buffer = new  index_t[buffer_len];
	if (!Dudley_checkPtr(DOF_buffer))
	{
	    /* fill DOF_buffer by the unset_dof marker to check if nodes are defined */
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < buffer_len; n++)
		DOF_buffer[n] = unset_dof;

	    /* fill the buffer by sending portions around in a circle */
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		if (p > 0)
		{		/* the initial send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    in->MPIInfo->msg_tag_counter++;
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
		dof_0 = distribution[buffer_rank];
		dof_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    k = in->globalDegreesOfFreedom[n];
		    if ((dof_0 <= k) && (k < dof_1))
		    {
			DOF_buffer[k - dof_0] = set_dof;
		    }
		}
	    }
	    /* count the entries in the DOF_buffer */
	    /* TODO: OMP parallel */
	    myNewDOFs = 0;
	    for (n = 0; n < myDOFs; ++n)
	    {
		if (DOF_buffer[n] == set_dof)
		{
		    DOF_buffer[n] = myNewDOFs;
		    myNewDOFs++;
		}
	    }
	    memset(loc_offsets, 0, in->MPIInfo->size * sizeof(dim_t));
	    loc_offsets[in->MPIInfo->rank] = myNewDOFs;
#ifdef ESYS_MPI
	    MPI_Allreduce(loc_offsets, offsets, in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm);
	    new_numGlobalDOFs = 0;
	    for (n = 0; n < in->MPIInfo->size; ++n)
	    {
		loc_offsets[n] = new_numGlobalDOFs;
		new_numGlobalDOFs += offsets[n];
	    }
#else
	    new_numGlobalDOFs = loc_offsets[0];
	    loc_offsets[0] = 0;
#endif
#pragma omp parallel
	    {
#pragma omp for private(n) schedule(static)
		for (n = 0; n < myDOFs; ++n)
		    DOF_buffer[n] += loc_offsets[in->MPIInfo->rank];
		/* now entries are collected from the buffer again by sending the entries around in a circle */
#pragma omp for private(n) schedule(static)
		for (n = 0; n < in->numNodes; ++n)
		    set_new_DOF[n] = TRUE;
	    }
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		dof_0 = distribution[buffer_rank];
		dof_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    k = in->globalDegreesOfFreedom[n];
		    if (set_new_DOF[n] && (dof_0 <= k) && (k < dof_1))
		    {
			in->globalDegreesOfFreedom[n] = DOF_buffer[k - dof_0];
			set_new_DOF[n] = FALSE;
		    }
		}
		if (p < in->MPIInfo->size - 1)
		{		/* the last send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    ESYS_MPI_INC_COUNTER(*(in->MPIInfo),1);
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
	    }
	}
	delete[] DOF_buffer;
    }
    delete[] distribution;
    delete[] loc_offsets;
    delete[] offsets;
    delete[] set_new_DOF;
    return new_numGlobalDOFs;
}

void Dudley_NodeFile_assignMPIRankToDOFs(Dudley_NodeFile * in, Esys_MPI_rank * mpiRankOfDOF, index_t * distribution)
{
    index_t min_DOF, max_DOF, k;
    dim_t n;
    Esys_MPI_rank p, p_min = in->MPIInfo->size, p_max = -1;
    /* first we calculate the min and max dof on this processor to reduce costs for searching */
    Dudley_NodeFile_setDOFRange(&min_DOF, &max_DOF, in);

    for (p = 0; p < in->MPIInfo->size; ++p)
    {
	if (distribution[p] <= min_DOF)
	    p_min = p;
	if (distribution[p] <= max_DOF)
	    p_max = p;
    }
#pragma omp parallel for private(n,k,p) schedule(static)
    for (n = 0; n < in->numNodes; ++n)
    {
	k = in->globalDegreesOfFreedom[n];
	for (p = p_min; p <= p_max; ++p)
	{
	    if (k < distribution[p + 1])
	    {
		mpiRankOfDOF[n] = p;
		break;
	    }
	}
    }
}

dim_t Dudley_NodeFile_createDenseReducedDOFLabeling(Dudley_NodeFile * in, index_t * reducedNodeMask)
{
    index_t min_dof, max_dof, unset_dof = -1, set_dof = 1, dof_0, dof_1, *DOF_buffer = NULL, k;
    Esys_MPI_rank buffer_rank, *distribution = NULL;
    dim_t p, buffer_len, n, myDOFs, *offsets = NULL, *loc_offsets = NULL, globalNumReducedDOFs = 0, myNewDOFs;
#ifdef ESYS_MPI
    Esys_MPI_rank dest, source;
    MPI_Status status;
#endif

    /* get the global range of node ids */
    Dudley_NodeFile_setGlobalDOFRange(&min_dof, &max_dof, in);

    distribution = new  index_t[in->MPIInfo->size + 1];
    offsets = new  dim_t[in->MPIInfo->size];
    loc_offsets = new  dim_t[in->MPIInfo->size];

    if (!(Dudley_checkPtr(distribution) || Dudley_checkPtr(offsets) || Dudley_checkPtr(loc_offsets)))
    {
	/* distribute the range of node ids */
	buffer_len = Esys_MPIInfo_setDistribution(in->MPIInfo, min_dof, max_dof, distribution);
	myDOFs = distribution[in->MPIInfo->rank + 1] - distribution[in->MPIInfo->rank];
	/* allocate buffers */
	DOF_buffer = new  index_t[buffer_len];
	if (!Dudley_checkPtr(DOF_buffer))
	{
	    /* fill DOF_buffer by the unset_dof marker to check if nodes are defined */
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < buffer_len; n++)
		DOF_buffer[n] = unset_dof;

	    /* fill the buffer by sending portions around in a circle */
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		if (p > 0)
		{		/* the initial send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    in->MPIInfo->msg_tag_counter++;
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
		dof_0 = distribution[buffer_rank];
		dof_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    if (reducedNodeMask[n] > -1)
		    {
			k = in->globalDegreesOfFreedom[n];
			if ((dof_0 <= k) && (k < dof_1))
			{
			    DOF_buffer[k - dof_0] = set_dof;
			}
		    }
		}
	    }
	    /* count the entries in the DOF_buffer */
	    /* TODO: OMP parallel */
	    myNewDOFs = 0;
	    for (n = 0; n < myDOFs; ++n)
	    {
		if (DOF_buffer[n] == set_dof)
		{
		    DOF_buffer[n] = myNewDOFs;
		    myNewDOFs++;
		}
	    }
	    memset(loc_offsets, 0, in->MPIInfo->size * sizeof(dim_t));
	    loc_offsets[in->MPIInfo->rank] = myNewDOFs;
#ifdef ESYS_MPI
	    MPI_Allreduce(loc_offsets, offsets, in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm);
	    globalNumReducedDOFs = 0;
	    for (n = 0; n < in->MPIInfo->size; ++n)
	    {
		loc_offsets[n] = globalNumReducedDOFs;
		globalNumReducedDOFs += offsets[n];
	    }
#else
	    globalNumReducedDOFs = loc_offsets[0];
	    loc_offsets[0] = 0;
#endif
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < myDOFs; ++n)
		DOF_buffer[n] += loc_offsets[in->MPIInfo->rank];
	    /* now entries are collected from the buffer again by sending the entries around in a circle */
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < in->numNodes; ++n)
		in->globalReducedDOFIndex[n] = loc_offsets[0] - 1;
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		dof_0 = distribution[buffer_rank];
		dof_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    if (reducedNodeMask[n] > -1)
		    {
			k = in->globalDegreesOfFreedom[n];
			if ((dof_0 <= k) && (k < dof_1))
			    in->globalReducedDOFIndex[n] = DOF_buffer[k - dof_0];
		    }
		}
		if (p < in->MPIInfo->size - 1)
		{		/* the last send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    ESYS_MPI_INC_COUNTER(*(in->MPIInfo),1);
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
	    }
	}
	delete[] DOF_buffer;
    }
    delete[] distribution;
    delete[] loc_offsets;
    delete[] offsets;
    return globalNumReducedDOFs;
}

dim_t Dudley_NodeFile_createDenseNodeLabeling(Dudley_NodeFile * in, index_t * node_distribution,
					      const index_t * dof_distribution)
{
    index_t myFirstDOF, myLastDOF, max_id, min_id, loc_max_id, loc_min_id, dof, id, itmp, nodeID_0, nodeID_1, dof_0,
	dof_1, *Node_buffer = NULL;
    dim_t n, my_buffer_len, buffer_len, globalNumNodes = 0, myNewNumNodes;
    Esys_MPI_rank p, buffer_rank;
    const index_t unset_nodeID = -1, set_nodeID = 1;
    const dim_t header_len = 2;
#ifdef ESYS_MPI
    Esys_MPI_rank dest, source;
    MPI_Status status;
#endif
    Esys_MPI_rank myRank = in->MPIInfo->rank;

    /* find the range of node ids controlled by me */

    myFirstDOF = dof_distribution[myRank];
    myLastDOF = dof_distribution[myRank + 1];
    max_id = -INDEX_T_MAX;
    min_id = INDEX_T_MAX;
#pragma omp parallel private(loc_max_id,loc_min_id)
    {
	loc_max_id = max_id;
	loc_min_id = min_id;
#pragma omp for private(n,dof,id) schedule(static)
	for (n = 0; n < in->numNodes; n++)
	{
	    dof = in->globalDegreesOfFreedom[n];
	    id = in->Id[n];
	    if ((myFirstDOF <= dof) && (dof < myLastDOF))
	    {
		loc_max_id = MAX(loc_max_id, id);
		loc_min_id = MIN(loc_min_id, id);
	    }
	}
#pragma omp critical
	{
	    max_id = MAX(loc_max_id, max_id);
	    min_id = MIN(loc_min_id, min_id);
	}
    }
    /* allocate a buffer */
    my_buffer_len = max_id >= min_id ? max_id - min_id + 1 : 0;

#ifdef ESYS_MPI
    MPI_Allreduce(&my_buffer_len, &buffer_len, 1, MPI_INT, MPI_MAX, in->MPIInfo->comm);
#else
    buffer_len = my_buffer_len;
#endif

    Node_buffer = new  index_t[buffer_len + header_len];
    if (!Dudley_checkPtr(Node_buffer))
    {
	/* mark and count the nodes in use */
#pragma omp parallel
	{
#pragma omp for private(n) schedule(static)
	    for (n = 0; n < buffer_len + header_len; n++)
		Node_buffer[n] = unset_nodeID;
#pragma omp for private(n) schedule(static)
	    for (n = 0; n < in->numNodes; n++)
		in->globalNodesIndex[n] = -1;
#pragma omp for private(n,dof,id) schedule(static)
	    for (n = 0; n < in->numNodes; n++)
	    {
		dof = in->globalDegreesOfFreedom[n];
		id = in->Id[n];
		if ((myFirstDOF <= dof) && (dof < myLastDOF))
		    Node_buffer[id - min_id + header_len] = set_nodeID;
	    }
	}
	myNewNumNodes = 0;
	for (n = 0; n < my_buffer_len; n++)
	{
	    if (Node_buffer[header_len + n] == set_nodeID)
	    {
		Node_buffer[header_len + n] = myNewNumNodes;
		myNewNumNodes++;
	    }
	}
	/* make the local number of nodes globally available */
#ifdef ESYS_MPI
	MPI_Allgather(&myNewNumNodes, 1, MPI_INT, node_distribution, 1, MPI_INT, in->MPIInfo->comm);
#else
	node_distribution[0] = myNewNumNodes;
#endif

	globalNumNodes = 0;
	for (p = 0; p < in->MPIInfo->size; ++p)
	{
	    itmp = node_distribution[p];
	    node_distribution[p] = globalNumNodes;
	    globalNumNodes += itmp;
	}
	node_distribution[in->MPIInfo->size] = globalNumNodes;

	/* offset nodebuffer */
	itmp = node_distribution[in->MPIInfo->rank];
#pragma omp for private(n) schedule(static)
	for (n = 0; n < my_buffer_len; n++)
	    Node_buffer[n + header_len] += itmp;

	/* now we send this buffer around to assign global node index: */
#ifdef ESYS_MPI
	dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	Node_buffer[0] = min_id;
	Node_buffer[1] = max_id;
	buffer_rank = in->MPIInfo->rank;
	for (p = 0; p < in->MPIInfo->size; ++p)
	{
	    nodeID_0 = Node_buffer[0];
	    nodeID_1 = Node_buffer[1];
	    dof_0 = dof_distribution[buffer_rank];
	    dof_1 = dof_distribution[buffer_rank + 1];
	    if (nodeID_0 <= nodeID_1)
	    {
#pragma omp for private(n,dof,id) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    dof = in->globalDegreesOfFreedom[n];
		    id = in->Id[n] - nodeID_0;
		    if ((dof_0 <= dof) && (dof < dof_1) && (id >= 0) && (id <= nodeID_1 - nodeID_0))
			in->globalNodesIndex[n] = Node_buffer[id + header_len];
		}
	    }
	    if (p < in->MPIInfo->size - 1)
	    {			/* the last send can be skipped */
#ifdef ESYS_MPI
		MPI_Sendrecv_replace(Node_buffer, buffer_len + header_len, MPI_INT,
				     dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
				     in->MPIInfo->comm, &status);
#endif
		ESYS_MPI_INC_COUNTER(*(in->MPIInfo),1);
	    }
	    buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
	}
    }
    delete[] Node_buffer;
    return globalNumNodes;
}

dim_t Dudley_NodeFile_createDenseReducedNodeLabeling(Dudley_NodeFile * in, index_t * reducedNodeMask)
{
    index_t min_node, max_node, unset_node = -1, set_node = 1, node_0, node_1, *Nodes_buffer = NULL, k;
    Esys_MPI_rank buffer_rank, *distribution = NULL;
    dim_t p, buffer_len, n, myNodes, *offsets = NULL, *loc_offsets = NULL, globalNumReducedNodes = 0, myNewNodes;
#ifdef ESYS_MPI
    Esys_MPI_rank dest, source;
    MPI_Status status;
#endif

    /* get the global range of node ids */
    Dudley_NodeFile_setGlobalNodeIDIndexRange(&min_node, &max_node, in);

    distribution = new  index_t[in->MPIInfo->size + 1];
    offsets = new  dim_t[in->MPIInfo->size];
    loc_offsets = new  dim_t[in->MPIInfo->size];

    if (!(Dudley_checkPtr(distribution) || Dudley_checkPtr(offsets) || Dudley_checkPtr(loc_offsets)))
    {
	/* distribute the range of node ids */
	buffer_len = Esys_MPIInfo_setDistribution(in->MPIInfo, min_node, max_node, distribution);
	myNodes = distribution[in->MPIInfo->rank + 1] - distribution[in->MPIInfo->rank];
	/* allocate buffers */
	Nodes_buffer = new  index_t[buffer_len];
	if (!Dudley_checkPtr(Nodes_buffer))
	{
	    /* fill Nodes_buffer by the unset_node marker to check if nodes are defined */
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < buffer_len; n++)
		Nodes_buffer[n] = unset_node;

	    /* fill the buffer by sending portions around in a circle */
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		if (p > 0)
		{		/* the initial send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(Nodes_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    in->MPIInfo->msg_tag_counter++;
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
		node_0 = distribution[buffer_rank];
		node_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    if (reducedNodeMask[n] > -1)
		    {
			k = in->globalNodesIndex[n];
			if ((node_0 <= k) && (k < node_1))
			{
			    Nodes_buffer[k - node_0] = set_node;
			}
		    }
		}
	    }
	    /* count the entries in the Nodes_buffer */
	    /* TODO: OMP parallel */
	    myNewNodes = 0;
	    for (n = 0; n < myNodes; ++n)
	    {
		if (Nodes_buffer[n] == set_node)
		{
		    Nodes_buffer[n] = myNewNodes;
		    myNewNodes++;
		}
	    }
	    memset(loc_offsets, 0, in->MPIInfo->size * sizeof(dim_t));
	    loc_offsets[in->MPIInfo->rank] = myNewNodes;
#ifdef ESYS_MPI
	    MPI_Allreduce(loc_offsets, offsets, in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm);
	    globalNumReducedNodes = 0;
	    for (n = 0; n < in->MPIInfo->size; ++n)
	    {
		loc_offsets[n] = globalNumReducedNodes;
		globalNumReducedNodes += offsets[n];
	    }
#else
	    globalNumReducedNodes = loc_offsets[0];
	    loc_offsets[0] = 0;
#endif
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < myNodes; ++n)
		Nodes_buffer[n] += loc_offsets[in->MPIInfo->rank];
	    /* now entries are collected from the buffer again by sending the entries around in a circle */
#pragma omp parallel for private(n) schedule(static)
	    for (n = 0; n < in->numNodes; ++n)
		in->globalReducedNodesIndex[n] = loc_offsets[0] - 1;
#ifdef ESYS_MPI
	    dest = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
	    source = Esys_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
#endif
	    buffer_rank = in->MPIInfo->rank;
	    for (p = 0; p < in->MPIInfo->size; ++p)
	    {
		node_0 = distribution[buffer_rank];
		node_1 = distribution[buffer_rank + 1];
#pragma omp parallel for private(n,k) schedule(static)
		for (n = 0; n < in->numNodes; n++)
		{
		    if (reducedNodeMask[n] > -1)
		    {
			k = in->globalNodesIndex[n];
			if ((node_0 <= k) && (k < node_1))
			    in->globalReducedNodesIndex[n] = Nodes_buffer[k - node_0];
		    }
		}
		if (p < in->MPIInfo->size - 1)
		{		/* the last send can be skipped */
#ifdef ESYS_MPI
		    MPI_Sendrecv_replace(Nodes_buffer, buffer_len, MPI_INT,
					 dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
					 in->MPIInfo->comm, &status);
#endif
		    ESYS_MPI_INC_COUNTER(*(in->MPIInfo),1);
		}
		buffer_rank = Esys_MPIInfo_mod(in->MPIInfo->size, buffer_rank - 1);
	    }
	}
	delete[] Nodes_buffer;
    }
    delete[] distribution;
    delete[] loc_offsets;
    delete[] offsets;
    return globalNumReducedNodes;
}
